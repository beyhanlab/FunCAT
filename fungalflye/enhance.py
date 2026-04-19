"""
funcat.enhance
~~~~~~~~~~~~~~~~~~
Assembly enhancement modules:
  1. Adaptive Flye parameter selection
  2. Iterative Medaka polishing with convergence detection
  3. Purge Duplicates (diploid haplotig removal)
  4. Contig confidence scoring
"""

import subprocess
import shutil
import json
import math
from pathlib import Path
from Bio import SeqIO


# ------------------------------------------------
# helpers
# ------------------------------------------------

def run(cmd):
    print(f"\n[funcat] Running: {cmd}\n")
    subprocess.run(cmd, shell=True, check=True)


def run_capture(cmd):
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return result.stdout.strip()


# ================================================
# MODULE 1 — Adaptive Flye parameter selection
# ================================================

def suggest_flye_params(reads, genome_size_bp, threads=8):
    """
    Analyse reads and return optimised Flye parameters as a dict.
    Prints a human-readable explanation of every decision made.
    """

    print("\n" + "=" * 60)
    print("🔬 Module 1 — Adaptive parameter selection")
    print("=" * 60)

    # --- read stats via seqkit ---
    print("\n[funcat] Analysing read characteristics...")

    stats_raw = run_capture(
        f"seqkit fx2tab -n -l -g {reads} 2>/dev/null"
    )

    lengths = []
    gcs = []

    for line in stats_raw.splitlines():
        parts = line.strip().split("\t")
        if len(parts) >= 3:
            try:
                lengths.append(int(parts[1]))
                gcs.append(float(parts[2]))
            except ValueError:
                pass

    if not lengths:
        print("[funcat] ⚠️  Could not parse read stats — using defaults")
        return {}

    total_bases = sum(lengths)
    coverage = total_bases / genome_size_bp
    mean_gc = sum(gcs) / len(gcs) if gcs else 50.0

    sorted_len = sorted(lengths, reverse=True)
    cumsum, read_n50 = 0, 0
    for L in sorted_len:
        cumsum += L
        if cumsum >= total_bases / 2:
            read_n50 = L
            break

    print(f"\n  Read N50       : {read_n50:,} bp")
    print(f"  Est. coverage  : {coverage:.1f}x")
    print(f"  Mean GC        : {mean_gc:.1f}%")

    params = {}
    reasoning = []

    # --- min-overlap ---
    # Rule: min_overlap ≈ 0.5 × read_n50, capped between 1000 and 10000
    raw_overlap = int(read_n50 * 0.5)
    min_overlap = max(1000, min(raw_overlap, 10000))
    params["min_overlap"] = min_overlap
    reasoning.append(
        f"  --min-overlap {min_overlap}  "
        f"(0.5 × read N50 {read_n50:,}, capped 1k–10k)"
    )

    # --- asm-coverage ---
    # Rule: cap at actual coverage if lower than default 60x
    if coverage < 40:
        asm_cov = max(int(coverage * 0.8), 20)
        reasoning.append(
            f"  --asm-coverage {asm_cov}  "
            f"(low coverage {coverage:.0f}x — reduced to avoid read starvation)"
        )
    elif coverage > 120:
        asm_cov = 80
        reasoning.append(
            f"  --asm-coverage 80  "
            f"(very high coverage {coverage:.0f}x — capped to reduce compute)"
        )
    else:
        asm_cov = 60
        reasoning.append(
            f"  --asm-coverage 60  (coverage {coverage:.0f}x — standard)"
        )
    params["asm_coverage"] = asm_cov

    # --- iterations ---
    # Rule: more iterations for lower coverage (graph needs more resolution)
    if coverage < 30:
        iters = 4
        reasoning.append(
            "  --iterations 4  (low coverage — extra iterations for graph resolution)"
        )
    else:
        iters = 3
        reasoning.append("  --iterations 3  (standard)")
    params["iterations"] = iters

    # --- AT-rich warning ---
    if mean_gc < 35:
        params["at_rich"] = True
        reasoning.append(
            f"  ⚠️  AT-rich genome ({mean_gc:.1f}% GC) — "
            "consider --min-overlap reduction if assembly fragments"
        )

    print("\n[funcat] Recommended Flye parameters:")
    for r in reasoning:
        print(r)

    return params


# ================================================
# MODULE 2 — Iterative Medaka polishing
# ================================================

def _count_changes(before_fasta, after_fasta, threads=4):
    """
    Count sequence-level changes between two assemblies using minimap2 asm5.
    Works with Medaka 2.x which no longer outputs VCF files.

    Aligns before → after, sums mismatches and indels from the cs tag.
    Returns total number of corrected bases (substitutions + indels).
    Falls back to 0 if minimap2 is unavailable or alignment fails.
    """
    import subprocess
    import re

    try:
        result = subprocess.run(
            [
                "minimap2", "-x", "asm5",
                "--cs", "-t", str(threads),
                str(before_fasta), str(after_fasta),
            ],
            capture_output=True, text=True, timeout=300,
        )
        changes = 0
        for line in result.stdout.splitlines():
            if line.startswith("S") or line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            # Parse cs tag — count mismatches (*) and indels (+ -)
            for col in cols:
                if col.startswith("cs:Z:"):
                    cs = col[5:]
                    # substitutions: *XY  insertions: +seq  deletions: -seq
                    changes += len(re.findall(r'\*[acgtACGT]{2}', cs))
                    changes += sum(len(m) - 1 for m in re.findall(r'\+[acgtACGT]+', cs))
                    changes += sum(len(m) - 1 for m in re.findall(r'-[acgtACGT]+', cs))
        return changes
    except Exception:
        return 0


def run_medaka_iterative(
    assembly,
    reads,
    outdir,
    threads,
    model,
    max_rounds=3,
    convergence_threshold=50,
):
    """
    Run Medaka up to max_rounds times.
    Stops early when variant calls drop below convergence_threshold
    (meaning the assembly has converged and further polishing won't help).
    Returns path to final polished FASTA.
    """

    print("\n" + "=" * 60)
    print("🔬 Module 2 — Iterative Medaka polishing")
    print(f"   Max rounds: {max_rounds}  |  Convergence threshold: {convergence_threshold} base corrections")
    print("=" * 60)

    current = Path(assembly)
    medaka_base = Path(outdir) / "medaka"
    medaka_base.mkdir(exist_ok=True)

    state_file = medaka_base / "polishing_state.json"

    # Resume support: load previous state
    if state_file.exists():
        state = json.loads(state_file.read_text())
        last_round = state.get("last_round", 0)
        last_variants = state.get("last_variants", None)
        current = Path(state.get("current_fasta", assembly))
        print(f"[funcat] Resuming from round {last_round}")
    else:
        state = {}
        last_round = 0
        last_variants = None

    final_fasta = None

    for round_num in range(last_round + 1, max_rounds + 1):

        round_dir = medaka_base / f"round_{round_num}"

        polished = round_dir / "consensus.fasta"

        if polished.exists():
            print(f"[funcat] Round {round_num} already done — skipping")
            current = polished
            final_fasta = polished
            continue

        print(f"\n[funcat] Medaka round {round_num} / {max_rounds}")

        round_dir.mkdir(exist_ok=True)

        # Snapshot the input so we can compare before/after for convergence
        current_before_polish = current

        run(
            f"medaka_consensus "
            f"-i {reads} "
            f"-d {current} "
            f"-o {round_dir} "
            f"-t {threads} "
            f"-m {model}"
        )

        if not polished.exists():
            raise RuntimeError(f"Medaka round {round_num} failed")

        # Count sequence changes between input and polished assembly.
        # Uses direct FASTA comparison via minimap2 asm5 — compatible with
        # Medaka 2.x which no longer outputs VCF files.
        changes_this_round = _count_changes(current_before_polish, polished, threads=threads)

        print(
            f"[funcat] Round {round_num} complete — "
            f"{changes_this_round} bases corrected"
        )

        # Save state for resumability
        state = {
            "last_round": round_num,
            "last_variants": changes_this_round,
            "current_fasta": str(polished),
        }
        state_file.write_text(json.dumps(state, indent=2))

        current = polished
        final_fasta = polished

        # Convergence check — stop if corrections fall below threshold
        if changes_this_round <= convergence_threshold:
            print(
                f"\n✅ Converged after round {round_num} "
                f"({changes_this_round} bases corrected ≤ threshold {convergence_threshold})"
            )
            break

        if last_variants is not None:
            improvement = last_variants - changes_this_round
            if improvement <= 0:
                print(
                    f"\n✅ No further improvement after round {round_num} — stopping"
                )
                break

        last_variants = changes_this_round

    if final_fasta is None:
        raise RuntimeError("Medaka polishing produced no output")

    print(f"\n[funcat] Final polished assembly: {final_fasta}")
    return final_fasta


# ================================================
# MODULE 3 — Purge Duplicates (diploid)
# ================================================

def run_purge_dups(assembly, reads, outdir, threads, minimap2_preset):
    """
    Run purge_dups to remove haplotig duplicates from diploid assemblies.
    Returns path to purged primary assembly.

    Requires: purge_dups, minimap2, split_fa, pbcstat, calcuts, get_seqs
    (all installed via: conda install -c bioconda purge_dups)
    """

    print("\n" + "=" * 60)
    print("🔬 Module 3 — Purge Duplicates (diploid haplotig removal)")
    print("=" * 60)

    purge_dir = Path(outdir) / "purge_dups"
    purge_dir.mkdir(exist_ok=True)

    purged_primary = purge_dir / "purged.fa"
    haplotigs = purge_dir / "hap.fa"

    if purged_primary.exists():
        print("[funcat] Existing purge_dups output detected — skipping")
        return purged_primary, haplotigs

    assembly = Path(assembly)

    # Check tools
    missing = [t for t in ["purge_dups", "split_fa", "pbcstat", "calcuts", "get_seqs"]
               if shutil.which(t) is None]
    if missing:
        print(f"\n⚠️  purge_dups tools not found: {missing}")
        print("   Install with: conda install -c bioconda purge_dups")
        print("   Skipping haplotig purging — assembly will be unpurged\n")
        return assembly, None

    print("\n[funcat] Step 1 — self-mapping assembly")
    self_paf = purge_dir / "self.paf"
    run(f"minimap2 -xasm5 -DP -t {threads} {assembly} {assembly} > {self_paf}")

    print("[funcat] Step 2 — mapping reads to assembly")
    read_paf = purge_dir / "reads.paf"
    run(f"minimap2 -x {minimap2_preset} -t {threads} {assembly} {reads} > {read_paf}")

    print("[funcat] Step 3 — coverage histogram")
    stat_file = purge_dir / "PB.stat"
    base_cov = purge_dir / "PB.base.cov"
    run(f"pbcstat {read_paf} -O {purge_dir}/")

    print("[funcat] Step 4 — calculating cutoffs")
    cutoffs_file = purge_dir / "cutoffs"
    run(f"calcuts {stat_file} > {cutoffs_file}")

    print("[funcat] Step 5 — splitting assembly")
    split_asm = purge_dir / "split.fa"
    run(f"split_fa {assembly} > {split_asm}")

    split_paf = purge_dir / "split.paf"
    run(f"minimap2 -xasm5 -DP -t {threads} {split_asm} {split_asm} > {split_paf}")

    print("[funcat] Step 6 — purging haplotigs")
    bed_file = purge_dir / "dups.bed"
    run(
        f"purge_dups -2 -T {cutoffs_file} -c {base_cov} "
        f"{self_paf} > {bed_file}"
    )

    run(f"get_seqs -e {bed_file} {assembly} -p {purge_dir}/purged")

    if not purged_primary.exists():
        print("⚠️  purge_dups did not produce output — using unpurged assembly")
        return assembly, None

    # Report
    orig_contigs = sum(1 for r in SeqIO.parse(str(assembly), "fasta"))
    purged_contigs = sum(1 for r in SeqIO.parse(str(purged_primary), "fasta"))
    removed = orig_contigs - purged_contigs

    print(f"\n✅ Purge Duplicates complete")
    print(f"   Original contigs : {orig_contigs}")
    print(f"   After purging    : {purged_contigs}  ({removed} haplotigs removed)")
    if haplotigs.exists():
        hap_contigs = sum(1 for r in SeqIO.parse(str(haplotigs), "fasta"))
        print(f"   Haplotigs saved  : {hap_contigs} → {haplotigs}")

    return purged_primary, haplotigs


# ================================================
# MODULE 4 — Contig confidence scoring
# ================================================

def score_contig_confidence(assembly, reads, outdir, threads, minimap2_preset):
    """
    Map reads back to assembly and compute per-contig metrics:
      - mean coverage
      - coverage uniformity (coefficient of variation)
      - a confidence label: GOOD / REVIEW / FLAG

    Writes a TSV report and prints a summary.
    Returns path to the TSV.
    """

    print("\n" + "=" * 60)
    print("🔬 Module 4 — Contig confidence scoring")
    print("=" * 60)

    score_dir = Path(outdir) / "confidence"
    score_dir.mkdir(exist_ok=True)

    bam = score_dir / "reads.bam"
    report = score_dir / "contig_confidence.tsv"

    if report.exists():
        print("[funcat] Existing confidence report detected — skipping")
        return report

    # Check samtools
    if shutil.which("samtools") is None:
        print("⚠️  samtools not found — skipping confidence scoring")
        print("   Install with: conda install -c bioconda samtools")
        return None

    print("[funcat] Mapping reads for coverage analysis...")

    run(
        f"minimap2 -ax {minimap2_preset} -t {threads} {assembly} {reads} "
        f"| samtools sort -@ {threads} -o {bam}"
    )
    run(f"samtools index {bam}")

    print("[funcat] Computing per-contig coverage statistics...")

    depth_raw = run_capture(f"samtools depth -a {bam}")

    # Build per-contig depth lists
    contig_depths = {}
    for line in depth_raw.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        ctg, depth = parts[0], int(parts[2])
        contig_depths.setdefault(ctg, []).append(depth)

    # Compute overall median coverage for relative thresholds
    all_depths = [d for depths in contig_depths.values() for d in depths]
    all_depths.sort()
    median_cov = all_depths[len(all_depths) // 2] if all_depths else 1

    rows = []
    flagged = 0
    review = 0

    for contig, depths in contig_depths.items():
        n = len(depths)
        mean_cov = sum(depths) / n
        variance = sum((d - mean_cov) ** 2 for d in depths) / n
        std_cov = math.sqrt(variance)
        cv = (std_cov / mean_cov) if mean_cov > 0 else 999
        low_cov_pct = 100 * sum(1 for d in depths if d < 5) / n

        # Confidence rules
        if mean_cov > median_cov * 1.8:
            label = "FLAG"       # likely collapsed repeat
            reason = "coverage >1.8x median (possible collapsed repeat)"
            flagged += 1
        elif mean_cov < median_cov * 0.2:
            label = "FLAG"       # possible contamination / misassembly
            reason = "coverage <0.2x median (possible contamination)"
            flagged += 1
        elif cv > 1.5 or low_cov_pct > 20:
            label = "REVIEW"     # uneven coverage
            reason = f"uneven coverage (CV={cv:.2f}, {low_cov_pct:.0f}% low-cov bases)"
            review += 1
        else:
            label = "GOOD"
            reason = "coverage looks uniform and expected"

        rows.append({
            "contig": contig,
            "length_bp": n,
            "mean_coverage": round(mean_cov, 1),
            "coverage_cv": round(cv, 3),
            "low_cov_pct": round(low_cov_pct, 1),
            "label": label,
            "reason": reason,
        })

    # Write TSV
    with open(report, "w") as f:
        headers = ["contig", "length_bp", "mean_coverage",
                   "coverage_cv", "low_cov_pct", "label", "reason"]
        f.write("\t".join(headers) + "\n")
        for row in rows:
            f.write("\t".join(str(row[h]) for h in headers) + "\n")

    total = len(rows)
    good = total - flagged - review

    print(f"\n✅ Confidence scoring complete ({total} contigs)")
    print(f"   GOOD   : {good}")
    print(f"   REVIEW : {review}  (check these before publishing)")
    print(f"   FLAG   : {flagged}  (likely collapsed repeats or contamination)")
    print(f"\n   Full report: {report}")

    if flagged > 0:
        print("\n   Flagged contigs:")
        for row in rows:
            if row["label"] == "FLAG":
                print(f"     {row['contig']} — {row['reason']}")

    return report


# ================================================
# MODULE 6 — Illumina polishing
# ================================================

def run_illumina_polishing(
    assembly,
    illumina_r1,
    illumina_r2,
    outdir,
    threads,
    polisher="polypolish"
):
    """
    Polish a long-read assembly using Illumina short reads.
    
    Uses either Polypolish (default, recommended) or Pilon for polishing.
    Significantly improves base accuracy and reduces internal stop codons.
    
    Returns path to polished assembly.
    """
    
    print("\n" + "=" * 60)
    print("🔬 Module 6 — Illumina polishing")
    print(f"   Polisher: {polisher}")
    print("=" * 60)
    
    polish_dir = Path(outdir) / "illumina_polish"
    polish_dir.mkdir(exist_ok=True)
    
    polished_assembly = polish_dir / f"polished_{polisher}.fasta"
    
    if polished_assembly.exists():
        print("[funcat] Existing Illumina polishing output detected — skipping")
        return polished_assembly
    
    assembly = Path(assembly)
    illumina_r1 = Path(illumina_r1)
    illumina_r2 = Path(illumina_r2)
    
    # Check required tools
    required_tools = ["bwa", "samtools"]
    if polisher == "polypolish":
        required_tools.extend(["polypolish", "polypolish_insert_filter"])
    elif polisher == "pilon":
        required_tools.append("pilon")
    
    missing_tools = [tool for tool in required_tools if shutil.which(tool) is None]
    
    if missing_tools:
        print(f"\n⚠️  Missing tools for Illumina polishing: {missing_tools}")
        if polisher == "polypolish":
            print("   Install with: conda install -c bioconda polypolish bwa samtools")
        else:
            print("   Install with: conda install -c bioconda pilon bwa samtools")
        print("   Skipping Illumina polishing — assembly will be unpolished\n")
        return assembly
    
    print("\n[funcat] Step 1 — Indexing assembly for BWA")
    run(f"bwa index {assembly}")
    
    print("[funcat] Step 2 — Mapping Illumina reads")
    sam_r1 = polish_dir / "alignments_1.sam"
    sam_r2 = polish_dir / "alignments_2.sam"
    
    run(f"bwa mem -t {threads} -a {assembly} {illumina_r1} > {sam_r1}")
    run(f"bwa mem -t {threads} -a {assembly} {illumina_r2} > {sam_r2}")
    
    if polisher == "polypolish":
        print("[funcat] Step 3 — Filtering alignments (Polypolish)")
        filtered_r1 = polish_dir / "filtered_1.sam"
        filtered_r2 = polish_dir / "filtered_2.sam"
        
        run(f"polypolish_insert_filter --in1 {sam_r1} --in2 {sam_r2} --out1 {filtered_r1} --out2 {filtered_r2}")
        
        print("[funcat] Step 4 — Polishing with Polypolish")
        run(f"polypolish {assembly} {filtered_r1} {filtered_r2} > {polished_assembly}")
        
    elif polisher == "pilon":
        print("[funcat] Step 3 — Converting to BAM and sorting")
        bam_file = polish_dir / "illumina_mapped.bam"
        
        run(f"samtools view -bS {sam_r1} | samtools sort -@ {threads} -o {bam_file}")
        run(f"samtools index {bam_file}")
        
        print("[funcat] Step 4 — Polishing with Pilon")
        run(f"pilon --genome {assembly} --frags {bam_file} --output polished_pilon --outdir {polish_dir} --changes --threads {threads}")
        
        # Pilon outputs with different name
        pilon_output = polish_dir / "polished_pilon.fasta"
        if pilon_output.exists():
            run(f"cp {pilon_output} {polished_assembly}")
    
    if not polished_assembly.exists():
        print(f"⚠️  {polisher} polishing failed — using original assembly")
        return assembly
    
    # Calculate improvement statistics
    try:
        original_size = sum(len(record.seq) for record in SeqIO.parse(str(assembly), "fasta"))
        polished_size = sum(len(record.seq) for record in SeqIO.parse(str(polished_assembly), "fasta"))
        size_change = polished_size - original_size
        
        print(f"\n✅ Illumina polishing complete")
        print(f"   Polisher used    : {polisher}")
        print(f"   Original size    : {original_size:,} bp")
        print(f"   Polished size    : {polished_size:,} bp")
        print(f"   Size change      : {size_change:+,} bp")
        print(f"   Output assembly  : {polished_assembly}")
        
        # Suggest running BUSCO to check improvement
        print(f"\n   💡 Tip: Run BUSCO on both assemblies to measure improvement:")
        print(f"      Internal stop codons should drop significantly (22% → <5%)")
        
    except Exception as e:
        print(f"✅ Illumina polishing complete (statistics calculation failed: {e})")
    
    return polished_assembly
