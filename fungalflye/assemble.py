import subprocess
from pathlib import Path
import shutil
import shutil as sh
import json
from Bio import SeqIO

from .enhance import (
    suggest_flye_params,
    run_medaka_iterative,
    run_purge_dups,
    score_contig_confidence,
    run_illumina_polishing,
)
from .scaffold import run_scaffold, run_telomere_scaffolding


# ------------------------------------------------
# Read type config
# ------------------------------------------------

READ_TYPE_CONFIGS = {
    "nano-hq": {
        "flye_flag": "--nano-hq",
        "minimap2_preset": "map-ont",
        "medaka_model": "r1041_e82_400bps_hac_g632",
        "label": "Nanopore HQ (R10.4+, Q20)",
    },
    "nano-raw": {
        "flye_flag": "--nano-raw",
        "minimap2_preset": "map-ont",
        "medaka_model": "r941_min_hac_g507",
        "label": "Nanopore raw (R9.4, standard)",
    },
    "pacbio-hifi": {
        "flye_flag": "--pacbio-hifi",
        "minimap2_preset": "map-hifi",
        "medaka_model": None,
        "label": "PacBio HiFi (CCS)",
    },
}

# Default enhancement feature set
DEFAULT_ENHANCEMENTS = {
    "adaptive_params":    True,
    "iterative_polish":   True,
    "purge_dups":         False,   # only meaningful for diploid
    "confidence_scoring": True,
    "scaffolding":        False,   # opt-in (adds time)
    "telo_scaffolding":   False,   # opt-in: attach telomeric fragments to chromosomes
    "illumina_polish":    False,   # opt-in: Illumina short-read polishing after Medaka
}


# ------------------------------------------------
# helpers
# ------------------------------------------------

def run(cmd):
    print(f"\n[fungalflye] Running: {cmd}\n")
    subprocess.run(cmd, shell=True, check=True)


def parse_genome_size(g):
    g = str(g).strip().lower()
    if g.endswith("g"):
        return int(float(g[:-1]) * 1_000_000_000)
    if g.endswith("m"):
        return int(float(g[:-1]) * 1_000_000)
    if g.endswith("k"):
        return int(float(g[:-1]) * 1_000)
    return int(g)


def check_dependencies(read_type="nano-hq", enhancements=None):
    if enhancements is None:
        enhancements = DEFAULT_ENHANCEMENTS

    tools = ["flye", "minimap2", "seqkit", "filtlong"]

    cfg = READ_TYPE_CONFIGS[read_type]
    if cfg["medaka_model"]:
        tools.append("medaka_consensus")
    else:
        tools.append("racon")

    missing = [t for t in tools if sh.which(t) is None]

    if missing:
        print("\n❌ Missing required tools:")
        for m in missing:
            print(f"  - {m}")
        print("\nInstall with:")
        print("  conda install -c bioconda flye minimap2 seqkit filtlong")
        if "medaka_consensus" in missing:
            print("  pip install medaka")
        print()
        raise SystemExit(1)


# ------------------------------------------------
# pruning helpers
# ------------------------------------------------

def write_prune_settings(path, d):
    with open(path, "w") as f:
        json.dump(d, f, indent=2)


def load_prune_settings(path):
    if not Path(path).exists():
        return None
    with open(path) as f:
        return json.load(f)


def prune_contained_contigs(fasta, out_fasta, threads=8,
                             min_identity=0.95, min_coverage=0.95):
    print(f"\n[fungalflye] Pruning contained contigs "
          f"(identity>={min_identity}, coverage>={min_coverage})\n")
    fasta, out_fasta = Path(fasta), Path(out_fasta)
    paf = out_fasta.with_suffix(".self.paf")
    run(f"minimap2 -x asm5 -t {threads} {fasta} {fasta} > {paf}")
    remove = set()
    with open(paf) as f:
        for line in f:
            cols = line.strip().split()
            q, t = cols[0], cols[5]
            if q == t:
                continue
            qlen, tlen = int(cols[1]), int(cols[6])
            matches, aln_len = int(cols[9]), int(cols[10])
            identity = matches / aln_len if aln_len > 0 else 0
            coverage = min(aln_len / qlen, 1.0) if qlen > 0 else 0
            if identity >= min_identity and coverage >= min_coverage and qlen < tlen:
                remove.add(q)
    kept, removed_count = [], 0
    for r in SeqIO.parse(str(fasta), "fasta"):
        if r.id in remove:
            removed_count += 1
        else:
            kept.append(r)
    SeqIO.write(kept, str(out_fasta), "fasta")
    paf.unlink(missing_ok=True)
    print(f"[fungalflye] Removed {removed_count} contained contigs")
    return removed_count


def prune_small_contigs(fasta, out_fasta, min_size=5000):
    print(f"\n[fungalflye] Removing contigs < {min_size} bp\n")
    fasta, out_fasta = Path(fasta), Path(out_fasta)
    kept, removed_count = [], 0
    for r in SeqIO.parse(str(fasta), "fasta"):
        if len(r.seq) < min_size:
            removed_count += 1
        else:
            kept.append(r)
    SeqIO.write(kept, str(out_fasta), "fasta")
    print(f"[fungalflye] Removed {removed_count} small contigs")
    return removed_count


# ------------------------------------------------
# AT-rich warning
# ------------------------------------------------

def _warn_if_at_rich(fasta, sample=3):
    try:
        total, gc = 0, 0
        for i, r in enumerate(SeqIO.parse(str(fasta), "fasta")):
            if i >= sample:
                break
            s = str(r.seq).upper()
            total += len(s)
            gc += s.count("G") + s.count("C")
        if total > 0:
            pct = 100 * gc / total
            if pct < 35:
                print(f"\n⚠️  Low GC content ({pct:.1f}%) — "
                      "AT-rich genome detected. Consider --asm-coverage 80.\n")
    except Exception:
        pass


# ------------------------------------------------
# Mito separation
# ------------------------------------------------

def _separate_mito(polished_fasta, assembly_info, outdir):
    mito_ids = set()
    try:
        with open(assembly_info) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 4:
                    continue
                if parts[3].strip() == "Y" and 20_000 <= int(parts[1]) <= 100_000:
                    mito_ids.add(parts[0])
    except Exception:
        return
    if not mito_ids:
        return
    nuclear, mito = [], []
    for r in SeqIO.parse(str(polished_fasta), "fasta"):
        (mito if r.id in mito_ids else nuclear).append(r)
    if mito:
        SeqIO.write(nuclear, str(outdir / "nuclear.fasta"), "fasta")
        SeqIO.write(mito,    str(outdir / "mitochondrial.fasta"), "fasta")
        print(f"\n[fungalflye] Separated {len(mito)} mitochondrial contig(s)")


# ------------------------------------------------
# main assembly pipeline
# ------------------------------------------------

def run_assembly(
    reads,
    genome_size,
    outdir,
    threads=8,
    min_read_len=0,
    downsample_cov=0,
    min_contig_size=5000,
    prune_identity=0.95,
    prune_coverage=0.95,
    read_type="nano-hq",
    ploidy="haploid",
    asm_coverage=60,
    enhancements=None,
    illumina_r1=None,
    illumina_r2=None,
    illumina_polisher="polypolish",
):
    if enhancements is None:
        enhancements = dict(DEFAULT_ENHANCEMENTS)

    # Purge dups only makes sense for diploid
    if ploidy == "haploid":
        enhancements["purge_dups"] = False

    if read_type not in READ_TYPE_CONFIGS:
        raise ValueError(f"Unknown read_type '{read_type}'")

    cfg = READ_TYPE_CONFIGS[read_type]
    check_dependencies(read_type, enhancements)

    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)
    reads = Path(reads)

    flye_dir        = outdir / "flye"
    filtered_reads  = outdir / "reads.filtered.fastq"
    downsampled_reads = outdir / "reads.downsampled.fastq"
    contained_fasta = outdir / "contained_pruned.fasta"
    pruned_fasta    = outdir / "pruned.fasta"
    final_fasta     = outdir / "final.fasta"
    prune_settings_path = outdir / "prune_settings.json"

    reads_used = reads

    print("\n" + "=" * 60)
    print("🧬 FungalFlye Assembly Pipeline")
    print("=" * 60)
    print(f"Read type : {cfg['label']}")
    print(f"Ploidy    : {ploidy}")
    print(f"Polisher  : {'Medaka (iterative)' if cfg['medaka_model'] and enhancements.get('iterative_polish') else 'Medaka (1 round)' if cfg['medaka_model'] else 'Racon'}")
    print(f"Enhancements active: {[k for k, v in enhancements.items() if v]}")
    print("=" * 60)

    # FILTERING
    if min_read_len > 0:
        if filtered_reads.exists():
            print("[fungalflye] Found filtered reads — skipping")
            reads_used = filtered_reads
        else:
            print("\n[fungalflye] Filtering reads")
            run(f"seqkit seq -m {min_read_len} {reads} > {filtered_reads}")
            reads_used = filtered_reads

    # DOWNSAMPLING
    if downsample_cov > 0:
        if downsampled_reads.exists():
            print("[fungalflye] Found downsampled reads — skipping")
            reads_used = downsampled_reads
        else:
            genome_bp = parse_genome_size(genome_size)
            target_bases = genome_bp * downsample_cov
            run(f"filtlong --min_length 1000 --target_bases {target_bases} "
                f"{reads_used} > {downsampled_reads}")
            reads_used = downsampled_reads

    if not reads_used.exists():
        raise RuntimeError("Reads missing after preprocessing")

    # MODULE 1 — Adaptive parameters
    flye_params = {}
    if enhancements.get("adaptive_params"):
        genome_bp = parse_genome_size(genome_size)
        flye_params = suggest_flye_params(reads_used, genome_bp, threads)
        if "asm_coverage" in flye_params:
            asm_coverage = flye_params["asm_coverage"]

    # FLYE
    assembly = flye_dir / "assembly.fasta"
    if assembly.exists():
        print("[fungalflye] Existing Flye assembly detected — skipping")
    else:
        print(f"\n[fungalflye] Running Flye ({ploidy} mode)")
        min_overlap = flye_params.get("min_overlap", "")
        iterations  = flye_params.get("iterations", 3)
        overlap_flag = f"--min-overlap {min_overlap}" if min_overlap else ""
        hap_flag = "--keep-haplotypes" if ploidy == "diploid" else ""
        run(
            f"flye {cfg['flye_flag']} {reads_used} "
            f"--genome-size {genome_size} "
            f"--threads {threads} "
            f"--iterations {iterations} "
            f"--asm-coverage {asm_coverage} "
            f"{overlap_flag} {hap_flag} "
            f"-o {flye_dir}"
        )

    if not assembly.exists():
        raise RuntimeError("Flye assembly failed")

    _warn_if_at_rich(assembly)

    # POLISHING
    if cfg["medaka_model"]:
        if enhancements.get("iterative_polish"):
            polished_fasta = run_medaka_iterative(
                assembly=assembly,
                reads=reads_used,
                outdir=outdir,
                threads=threads,
                model=cfg["medaka_model"],
                max_rounds=3,
                convergence_threshold=50,
            )
        else:
            # Single-round Medaka fallback
            from .enhance import run_medaka_iterative as _med
            polished_fasta = _med(
                assembly=assembly, reads=reads_used, outdir=outdir,
                threads=threads, model=cfg["medaka_model"],
                max_rounds=1, convergence_threshold=0,
            )
    else:
        # PacBio — use Racon
        paf = outdir / "reads.paf"
        racon_fasta = outdir / "racon.fasta"
        if not racon_fasta.exists():
            run(f"minimap2 -x {cfg['minimap2_preset']} -t {threads} "
                f"{assembly} {reads_used} > {paf}")
            run(f"racon -t {threads} {reads_used} {paf} {assembly} > {racon_fasta}")
            paf.unlink(missing_ok=True)
        polished_fasta = racon_fasta

    # MODULE 6 — Illumina polishing (optional)
    if illumina_r1 and illumina_r2:
        polished_fasta = run_illumina_polishing(
            assembly=polished_fasta,
            illumina_r1=illumina_r1,
            illumina_r2=illumina_r2,
            outdir=outdir,
            threads=threads,
            polisher=illumina_polisher,
        )

    # MODULE 3 — Purge Duplicates
    if enhancements.get("purge_dups") and ploidy == "diploid":
        polished_fasta, _ = run_purge_dups(
            assembly=polished_fasta,
            reads=reads_used,
            outdir=outdir,
            threads=threads,
            minimap2_preset=cfg["minimap2_preset"],
        )

    # MODULE 5 — Scaffolding
    if enhancements.get("scaffolding"):
        polished_fasta = run_scaffold(
            assembly=polished_fasta,
            reads=reads_used,
            outdir=outdir,
            threads=threads,
            minimap2_preset=cfg["minimap2_preset"],
            min_support=3,
            end_window=2000,
        )

    # MODULE 6 — Telomere-guided scaffolding (runs after general scaffolding, before pruning)
    if enhancements.get("telo_scaffolding"):
        telo_motif = enhancements.get("telo_motif", "TTAGGG")
        telo_result = run_telomere_scaffolding(
            assembly=polished_fasta,
            reads=reads_used,
            outdir=outdir,
            threads=threads,
            minimap2_preset=cfg["minimap2_preset"],
            telomere_motif=telo_motif,
            min_support=2,
        )
        polished_fasta = telo_result

    # Mito separation
    assembly_info = flye_dir / "assembly_info.txt"
    if assembly_info.exists():
        _separate_mito(polished_fasta, assembly_info, outdir)

    # PRUNING
    current_settings = {
        "min_contig_size": int(min_contig_size),
        "prune_identity":  float(prune_identity),
        "prune_coverage":  float(prune_coverage),
    }
    old_settings = load_prune_settings(prune_settings_path)
    rerun = (
        not contained_fasta.exists()
        or not pruned_fasta.exists()
        or old_settings != current_settings
    )
    if rerun:
        prune_contained_contigs(
            polished_fasta, contained_fasta, threads, prune_identity, prune_coverage
        )
        prune_small_contigs(contained_fasta, pruned_fasta, min_contig_size)
        write_prune_settings(prune_settings_path, current_settings)
    else:
        print("[fungalflye] Existing pruned assembly with same settings — skipping")

    # FINAL
    if final_fasta.exists() and old_settings != current_settings:
        shutil.copy(pruned_fasta, final_fasta)
    elif not final_fasta.exists():
        shutil.copy(pruned_fasta, final_fasta)
    else:
        print("[fungalflye] Final assembly already exists")

    # MODULE 4 — Confidence scoring
    if enhancements.get("confidence_scoring"):
        score_contig_confidence(
            assembly=final_fasta,
            reads=reads_used,
            outdir=outdir,
            threads=threads,
            minimap2_preset=cfg["minimap2_preset"],
        )

    n_contigs = sum(1 for l in open(final_fasta) if l.startswith(">"))

    print("\n" + "=" * 60)
    print("✅ Assembly Complete")
    print("=" * 60)
    print(f"Final assembly : {final_fasta}")
    print(f"Contigs        : {n_contigs}")
    print(f"Contig cutoff  : {min_contig_size} bp")
    print(f"Ploidy         : {ploidy}")
    print("=" * 60 + "\n")

    return str(final_fasta)
