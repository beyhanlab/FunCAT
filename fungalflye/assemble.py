import subprocess
from pathlib import Path
import shutil
import shutil as sh
import json
from Bio import SeqIO


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


# ------------------------------------------------
# helper runner
# ------------------------------------------------

def run(cmd):
    print(f"\n[fungalflye] Running: {cmd}\n")
    subprocess.run(cmd, shell=True, check=True)


# ------------------------------------------------
# genome size parser (handles 40m, 40M, 1g, 1G, raw bp)
# ------------------------------------------------

def parse_genome_size(g):
    """Return genome size in base pairs as an integer."""
    g = str(g).strip().lower()
    if g.endswith("g"):
        return int(float(g[:-1]) * 1_000_000_000)
    if g.endswith("m"):
        return int(float(g[:-1]) * 1_000_000)
    if g.endswith("k"):
        return int(float(g[:-1]) * 1_000)
    return int(g)


# ------------------------------------------------
# dependency check
# ------------------------------------------------

def check_dependencies(read_type="nano-hq", use_medaka=True):

    tools = ["flye", "minimap2", "seqkit", "filtlong"]

    if use_medaka and READ_TYPE_CONFIGS[read_type]["medaka_model"] is not None:
        tools.append("medaka_consensus")
    else:
        tools.append("racon")

    missing = []

    for t in tools:
        if sh.which(t) is None:
            missing.append(t)

    if missing:
        print("\n❌ Missing required tools:")
        for m in missing:
            print(f"  - {m}")

        print("\nInstall with:")
        print("  conda install -c bioconda flye minimap2 seqkit filtlong")
        if "medaka_consensus" in missing:
            print("  pip install medaka")
        if "racon" in missing:
            print("  conda install -c bioconda racon")
        print()

        raise SystemExit(1)


# ------------------------------------------------
# pruning helpers
# ------------------------------------------------

def write_prune_settings(settings_path, settings_dict):
    with open(settings_path, "w") as fh:
        json.dump(settings_dict, fh, indent=2)


def load_prune_settings(settings_path):
    if not Path(settings_path).exists():
        return None
    with open(settings_path) as fh:
        return json.load(fh)


def prune_contained_contigs(
    fasta,
    out_fasta,
    threads=8,
    min_identity=0.95,
    min_coverage=0.95
):

    print(
        f"\n[fungalflye] Pruning contained contigs "
        f"(identity >= {min_identity}, coverage >= {min_coverage})\n"
    )

    fasta = Path(fasta)
    out_fasta = Path(out_fasta)

    paf = out_fasta.with_suffix(".self.paf")

    # asm5 is correct for self-alignment contained-contig detection
    run(f"minimap2 -x asm5 -t {threads} {fasta} {fasta} > {paf}")

    remove = set()

    with open(paf) as f:
        for line in f:
            cols = line.strip().split()

            q = cols[0]
            t = cols[5]

            if q == t:
                continue

            qlen = int(cols[1])
            tlen = int(cols[6])

            matches = int(cols[9])
            aln_len = int(cols[10])

            identity = matches / aln_len if aln_len > 0 else 0
            # cap coverage at 1.0 to avoid false positives from indels
            coverage = min(aln_len / qlen, 1.0) if qlen > 0 else 0

            if identity >= min_identity and coverage >= min_coverage:
                if qlen < tlen:
                    remove.add(q)

    kept_records = []
    removed_count = 0

    for record in SeqIO.parse(str(fasta), "fasta"):
        if record.id in remove:
            removed_count += 1
        else:
            kept_records.append(record)

    SeqIO.write(kept_records, str(out_fasta), "fasta")

    paf.unlink(missing_ok=True)

    print(f"[fungalflye] Removed {removed_count} contained contigs")

    return removed_count


def prune_small_contigs(fasta, out_fasta, min_size=5000):

    print(f"\n[fungalflye] Removing contigs smaller than {min_size} bp\n")

    fasta = Path(fasta)
    out_fasta = Path(out_fasta)

    kept_records = []
    removed_count = 0

    for record in SeqIO.parse(str(fasta), "fasta"):
        if len(record.seq) < min_size:
            removed_count += 1
        else:
            kept_records.append(record)

    SeqIO.write(kept_records, str(out_fasta), "fasta")

    print(f"[fungalflye] Removed {removed_count} small contigs")

    return removed_count


# ------------------------------------------------
# Medaka polishing
# ------------------------------------------------

def run_medaka(assembly, reads, outdir, threads, model):
    """Run Medaka consensus polishing. Returns path to polished FASTA."""

    medaka_dir = outdir / "medaka"

    polished = medaka_dir / "consensus.fasta"

    if polished.exists():
        print("[fungalflye] Existing Medaka polish detected — skipping")
        return polished

    print(f"\n[fungalflye] Running Medaka polishing (model: {model})\n")

    medaka_dir.mkdir(exist_ok=True)

    run(
        f"medaka_consensus "
        f"-i {reads} "
        f"-d {assembly} "
        f"-o {medaka_dir} "
        f"-t {threads} "
        f"-m {model}"
    )

    if not polished.exists():
        raise RuntimeError("Medaka polishing failed — consensus.fasta not found")

    return polished


# ------------------------------------------------
# Racon polishing (fallback for PacBio HiFi or if medaka unavailable)
# ------------------------------------------------

def run_racon(assembly, reads, outdir, threads, minimap2_preset):
    """Run one round of Racon polishing. Returns path to polished FASTA."""

    paf = outdir / "reads.paf"
    racon_fasta = outdir / "racon.fasta"

    if racon_fasta.exists():
        print("[fungalflye] Existing Racon polish detected — skipping")
        return racon_fasta

    print("\n[fungalflye] Mapping reads for Racon")

    run(
        f"minimap2 -x {minimap2_preset} -t {threads} "
        f"{assembly} {reads} > {paf}"
    )

    print("\n[fungalflye] Running Racon polishing")

    # Note: -t threads flag added
    run(f"racon -t {threads} {reads} {paf} {assembly} > {racon_fasta}")

    paf.unlink(missing_ok=True)

    if not racon_fasta.exists():
        raise RuntimeError("Racon polishing failed")

    return racon_fasta


# ------------------------------------------------
# main assembly pipeline (RESUMABLE)
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
):

    if read_type not in READ_TYPE_CONFIGS:
        raise ValueError(
            f"Unknown read_type '{read_type}'. "
            f"Choose from: {list(READ_TYPE_CONFIGS)}"
        )

    cfg = READ_TYPE_CONFIGS[read_type]

    # PacBio HiFi doesn't use Medaka — use Racon instead
    use_medaka = (cfg["medaka_model"] is not None)

    check_dependencies(read_type=read_type, use_medaka=use_medaka)

    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)

    reads = Path(reads)

    flye_dir = outdir / "flye"

    filtered_reads = outdir / "reads.filtered.fastq"
    downsampled_reads = outdir / "reads.downsampled.fastq"

    contained_fasta = outdir / "contained_pruned.fasta"
    pruned_fasta = outdir / "pruned.fasta"
    final_fasta = outdir / "final.fasta"

    prune_settings_path = outdir / "prune_settings.json"

    reads_used = reads

    print("\n" + "=" * 60)
    print("🧬 FungalFlye Assembly Pipeline")
    print("=" * 60)
    print(f"Read type : {cfg['label']}")
    print(f"Ploidy    : {ploidy}")
    print(f"Polisher  : {'Medaka' if use_medaka else 'Racon'}")
    print("=" * 60)

    # ------------------------------------------------
    # FILTERING
    # ------------------------------------------------

    if min_read_len > 0:

        if filtered_reads.exists():
            print("[fungalflye] Found filtered reads — skipping")
            reads_used = filtered_reads
        else:
            print("\n[fungalflye] Filtering reads")
            run(f"seqkit seq -m {min_read_len} {reads} > {filtered_reads}")
            reads_used = filtered_reads

    # ------------------------------------------------
    # DOWNSAMPLING
    # ------------------------------------------------

    if downsample_cov > 0:

        if downsampled_reads.exists():
            print("[fungalflye] Found downsampled reads — skipping")
            reads_used = downsampled_reads
        else:
            print("\n[fungalflye] Downsampling reads")

            genome_bp = parse_genome_size(genome_size)
            target_bases = genome_bp * downsample_cov

            run(
                f"filtlong --min_length 1000 "
                f"--target_bases {target_bases} "
                f"{reads_used} > {downsampled_reads}"
            )

            reads_used = downsampled_reads

    if not reads_used.exists():
        raise RuntimeError("Reads missing after preprocessing")

    # ------------------------------------------------
    # FLYE
    # ------------------------------------------------

    assembly = flye_dir / "assembly.fasta"

    if assembly.exists():
        print("[fungalflye] Existing Flye assembly detected — skipping")
    else:
        print(f"\n[fungalflye] Running Flye assembly ({ploidy} mode)")

        flye_cmd = (
            f"flye {cfg['flye_flag']} {reads_used} "
            f"--genome-size {genome_size} "
            f"--threads {threads} "
            f"--iterations 3 "
            f"--asm-coverage {asm_coverage} "
            f"-o {flye_dir}"
        )

        # Only keep haplotypes for diploid assemblies
        if ploidy == "diploid":
            flye_cmd += " --keep-haplotypes"

        run(flye_cmd)

    if not assembly.exists():
        raise RuntimeError("Flye assembly failed")

    # ------------------------------------------------
    # GC WARNING for AT-rich genomes
    # ------------------------------------------------

    _warn_if_at_rich(assembly)

    # ------------------------------------------------
    # POLISHING
    # ------------------------------------------------

    if use_medaka:
        polished_fasta = run_medaka(
            assembly=assembly,
            reads=reads_used,
            outdir=outdir,
            threads=threads,
            model=cfg["medaka_model"],
        )
    else:
        polished_fasta = run_racon(
            assembly=assembly,
            reads=reads_used,
            outdir=outdir,
            threads=threads,
            minimap2_preset=cfg["minimap2_preset"],
        )

    # ------------------------------------------------
    # MITO SEPARATION
    # ------------------------------------------------

    assembly_info = flye_dir / "assembly_info.txt"
    if assembly_info.exists():
        _separate_mito(polished_fasta, assembly_info, outdir)

    # ------------------------------------------------
    # PRUNING
    # ------------------------------------------------

    current_prune_settings = {
        "min_contig_size": int(min_contig_size),
        "prune_identity": float(prune_identity),
        "prune_coverage": float(prune_coverage),
    }

    old_prune_settings = load_prune_settings(prune_settings_path)

    rerun_pruning = (
        (not contained_fasta.exists())
        or (not pruned_fasta.exists())
        or (old_prune_settings != current_prune_settings)
    )

    if rerun_pruning:
        print("\n[fungalflye] Running pruning workflow")

        prune_contained_contigs(
            polished_fasta,
            contained_fasta,
            threads=threads,
            min_identity=prune_identity,
            min_coverage=prune_coverage,
        )

        prune_small_contigs(
            contained_fasta,
            pruned_fasta,
            min_size=min_contig_size,
        )

        write_prune_settings(prune_settings_path, current_prune_settings)

    else:
        print("[fungalflye] Existing pruned assembly with same settings — skipping")

    # ------------------------------------------------
    # FINAL
    # ------------------------------------------------

    if final_fasta.exists():
        if old_prune_settings != current_prune_settings:
            shutil.copy(pruned_fasta, final_fasta)
        else:
            print("[fungalflye] Final assembly already exists")
    else:
        shutil.copy(pruned_fasta, final_fasta)

    n_contigs = sum(
        1 for line in open(final_fasta) if line.startswith(">")
    )

    print("\n" + "=" * 60)
    print("✅ Assembly Complete")
    print("=" * 60)
    print(f"Final assembly : {final_fasta}")
    print(f"Contigs        : {n_contigs}")
    print(f"Contig cutoff  : {min_contig_size} bp")
    print(f"Ploidy mode    : {ploidy}")
    print("=" * 60 + "\n")

    return str(final_fasta)


# ------------------------------------------------
# AT-rich genome warning
# ------------------------------------------------

def _warn_if_at_rich(fasta, sample_contigs=3):
    """Scan the first few contigs and warn if GC < 35%."""
    try:
        total, gc = 0, 0
        for i, rec in enumerate(SeqIO.parse(str(fasta), "fasta")):
            if i >= sample_contigs:
                break
            s = str(rec.seq).upper()
            total += len(s)
            gc += s.count("G") + s.count("C")
        if total > 0:
            gc_pct = 100 * gc / total
            if gc_pct < 35:
                print(
                    f"\n⚠️  Low GC content detected ({gc_pct:.1f}%). "
                    "AT-rich genomes may require adjusted Flye parameters. "
                    "Consider --asm-coverage 80 or higher.\n"
                )
    except Exception:
        pass


# ------------------------------------------------
# Mitochondrial contig separation
# ------------------------------------------------

def _separate_mito(polished_fasta, assembly_info, outdir):
    """
    Parse Flye assembly_info.txt and separate circular contigs
    in the mito size range (20–100 kb) into a separate FASTA.
    """
    mito_ids = set()

    try:
        with open(assembly_info) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 4:
                    continue
                name = parts[0]
                length = int(parts[1])
                circular = parts[3].strip() == "Y"
                if circular and 20_000 <= length <= 100_000:
                    mito_ids.add(name)
    except Exception:
        return

    if not mito_ids:
        return

    nuclear_fasta = outdir / "nuclear.fasta"
    mito_fasta = outdir / "mitochondrial.fasta"

    nuclear = []
    mito = []

    for rec in SeqIO.parse(str(polished_fasta), "fasta"):
        if rec.id in mito_ids:
            mito.append(rec)
        else:
            nuclear.append(rec)

    if mito:
        SeqIO.write(nuclear, str(nuclear_fasta), "fasta")
        SeqIO.write(mito, str(mito_fasta), "fasta")
        print(
            f"\n[fungalflye] Separated {len(mito)} putative mitochondrial "
            f"contig(s) → {mito_fasta}"
        )
        print(f"[fungalflye] Nuclear contigs → {nuclear_fasta}\n")
