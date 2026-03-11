import subprocess
from pathlib import Path
import shutil
import shutil as sh
import json
from Bio import SeqIO


# ------------------------------------------------
# helper runner
# ------------------------------------------------

def run(cmd):
    print(f"\n[fungalflye] Running: {cmd}\n")
    subprocess.run(cmd, shell=True, check=True)


# ------------------------------------------------
# dependency check
# ------------------------------------------------

def check_dependencies():

    tools = ["flye", "minimap2", "racon", "seqkit", "filtlong"]

    missing = []

    for t in tools:
        if sh.which(t) is None:
            missing.append(t)

    if missing:
        print("\n❌ Missing required tools:")
        for m in missing:
            print(f"  - {m}")

        print("\nInstall with:")
        print("conda install -c bioconda flye minimap2 racon seqkit filtlong\n")

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

    run(f"minimap2 -x asm10 -t {threads} {fasta} {fasta} > {paf}")

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
            coverage = aln_len / qlen if qlen > 0 else 0

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


def prune_small_contigs(fasta, out_fasta, min_size=20000):

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
# main assembly pipeline (RESUMABLE)
# ------------------------------------------------

def run_assembly(
    reads,
    genome_size,
    outdir,
    threads=8,
    min_read_len=0,
    downsample_cov=0,
    min_contig_size=20000,
    prune_identity=0.95,
    prune_coverage=0.95,
):

    check_dependencies()

    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)

    reads = Path(reads)

    flye_dir = outdir / "flye"

    filtered_reads = outdir / "reads.filtered.fastq"
    downsampled_reads = outdir / "reads.downsampled.fastq"

    paf = outdir / "reads.paf"

    racon_fasta = outdir / "racon.fasta"
    contained_fasta = outdir / "contained_pruned.fasta"
    pruned_fasta = outdir / "pruned.fasta"
    final_fasta = outdir / "final.fasta"

    prune_settings_path = outdir / "prune_settings.json"

    reads_used = reads

    print("\n" + "=" * 60)
    print("🧬 FungalFlye Assembly Pipeline")
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

            g = str(genome_size).lower()

            if "m" in g:
                genome_bp = int(g.replace("m", "")) * 1_000_000
            else:
                genome_bp = int(g)

            target_bases = genome_bp * downsample_cov

            run(
                f"filtlong --target_bases {target_bases} "
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
        print("\n[fungalflye] Running Flye assembly")
        run(
            f"flye --nano-hq {reads_used} "
            f"--genome-size {genome_size} "
            f"--threads {threads} "
            f"--iterations 3 "
            f"--asm-coverage 60 "
            f"--keep-haplotypes "
            f"-o {flye_dir}"
        )

    if not assembly.exists():
        raise RuntimeError("Flye assembly failed")

    # ------------------------------------------------
    # MAPPING
    # ------------------------------------------------

    if paf.exists():
        print("[fungalflye] Existing read mapping detected — skipping")
    else:
        print("\n[fungalflye] Mapping reads")
        run(
            f"minimap2 -x map-ont -t {threads} "
            f"{assembly} {reads_used} > {paf}"
        )

    # ------------------------------------------------
    # RACON
    # ------------------------------------------------

    if racon_fasta.exists():
        print("[fungalflye] Existing racon polish detected — skipping")
    else:
        print("\n[fungalflye] Running Racon polishing")
        run(f"racon {reads_used} {paf} {assembly} > {racon_fasta}")

    if not racon_fasta.exists():
        raise RuntimeError("Racon polishing failed")

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
            racon_fasta,
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
        print("[fungalflye] Existing pruned assembly detected with same settings — skipping")

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
    print(f"Final assembly: {final_fasta}")
    print(f"Contigs: {n_contigs}")
    print(f"Small-contig cutoff used: {min_contig_size} bp")
    print("=" * 60 + "\n")

    return str(final_fasta)