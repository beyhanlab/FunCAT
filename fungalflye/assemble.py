# assemble.py

import subprocess
from pathlib import Path
import shutil
import shutil as sh


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

    tools = ["flye", "minimap2", "racon", "seqkit"]

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
# prune redundant contigs
# ------------------------------------------------

def prune_contained_contigs(fasta, out_fasta, threads=8):

    print("\n[fungalflye] Pruning redundant contigs (>95% contained)\n")

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

            if identity > 0.95 and coverage > 0.95:
                if qlen < tlen:
                    remove.add(q)

    with open(out_fasta, "w") as out:
        keep = True
        with open(fasta) as f:
            for line in f:
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    keep = name not in remove
                if keep:
                    out.write(line)

    paf.unlink(missing_ok=True)

    print(f"[fungalflye] Removed {len(remove)} redundant contigs")


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
    pruned_fasta = outdir / "pruned.fasta"
    final_fasta = outdir / "final.fasta"

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

    # safety
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
    # RACON POLISHING
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

    if pruned_fasta.exists():
        print("[fungalflye] Existing pruned assembly detected — skipping")
    else:
        prune_contained_contigs(racon_fasta, pruned_fasta, threads)

    # ------------------------------------------------
    # FINAL
    # ------------------------------------------------

    if final_fasta.exists():
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
    print("=" * 60 + "\n")

    return str(final_fasta)