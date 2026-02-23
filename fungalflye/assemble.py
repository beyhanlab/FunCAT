# assemble.py

import subprocess
from pathlib import Path
import shutil
import sys

# -------------------------

# helper runner

# -------------------------

def run(cmd):
print(f"\n[fungalflye] Running: {cmd}\n")
subprocess.run(cmd, shell=True, check=True)

# -------------------------

# dependency check

# -------------------------

def check_dependencies():

```
tools = ["flye", "minimap2", "racon", "seqkit"]

missing = []

for t in tools:
    if shutil.which(t) is None:
        missing.append(t)

if missing:
    print("\n❌ Missing required tools:")
    for m in missing:
        print(f"  - {m}")

    print("\nInstall with:")
    print("  conda env create -f environment.yml")
    print("  conda activate fungalflye\n")

    sys.exit(1)
```

# -------------------------

# containment pruning

# -------------------------

def prune_contained_contigs(fasta, out_fasta, threads=8):

```
print("\n[fungalflye] Pruning redundant contigs (>95% contained)\n")

fasta = Path(fasta)
out_fasta = Path(out_fasta)

paf = out_fasta.with_suffix(".self.paf")

run(
    f"minimap2 -x asm10 -t {threads} {fasta} {fasta} > {paf}"
)

remove = set()
lengths = {}

# collect lengths
with open(fasta) as f:
    name = None
    seq = []
    for line in f:
        if line.startswith(">"):
            if name:
                lengths[name] = len("".join(seq))
            name = line[1:].split()[0]
            seq = []
        else:
            seq.append(line.strip())
    if name:
        lengths[name] = len("".join(seq))

# parse alignments
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

# write filtered fasta
with open(out_fasta, "w") as out:

    keep = True

    with open(fasta) as f:
        for line in f:

            if line.startswith(">"):
                name = line[1:].split()[0]
                keep = name not in remove

            if keep:
                out.write(line)

print(f"[fungalflye] Removed {len(remove)} redundant contigs")

paf.unlink(missing_ok=True)
```

# -------------------------

# main assembly pipeline

# -------------------------

def run_assembly(
reads,
genome_size,
outdir,
threads=8,
min_read_len=3000,
downsample_cov=0,
):

```
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

# -------------------------
# resume detection
# -------------------------

if final_fasta.exists():
    print("\n[fungalflye] Existing final assembly detected.")
    print(f"Path: {final_fasta}")
    print("Skipping assembly (resume mode)\n")
    return str(final_fasta)

# -------------------------
# read filtering
# -------------------------

print("\n[fungalflye] Preparing reads\n")

reads_used = reads

if min_read_len > 0:

    if str(reads) != str(filtered_reads):

        run(
            f"seqkit seq -m {min_read_len} {reads} > {filtered_reads}"
        )

        reads_used = filtered_reads

# -------------------------
# optional downsampling
# -------------------------

if downsample_cov > 0:

    print("\n[fungalflye] Downsampling reads\n")

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

if not reads_used.exists() or reads_used.stat().st_size == 0:
    raise RuntimeError(
        "Reads file is empty after filtering/downsampling."
    )

# -------------------------
# Flye assembly
# -------------------------

print("\n[fungalflye] Running Flye assembly\n")

run(
    f"flye --nano-hq {reads_used} "
    f"--genome-size {genome_size} "
    f"--threads {threads} "
    f"--iterations 3 "
    f"--asm-coverage 60 "
    f"--keep-haplotypes "
    f"-o {flye_dir}"
)

assembly = flye_dir / "assembly.fasta"

if not assembly.exists():
    raise RuntimeError("Flye assembly not produced.")

# -------------------------
# mapping reads
# -------------------------

print("\n[fungalflye] Mapping reads\n")

run(
    f"minimap2 -x map-ont -t {threads} "
    f"{assembly} {reads_used} > {paf}"
)

# -------------------------
# racon polishing
# -------------------------

print("\n[fungalflye] Running Racon polishing\n")

run(
    f"racon {reads_used} {paf} {assembly} > {racon_fasta}"
)

# -------------------------
# prune redundant contigs
# -------------------------

prune_contained_contigs(
    racon_fasta,
    pruned_fasta,
    threads=threads
)

shutil.copy(pruned_fasta, final_fasta)

n_contigs = sum(
    1 for line in open(final_fasta) if line.startswith(">")
)

print("\n" + "=" * 50)
print("🧬 Assembly Complete")
print("=" * 50)
print(f"Final assembly: {final_fasta}")
print(f"Contigs: {n_contigs}")
print("=" * 50 + "\n")

return str(final_fasta)
```
