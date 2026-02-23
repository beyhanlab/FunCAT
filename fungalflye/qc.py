import subprocess
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from Bio import SeqIO
from collections import Counter
import shutil

# -------------------------

# helper: run shell command

# -------------------------

def run(cmd):
print(f"\n[fungalflye] Running: {cmd}\n")
subprocess.run(cmd, shell=True, check=True)

# -------------------------

# dependency check

# -------------------------

def check_dependencies():

```
tools = ["seqkit"]

missing = []

for t in tools:
    if shutil.which(t) is None:
        missing.append(t)

if missing:
    print("\n❌ Missing required tools:")
    for m in missing:
        print(f"  - {m}")

    print("\nInstall with:")
    print("  conda install -c bioconda seqkit\n")

    raise SystemExit(1)
```

# -------------------------

# telomere utilities

# -------------------------

def revcomp(s):
table = str.maketrans("ACGT", "TGCA")
return s.translate(table)[::-1]

def hamming(a, b):
return sum(x != y for x, y in zip(a, b))

def scan_window(seq, motif, max_mismatch):

```
m = len(motif)
hits = 0
best = m

for i in range(len(seq) - m + 1):

    w = seq[i:i+m]
    d = hamming(w, motif)

    if d <= max_mismatch:
        hits += 1
        best = min(best, d)

return hits, best if hits else None
```

def analyze_end(seq, motif, window, max_mismatch):

```
rc = revcomp(motif)

hits1, best1 = scan_window(seq, motif, max_mismatch)
hits2, best2 = scan_window(seq, rc, max_mismatch)

hits = hits1 + hits2

best = min(
    x for x in [best1, best2] if x is not None
) if hits else None

density = round(hits / (window / 1000), 2)

return hits, best, density
```

# -------------------------

# tandem repeat detection

# -------------------------

def max_tandem_run(seq, motif, max_mismatch=0):

```
seq = seq.upper()
motif = motif.upper()
m = len(motif)

best = 0

for i in range(0, len(seq) - m + 1):

    run = 0
    j = i

    while j + m <= len(seq):

        unit = seq[j:j+m]

        if hamming(unit, motif) <= max_mismatch:
            run += 1
            j += m
        else:
            break

    if run > best:
        best = run

return best
```

def tandem_metrics(seq, motif, max_mismatch=0):

```
rc = revcomp(motif)

r1 = max_tandem_run(seq, motif, max_mismatch)
r2 = max_tandem_run(seq, rc, max_mismatch)

if r1 >= r2:
    return r1, "motif", r1 * len(motif)
else:
    return r2, "revcomp", r2 * len(motif)
```

# -------------------------

# telomere discovery

# -------------------------

def discover_telomere_motif(fasta, k=6, window=3000):

```
print("\n[fungalflye] Discovering telomere motif...")

kmer_counts = Counter()

for record in SeqIO.parse(fasta, "fasta"):

    seq = str(record.seq).upper()
    ends = seq[:window] + seq[-window:]

    for i in range(len(ends) - k + 1):

        kmer = ends[i:i+k]

        if "N" not in kmer:
            kmer_counts[kmer] += 1

most_common = kmer_counts.most_common(10)

print("\nTop candidate telomere kmers:")
for kmer, count in most_common:
    print(f"  {kmer}  ({count})")

best = most_common[0][0]

print(f"\n[fungalflye] Selected motif: {best}")

return best
```

# -------------------------

# telomere scan

# -------------------------

def scan_telomeres(
fasta,
motif,
window=5000,
max_mismatch=2,
tandem_mismatch=2,
min_tandem_repeats=3
):

```
rows = []

for record in SeqIO.parse(fasta, "fasta"):

    seq = str(record.seq).upper()

    for side, piece in [
        ("start", seq[:window]),
        ("end", seq[-window:])
    ]:

        hits, best, density = analyze_end(
            piece, motif, window, max_mismatch
        )

        max_rep, orient, run_bp = tandem_metrics(
            piece, motif, tandem_mismatch
        )

        telomeric = "YES" if max_rep >= min_tandem_repeats else "NO"

        rows.append([
            record.id,
            side,
            hits,
            best,
            density,
            max_rep,
            run_bp,
            orient,
            telomeric
        ])

df = pd.DataFrame(
    rows,
    columns=[
        "contig",
        "side",
        "hits",
        "best_mismatch",
        "density_per_kb",
        "max_consecutive_repeats",
        "max_run_bp",
        "orientation",
        "telomeric"
    ]
)

return df
```

# -------------------------

# pretty report

# -------------------------

def print_assembly_report(fasta, lengths, telo_df=None):

```
total_size = sum(lengths)
contigs = len(lengths)
largest = max(lengths)

sorted_lengths = sorted(lengths, reverse=True)

cumsum = 0
n50 = 0
l50 = 0

for i, L in enumerate(sorted_lengths, 1):
    cumsum += L
    if cumsum >= total_size / 2:
        n50 = L
        l50 = i
        break

print("\n" + "=" * 60)
print("🧬 FungalFlye Assembly Report")
print("=" * 60)

print(f"\nAssembly file: {fasta}\n")

print(f"Contigs:                {contigs}")
print(f"Total size:             {total_size:,} bp")
print(f"Largest contig:         {largest:,} bp")
print(f"N50:                    {n50:,} bp")
print(f"L50:                    {l50}")

if telo_df is not None:

    total_ends = len(telo_df)
    telomeric = (telo_df["telomeric"] == "YES").sum()

    both = 0
    for c in telo_df["contig"].unique():
        sub = telo_df[telo_df["contig"] == c]
        if all(sub["telomeric"] == "YES"):
            both += 1

    density = telo_df["density_per_kb"].mean()

    print(f"\nTelomeric ends:         {telomeric} / {total_ends}")
    print(f"Chromosomes complete:   {both}")
    print(f"Mean telomere density:  {density:.2f} hits/kb")

print("\n" + "=" * 60)

if contigs <= 12:
    print("✅ Assembly appears chromosome-level complete")
else:
    print("⚠️ Assembly fragmented — consider tuning parameters")

print("=" * 60 + "\n")
```

# -------------------------

# main QC

# -------------------------

def run_qc(fasta, telomere=None, run_telomeres=False):

```
check_dependencies()

fasta = Path(fasta)

print("\n[fungalflye] Starting QC...")

# output directory next to fasta
outdir = fasta.parent / "fungalflye_qc"
outdir.mkdir(exist_ok=True)

stats_file = outdir / "stats.txt"
lengths_file = outdir / "lengths.tsv"

# stats
run(f"seqkit stats {fasta} > {stats_file}")
run(f"seqkit fx2tab -n -l {fasta} > {lengths_file}")

df = pd.read_csv(lengths_file, sep="\t", header=None)
lengths = df[1].tolist()

# histogram
plt.figure(figsize=(6, 4))
plt.hist(lengths, bins=30)
plt.xlabel("Contig length (bp)")
plt.ylabel("Count")
plt.title("Contig length distribution")
plt.tight_layout()
plt.savefig(outdir / "length_histogram.png")
plt.close()

telo_df = None

# telomeres
if run_telomeres and telomere:

    print(f"\n[fungalflye] Scanning telomeres using motif: {telomere}")

    telo_df = scan_telomeres(fasta, telomere)

    telo_df.to_csv(outdir / "telomeres.tsv", sep="\t", index=False)

elif run_telomeres:
    print("\n[fungalflye] WARNING: Telomere scan requested but motif missing")

# final report
print_assembly_report(fasta, lengths, telo_df)

print(f"\n✅ QC complete. Results in: {outdir}\n")
```
