import subprocess
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from Bio import SeqIO
from collections import Counter
import shutil


def run(cmd):
    print(f"\n[fungalflye] Running: {cmd}\n")
    subprocess.run(cmd, shell=True, check=True)


def check_dependencies():
    tools = ["seqkit"]
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        print("\n❌ Missing required tools:")
        for m in missing:
            print(f"  - {m}")
        print("\nInstall with: conda install -c bioconda seqkit\n")
        raise SystemExit(1)


def revcomp(s):
    table = str.maketrans("ACGT", "TGCA")
    return s.translate(table)[::-1]


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def scan_window(seq, motif, max_mismatch):
    m = len(motif)
    hits, best = 0, m
    for i in range(len(seq) - m + 1):
        d = hamming(seq[i:i+m], motif)
        if d <= max_mismatch:
            hits += 1
            best = min(best, d)
    return hits, best if hits else None


def analyze_end(seq, motif, window, max_mismatch):
    rc = revcomp(motif)
    hits1, best1 = scan_window(seq, motif, max_mismatch)
    hits2, best2 = scan_window(seq, rc,    max_mismatch)
    hits = hits1 + hits2
    best = min(x for x in [best1, best2] if x is not None) if hits else None
    density = round(hits / (window / 1000), 2)
    return hits, best, density


def max_tandem_run(seq, motif, max_mismatch=0):
    seq, motif = seq.upper(), motif.upper()
    m, best = len(motif), 0
    for i in range(0, len(seq) - m + 1):
        run, j = 0, i
        while j + m <= len(seq):
            if hamming(seq[j:j+m], motif) <= max_mismatch:
                run += 1; j += m
            else:
                break
        if run > best:
            best = run
    return best


def tandem_metrics(seq, motif, max_mismatch=0):
    rc = revcomp(motif)
    r1 = max_tandem_run(seq, motif, max_mismatch)
    r2 = max_tandem_run(seq, rc,    max_mismatch)
    if r1 >= r2:
        return r1, "motif",   r1 * len(motif)
    return r2, "revcomp", r2 * len(motif)


def discover_telomere_motif(fasta, k=6, window=3000):
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


def scan_telomeres(fasta, motif, window=5000, max_mismatch=2,
                   tandem_mismatch=2, min_tandem_repeats=3):
    rows = []
    for record in SeqIO.parse(fasta, "fasta"):
        seq = str(record.seq).upper()
        for side, piece in [("start", seq[:window]), ("end", seq[-window:])]:
            hits, best, density = analyze_end(piece, motif, window, max_mismatch)
            max_rep, orient, run_bp = tandem_metrics(piece, motif, tandem_mismatch)
            telomeric = "YES" if max_rep >= min_tandem_repeats else "NO"
            rows.append([record.id, side, hits, best, density,
                         max_rep, run_bp, orient, telomeric])
    return pd.DataFrame(rows, columns=[
        "contig", "side", "hits", "best_mismatch", "density_per_kb",
        "max_consecutive_repeats", "max_run_bp", "orientation", "telomeric"
    ])


def print_assembly_report(fasta, lengths, telo_df=None):
    total_size = sum(lengths)
    sorted_l   = sorted(lengths, reverse=True)
    cumsum, n50, l50 = 0, 0, 0
    for i, L in enumerate(sorted_l, 1):
        cumsum += L
        if cumsum >= total_size / 2:
            n50, l50 = L, i
            break
    print("\n" + "=" * 60)
    print("🧬 FungalFlye Assembly Report")
    print("=" * 60)
    print(f"\nAssembly file : {fasta}\n")
    print(f"Contigs       : {len(lengths)}")
    print(f"Total size    : {total_size:,} bp")
    print(f"Largest contig: {max(lengths):,} bp")
    print(f"N50           : {n50:,} bp")
    print(f"L50           : {l50}")
    if telo_df is not None:
        total_ends = len(telo_df)
        telomeric  = (telo_df["telomeric"] == "YES").sum()
        both = sum(
            1 for c in telo_df["contig"].unique()
            if all(telo_df[telo_df["contig"] == c]["telomeric"] == "YES")
        )
        print(f"\nTelomeric ends       : {telomeric} / {total_ends}")
        print(f"Chromosomes complete : {both}")
        print(f"Mean telomere density: {telo_df['density_per_kb'].mean():.2f} hits/kb")
    print("\n" + "=" * 60)
    verdict = "✅ Assembly appears chromosome-level complete" if len(lengths) <= 50 \
              else "⚠️  Assembly fragmented — consider tuning parameters or enabling scaffolding"
    print(verdict)
    print("=" * 60 + "\n")


def run_qc(fasta, telomere=None, run_telomeres=False,
           report=True, run_metadata=None):
    """
    Parameters
    ----------
    fasta         : path to assembly FASTA
    telomere      : motif string, or None (auto-discovered if run_telomeres=True)
    run_telomeres : run telomere scan
    report        : generate HTML report (default True)
    run_metadata  : dict passed to report generator {assembly_name, ploidy, ...}
    """
    check_dependencies()
    fasta = Path(fasta)
    print("\n[fungalflye] Starting QC...")

    outdir = fasta.parent / "fungalflye_qc"
    outdir.mkdir(exist_ok=True)

    run(f"seqkit stats {fasta} > {outdir / 'stats.txt'}")
    run(f"seqkit fx2tab -n -l {fasta} > {outdir / 'lengths.tsv'}")

    df = pd.read_csv(outdir / "lengths.tsv", sep="\t", header=None)
    lengths = df[1].tolist()

    plt.figure(figsize=(6, 4))
    plt.hist(lengths, bins=30)
    plt.xlabel("Contig length (bp)")
    plt.ylabel("Count")
    plt.title("Contig length distribution")
    plt.tight_layout()
    plt.savefig(outdir / "length_histogram.png")
    plt.close()

    telo_df = None
    if run_telomeres:
        if telomere is None:
            telomere = discover_telomere_motif(str(fasta))
        print(f"\n[fungalflye] Scanning telomeres using motif: {telomere}")
        telo_df = scan_telomeres(str(fasta), telomere)
        telo_df.to_csv(outdir / "telomeres.tsv", sep="\t", index=False)

    print_assembly_report(fasta, lengths, telo_df)
    print(f"\n✅ QC complete. Results in: {outdir}\n")

    # HTML report
    if report:
        confidence_tsv = fasta.parent / "confidence" / "contig_confidence.tsv"
        from .report import generate_report
        generate_report(
            fasta=fasta,
            outdir=fasta.parent,
            run_metadata=run_metadata or {},
            telo_df=telo_df,
            confidence_tsv=confidence_tsv if confidence_tsv.exists() else None,
        )
