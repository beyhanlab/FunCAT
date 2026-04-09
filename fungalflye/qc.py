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
    print(f"\n[funcat] Running: {cmd}\n")
    subprocess.run(cmd, shell=True, check=True)


def check_dependencies(assembly_mode=False):
    """
    Check all required and optional dependencies.
    assembly_mode=True checks assembly tools (flye, medaka, minimap2, etc.)
    Always checks core tools needed for QC.
    """

    # tool -> (install command, required)
    CORE = {
        "seqkit":   ("conda install -c bioconda seqkit",         True),
        "minimap2": ("conda install -c bioconda minimap2",        True),
        "samtools": ("conda install -c bioconda samtools",        True),
    }

    ASSEMBLY = {
        "flye":     ("conda install -c bioconda flye",            True),
        "filtlong": ("conda install -c bioconda filtlong",        True),
        "racon":    ("conda install -c bioconda racon",           False),  # only for PacBio
        "medaka_consensus": ("pip install medaka",                True),
    }

    OPTIONAL = {
        "purge_dups": ("conda install -c bioconda purge_dups",   False),
        "nucmer":     ("conda install -c bioconda mummer",        False),
    }

    to_check = dict(CORE)
    if assembly_mode:
        to_check.update(ASSEMBLY)

    missing_required = []
    missing_optional = []

    for tool, (install_cmd, required) in to_check.items():
        if shutil.which(tool) is None:
            if required:
                missing_required.append((tool, install_cmd))
            else:
                missing_optional.append((tool, install_cmd))

    # Also check optional tools and warn (don't fail)
    for tool, (install_cmd, _) in OPTIONAL.items():
        if shutil.which(tool) is None:
            missing_optional.append((tool, install_cmd))

    if missing_optional:
        print("\n⚠️  Optional tools not found (some features unavailable):")
        for tool, cmd in missing_optional:
            print(f"  - {tool:20s}  →  {cmd}")

    if missing_required:
        print("\n❌ Missing required tools — install before running FungalFlye:\n")
        for tool, cmd in missing_required:
            print(f"  {tool:20s}  →  {cmd}")
        print()
        raise SystemExit(1)

    # Check Python packages
    py_packages = {"Bio": "biopython", "pandas": "pandas", "matplotlib": "matplotlib"}
    missing_py = []
    for mod, pkg in py_packages.items():
        try:
            __import__(mod)
        except ImportError:
            missing_py.append(pkg)
    if missing_py:
        print(f"\n❌ Missing Python packages: {', '.join(missing_py)}")
        print(f"   Install with: pip install {' '.join(missing_py)}\n")
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

    def is_junk(kmer):
        if len(set(kmer)) == 1:
            return True  # pure homopolymer: TTTTTT, AAAAAA
        if len(set(kmer)) == 2 and kmer == kmer[:2] * (len(kmer) // 2):
            return True  # dinucleotide repeat: ATATAT, TATATA
        return False

    filtered = {k: v for k, v in kmer_counts.items() if not is_junk(k)}
    if not filtered:
        print("[fungalflye] Warning: could not find non-homopolymer motif — defaulting to TTAGGG")
        return "TTAGGG"

    most_common = sorted(filtered.items(), key=lambda x: x[1], reverse=True)[:10]
    print("\nTop candidate telomere kmers (homopolymers excluded):")
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
    check_dependencies()
    fasta = Path(fasta)
    print("\n[fungalflye] Starting QC...")
    outdir = fasta.parent / "fungalflye_qc"
    outdir.mkdir(parents=True, exist_ok=True)
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
