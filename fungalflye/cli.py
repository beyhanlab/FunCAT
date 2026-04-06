import typer
import subprocess
from pathlib import Path
from typing import Optional
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd
import matplotlib.pyplot as plt

from .assemble import run_assembly, READ_TYPE_CONFIGS
from .qc import run_qc
from .compare import run_snp_analysis
from .dotplot_run import run_dotplot

app = typer.Typer(help="FungalFlye — long-read fungal genome assembly toolkit")


# ------------------------------------------------
# helper runner
# ------------------------------------------------

def run(cmd):
    print(f"\n[fungalflye] Running: {cmd}\n")
    subprocess.run(cmd, shell=True, check=True)


# ------------------------------------------------
# read length analysis
# ------------------------------------------------

def analyze_reads(reads, outdir):

    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)

    lengths_file = outdir / "read_lengths.tsv"

    run(f"seqkit fx2tab -n -l {reads} > {lengths_file}")

    df = pd.read_csv(lengths_file, sep="\t", header=None)
    lengths = df[1]

    total_reads = len(lengths)
    total_bases = lengths.sum()

    sorted_lengths = sorted(lengths, reverse=True)
    cumsum = 0
    n50 = 0

    for L in sorted_lengths:
        cumsum += L
        if cumsum >= total_bases / 2:
            n50 = L
            break

    plt.figure(figsize=(6, 4))
    plt.hist(lengths, bins=100)
    plt.xlabel("Read length (bp)")
    plt.ylabel("Count")
    plt.title("Read length distribution")
    plt.tight_layout()
    plt.savefig(outdir / "read_length_histogram.png")
    plt.close()

    return total_reads, total_bases, n50, lengths


# ------------------------------------------------
# filtering preview
# ------------------------------------------------

def preview_filter(lengths, cutoff):
    kept = lengths[lengths >= cutoff]
    removed = lengths[lengths < cutoff]
    return len(kept), len(removed)


# ------------------------------------------------
# Assembly command
# ------------------------------------------------

@app.command()
def assemble(
    reads: str = typer.Argument(..., help="Path to input reads (FASTQ/FASTQ.gz)"),
    gsize: str = typer.Argument(..., help="Estimated genome size, e.g. 40m, 1.2g"),
    outdir: str = typer.Option("fungalflye_out", help="Output directory"),
    threads: int = typer.Option(8, help="Number of CPU threads"),
    min_read_len: int = typer.Option(0, help="Minimum read length filter (0 = off)"),
    downsample_cov: int = typer.Option(0, help="Downsample to this coverage depth (0 = off)"),
    min_contig_size: int = typer.Option(5000, help="Minimum contig size to keep after pruning"),
    read_type: str = typer.Option(
        "nano-hq",
        help="Read type: nano-hq | nano-raw | pacbio-hifi",
    ),
    ploidy: str = typer.Option(
        "haploid",
        help="Ploidy mode: haploid | diploid. Diploid enables --keep-haplotypes in Flye.",
    ),
    asm_coverage: int = typer.Option(
        60,
        help="Target coverage for Flye assembly (--asm-coverage)",
    ),
):
    """
    Run the full FungalFlye assembly pipeline.

    Read types:
      nano-hq      Nanopore HQ reads (R10.4+, Q20) — uses Medaka polishing
      nano-raw     Nanopore raw reads (R9.4, standard) — uses Medaka polishing
      pacbio-hifi  PacBio HiFi / CCS reads — uses Racon polishing

    Ploidy:
      haploid   Most filamentous fungi (Aspergillus, Neurospora, Fusarium...)
      diploid   Diploid/heterozygous species (Histoplasma, Candida, Cryptococcus...)
    """

    final = run_assembly(
        reads=reads,
        genome_size=gsize,
        outdir=outdir,
        threads=threads,
        min_read_len=min_read_len,
        downsample_cov=downsample_cov,
        min_contig_size=min_contig_size,
        read_type=read_type,
        ploidy=ploidy,
        asm_coverage=asm_coverage,
    )

    typer.echo(f"\n✅ Final assembly: {final}\n")


# ------------------------------------------------
# QC — telomere bug fixed: run_telomeres now
# exposed and correctly forwarded
# ------------------------------------------------

@app.command()
def qc(
    fasta: str = typer.Argument(..., help="Path to assembly FASTA"),
    telomere: Optional[str] = typer.Option(
        None, help="Telomere motif sequence (e.g. TTAGGG). Auto-discovered if omitted."
    ),
    run_telomeres: bool = typer.Option(
        True, help="Run telomere analysis (default: on)"
    ),
):
    """
    Run assembly QC: contig stats, length histogram, and telomere analysis.
    """
    run_qc(fasta, telomere=telomere, run_telomeres=run_telomeres)


# ------------------------------------------------
# SNP comparison
# ------------------------------------------------

@app.command()
def snps(
    reference: str = typer.Argument(..., help="Reference genome FASTA"),
    query: str = typer.Argument(..., help="Query genome FASTA"),
    outdir: str = typer.Option("fungalflye_snps", help="Output directory"),
):
    """
    Detect SNPs between two genome assemblies using NUCmer.
    """
    run_snp_analysis(reference, query, outdir)


# ------------------------------------------------
# Dotplot
# ------------------------------------------------

@app.command()
def dotplot(
    reference: str = typer.Argument(..., help="Reference genome FASTA"),
    query: str = typer.Argument(..., help="Query genome FASTA"),
    outdir: str = typer.Option("fungalflye_dotplots", help="Output directory"),
):
    """
    Generate a whole-genome dotplot between two assemblies.
    """
    run_dotplot(reference, query, outdir)


# ------------------------------------------------
# Batch comparisons — parallelised
# ------------------------------------------------

def _compare_pair(args):
    g1, g2, outdir = args
    pair_name = f"{g1.stem}_vs_{g2.stem}"
    pair_dir = outdir / pair_name
    pair_dir.mkdir(exist_ok=True)
    run_snp_analysis(g1, g2, pair_dir)
    run_dotplot(g1, g2, pair_dir)
    return pair_name


@app.command()
def compare_folder(
    folder: str = typer.Argument(..., help="Folder containing genome FASTA files"),
    outdir: str = typer.Option("fungalflye_comparisons", help="Output directory"),
    threads: int = typer.Option(4, help="Parallel comparison workers"),
):
    """
    Run SNP + dotplot comparisons for all genome pairs in a folder (parallelised).
    """

    folder = Path(folder)
    genomes = list(folder.glob("*.fasta")) + list(folder.glob("*.fa"))

    if len(genomes) < 2:
        typer.echo("Need at least two genomes in the folder.")
        raise typer.Exit()

    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)

    pairs = [
        (genomes[i], genomes[j], outdir)
        for i in range(len(genomes))
        for j in range(i + 1, len(genomes))
    ]

    typer.echo(f"\nFound {len(genomes)} genomes → {len(pairs)} pairs\n")

    with ProcessPoolExecutor(max_workers=threads) as pool:
        futures = {pool.submit(_compare_pair, p): p for p in pairs}
        for fut in as_completed(futures):
            try:
                name = fut.result()
                typer.echo(f"  ✅ {name}")
            except Exception as exc:
                typer.echo(f"  ❌ {futures[fut][0].stem}_vs_{futures[fut][1].stem}: {exc}")

    typer.echo("\n🎉 All comparisons finished\n")


if __name__ == "__main__":
    app()
