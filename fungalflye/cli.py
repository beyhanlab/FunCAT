import typer
import subprocess
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

from .assemble import run_assembly
from .qc import run_qc
from .compare import run_snp_analysis
from .dotplot_run import run_dotplot

app = typer.Typer(help="FungalFlye Long-read fungal genome assembly toolkit")


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

    # calculate N50
    sorted_lengths = sorted(lengths, reverse=True)
    cumsum = 0
    n50 = 0

    for L in sorted_lengths:
        cumsum += L
        if cumsum >= total_bases / 2:
            n50 = L
            break

    # histogram
    plt.figure(figsize=(6, 4))
    plt.hist(lengths, bins=100)
    plt.xlabel("Read length (bp)")
    plt.ylabel("Count")
    plt.title("Read length distribution")
    plt.tight_layout()
    plt.savefig(outdir / "read_length_histogram.png")

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
    reads: str,
    gsize: str,
    outdir: str = "fungalflye_out",
    threads: int = 8,
    min_read_len: int = 0,
    downsample_cov: int = 0,
    min_contig_size: int = 20000
):
    """
    Run full FungalFlye assembly pipeline
    """

    final = run_assembly(
        reads=reads,
        genome_size=gsize,
        outdir=outdir,
        threads=threads,
        min_read_len=min_read_len,
        downsample_cov=downsample_cov,
        min_contig_size=min_contig_size
    )

    typer.echo(f"\n✅ Final assembly: {final}\n")


# ------------------------------------------------
# QC
# ------------------------------------------------
@app.command()
def qc(
    fasta: str,
    telomere: str | None = None
):
    """
    Run assembly QC and telomere analysis
    """
    run_qc(fasta, telomere)


# ------------------------------------------------
# SNP comparison
# ------------------------------------------------
@app.command()
def snps(
    reference: str,
    query: str,
    outdir: str = "fungalflye_snps"
):
    """
    Detect SNPs between two genomes
    """

    run_snp_analysis(reference, query, outdir)


# ------------------------------------------------
# Dotplot
# ------------------------------------------------
@app.command()
def dotplot(
    reference: str,
    query: str,
    outdir: str = "fungalflye_dotplots"
):
    """
    Generate genome dotplot
    """

    run_dotplot(reference, query, outdir)


# ------------------------------------------------
# Batch comparisons
# ------------------------------------------------
@app.command()
def compare_folder(
    folder: str,
    outdir: str = "fungalflye_comparisons"
):
    """
    Run SNP + dotplot comparisons for all genome pairs in a folder
    """

    folder = Path(folder)
    genomes = list(folder.glob("*.fasta")) + list(folder.glob("*.fa"))

    if len(genomes) < 2:
        typer.echo("Need at least two genomes")
        raise typer.Exit()

    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)

    typer.echo(f"\nFound {len(genomes)} genomes\n")

    for i in range(len(genomes)):
        for j in range(i + 1, len(genomes)):

            g1 = genomes[i]
            g2 = genomes[j]

            pair_name = f"{g1.stem}_vs_{g2.stem}"

            typer.echo(f"\n🧬 Comparing {pair_name}\n")

            pair_dir = outdir / pair_name
            pair_dir.mkdir(exist_ok=True)

            run_snp_analysis(g1, g2, pair_dir)
            run_dotplot(g1, g2, pair_dir)

    typer.echo("\n🎉 All comparisons finished\n")


# ------------------------------------------------
# INTERACTIVE MODE
# ------------------------------------------------
@app.command()
def interactive():

    typer.echo("\n🧬 Welcome to FungalFlye")
    typer.echo("Long-read fungal genome assembly assistant\n")

    reads = typer.prompt("Enter path to Nanopore reads")
    gsize_input = typer.prompt("Estimated genome size (e.g. 40m)")
    threads = typer.prompt("Threads", default=8)
    outdir = typer.prompt("Output directory", default="fungalflye_out")

    Path(outdir).mkdir(exist_ok=True)

    typer.echo("\n🔎 Analyzing read lengths...\n")

    total_reads, total_bases, read_n50, lengths = analyze_reads(reads, outdir)

    gsize_numeric = gsize_input.lower().replace("m", "000000")
    gsize_numeric = int(gsize_numeric)

    coverage = total_bases / gsize_numeric

    typer.echo(f"Total reads: {total_reads:,}")
    typer.echo(f"Total bases: {total_bases:,}")
    typer.echo(f"Read N50: {read_n50:,} bp")
    typer.echo(f"Estimated coverage: {coverage:.1f}×")

    suggested_cutoff = int(read_n50 * 0.7)

    typer.echo(f"\nSuggested minimum read length: {suggested_cutoff} bp")

    apply_filter = typer.confirm("Apply read filtering?", default=True)

    filtered_reads = reads

    if apply_filter:

        cutoff = typer.prompt(
            "Minimum read length",
            default=suggested_cutoff
        )

        kept, removed = preview_filter(lengths, cutoff)

        typer.echo(f"\nReads kept: {kept:,}")
        typer.echo(f"Reads removed: {removed:,}")

        confirm = typer.confirm("Continue with filtering?", default=True)

        if confirm:

            filtered_reads = f"{outdir}/filtered.fastq"

            run(
                f"seqkit seq -m {cutoff} {reads} > {filtered_reads}"
            )

    typer.echo("\n🚀 Starting assembly pipeline\n")

    assemble(
        reads=filtered_reads,
        gsize=gsize_input,
        outdir=outdir,
        threads=threads
    )

    typer.echo("\n🎉 Pipeline complete!\n")


if __name__ == "__main__":
    app()
