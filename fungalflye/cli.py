import typer
import subprocess
from pathlib import Path
from typing import Optional
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .assemble import run_assembly, READ_TYPE_CONFIGS
from .qc import run_qc, scan_telomeres, discover_telomere_motif
from .compare import run_snp_analysis
from .dotplot_run import run_dotplot
from .report import generate_report

def _check_dependencies():
    """Check all required tools are installed and warn about optional ones."""
    import shutil

    required = {
        "flye":      "conda install -c bioconda flye",
        "minimap2":  "conda install -c bioconda minimap2",
        "seqkit":    "conda install -c bioconda seqkit",
        "samtools":  "conda install -c bioconda samtools",
        "filtlong":  "conda install -c bioconda filtlong",
    }
    optional = {
        "medaka":       "pip install medaka",
        "racon":        "conda install -c bioconda racon",
        "purge_dups":   "conda install -c bioconda purge_dups",
        "nucmer":       "conda install -c bioconda mummer",
        "polypolish":   "conda install -c bioconda polypolish",
        "polypolish_insert_filter": "conda install -c bioconda polypolish",
        "pilon":        "conda install -c bioconda pilon",
        "bwa":          "conda install -c bioconda bwa",
    }

    missing_required = []
    missing_optional = []

    for tool, install_cmd in required.items():
        if shutil.which(tool) is None:
            missing_required.append((tool, install_cmd))

    for tool, install_cmd in optional.items():
        if shutil.which(tool) is None:
            missing_optional.append((tool, install_cmd))

    if missing_required:
        typer.echo("\n❌  Missing required dependencies:\n", err=True)
        for tool, cmd in missing_required:
            typer.echo(f"   {tool:15s}  →  {cmd}", err=True)
        typer.echo("\nInstall the above tools and re-run FunCAT.\n", err=True)
        raise typer.Exit(1)

    if missing_optional:
        typer.echo("\nℹ️   Some optional modules are not installed (core assembly works fine):")
        for tool, cmd in missing_optional:
            typer.echo(f"   {tool:30s}  →  {cmd}")
        typer.echo("   Install any of the above to unlock those specific features.\n")


def _default_wizard():
    """Launch the wizard when funcat is run with no subcommand."""
    _check_dependencies()
    from .wizard import wizard
    wizard()

app = typer.Typer(
    help=(
        "FunCAT — Fungal Chromosome Assembly Tool\n\n"
        "Developed by Jacob Durazo, Beyhan Lab, J. Craig Venter Institute\n\n"
        "GETTING STARTED:\n\n"
        "  Run the interactive wizard (recommended):\n\n"
        "    funcat\n\n"
        "  This guides you step-by-step through assembly, polishing,\n"
        "  scaffolding, QC, and HTML report generation.\n\n"
        "COMMON COMMANDS:\n\n"
        "  funcat              Launch interactive wizard (full pipeline)\n"
        "  funcat qc           Run QC + HTML report on any assembly FASTA\n"
        "  funcat report       Generate a standalone HTML report\n"
        "  funcat telo-scaffold  Attach telomeric fragments to chromosome ends\n\n"
        "ILLUMINA POLISHING:\n\n"
        "  Add Illumina reads for dramatically improved base accuracy:\n\n"
        "    funcat assemble long_reads.fastq.gz 40m \\\n"
        "      --illumina-r1 illumina_R1.fastq.gz \\\n"
        "      --illumina-r2 illumina_R2.fastq.gz\n\n"
        "Citation: Durazo J, et al. FunCAT (2025) [Manuscript in preparation]"
    ),
    invoke_without_command=True,
)

@app.callback()
def main(ctx: typer.Context):
    if ctx.invoked_subcommand is None:
        _default_wizard()


def run(cmd):
    print(f"\n[funcat] Running: {cmd}\n")
    subprocess.run(cmd, shell=True, check=True)


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
    cumsum, n50 = 0, 0
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


def preview_filter(lengths, cutoff):
    kept    = lengths[lengths >= cutoff]
    removed = lengths[lengths  < cutoff]
    return len(kept), len(removed)


@app.command()
def assemble(
    reads:          str = typer.Argument(..., help="Path to input reads (FASTQ/FASTQ.gz)"),
    gsize:          str = typer.Argument(..., help="Estimated genome size, e.g. 40m, 1.2g"),
    outdir:         str = typer.Option("funcat_out",  help="Output directory"),
    threads:        int = typer.Option(8,                 help="CPU threads"),
    min_read_len:   int = typer.Option(0,                 help="Min read length filter (0=off)"),
    downsample_cov: int = typer.Option(0,                 help="Downsample coverage (0=off)"),
    min_contig_size:int = typer.Option(5000,              help="Min contig size after pruning"),
    read_type:      str = typer.Option("nano-hq",         help="nano-hq | nano-raw | pacbio-hifi"),
    ploidy:         str = typer.Option("haploid",         help="haploid | diploid"),
    asm_coverage:   int = typer.Option(60,                help="Flye --asm-coverage"),
    illumina_r1:    Optional[str] = typer.Option(None,    help="Illumina R1 reads for polishing (optional)"),
    illumina_r2:    Optional[str] = typer.Option(None,    help="Illumina R2 reads for polishing (optional)"),
    illumina_polisher: str = typer.Option("polypolish",   help="polypolish | pilon"),
):
    """
    Run the full FunCAT assembly pipeline (non-interactive).
    For guided setup, just run: funcat
    
    NEW: Illumina polishing support!
    Use --illumina-r1 and --illumina-r2 to provide Illumina reads
    for polishing after Medaka. This dramatically improves base
    accuracy and reduces internal stop codons.
    """
    _check_dependencies()
    final = run_assembly(
        reads=reads, genome_size=gsize, outdir=outdir, threads=threads,
        min_read_len=min_read_len, downsample_cov=downsample_cov,
        min_contig_size=min_contig_size, read_type=read_type,
        ploidy=ploidy, asm_coverage=asm_coverage,
        illumina_r1=illumina_r1, illumina_r2=illumina_r2,
        illumina_polisher=illumina_polisher,
    )
    typer.echo(f"\n✅ Final assembly: {final}\n")


@app.command()
def qc(
    fasta:         str           = typer.Argument(..., help="Path to assembly FASTA"),
    telomere:      Optional[str] = typer.Option(None,  help="Telomere motif (auto if omitted)"),
    run_telomeres: bool          = typer.Option(True,  help="Run telomere analysis"),
    html_report:   bool          = typer.Option(True,  help="Generate HTML report"),
    name:          str           = typer.Option("",    help="Assembly name for report"),
):
    """
    Run QC on any assembly FASTA: stats, contig histogram,
    telomere completeness scan, and self-contained HTML report.

    Example:
      funcat qc assembly.fasta --telomere TTAGGG --name MyStrain
    """
    run_qc(
        fasta=fasta, telomere=telomere, run_telomeres=run_telomeres,
        report=html_report,
        run_metadata={"assembly_name": name or Path(fasta).stem},
    )


@app.command()
def report(
    fasta:      str           = typer.Argument(..., help="Path to assembly FASTA"),
    outdir:     str           = typer.Option("",   help="Output dir (defaults to FASTA dir)"),
    name:       str           = typer.Option("",   help="Assembly name"),
    confidence: Optional[str] = typer.Option(None, help="Path to contig_confidence.tsv"),
    telomere:   Optional[str] = typer.Option(None, help="Telomere motif (auto if omitted)"),
):
    """
    Generate a standalone HTML report for any assembly FASTA.
    Works on any FASTA — does not need to be from a full funcat run.
    """
    fasta_path = Path(fasta)
    out_path   = Path(outdir) if outdir else fasta_path.parent

    telo_df = None
    if telomere:
        telo_df = scan_telomeres(str(fasta_path), telomere)
    elif typer.confirm("Run telomere auto-discovery?", default=True):
        motif   = discover_telomere_motif(str(fasta_path))
        telo_df = scan_telomeres(str(fasta_path), motif)

    html_path = generate_report(
        fasta=fasta_path,
        outdir=out_path,
        run_metadata={"assembly_name": name or fasta_path.stem},
        telo_df=telo_df,
        confidence_tsv=confidence,
    )
    typer.echo(f"\n✅ Report: {html_path}\n")


@app.command()
def telo_scaffold(
    fasta:    str = typer.Argument(..., help="Path to assembly FASTA"),
    reads:    str = typer.Argument(..., help="Path to reads FASTQ"),
    outdir:   str = typer.Option("funcat_telo", help="Output directory"),
    motif:    str = typer.Option("TTAGGG", help="Telomere motif"),
    threads:  int = typer.Option(8, help="CPU threads"),
    support:  int = typer.Option(5, help="Min supporting reads for a bridge"),
):
    """
    Run telomere-guided scaffolding on an existing assembly.
    Attaches small telomeric fragments to uncapped chromosome ends.
    """
    from .scaffold import run_telomere_scaffolding
    result = run_telomere_scaffolding(
        assembly=fasta, reads=reads, outdir=outdir, threads=threads,
        minimap2_preset="map-ont", telomere_motif=motif, min_support=support,
    )
    typer.echo(f"\n✅ Telomere-scaffolded assembly: {result}\n")


@app.command()
def snps(
    reference: str = typer.Argument(..., help="Reference genome FASTA"),
    query:     str = typer.Argument(..., help="Query genome FASTA"),
    outdir:    str = typer.Option("funcat_snps", help="Output directory"),
):
    """Detect SNPs between two genome assemblies using NUCmer."""
    run_snp_analysis(reference, query, outdir)


@app.command()
def dotplot(
    reference: str = typer.Argument(..., help="Reference genome FASTA"),
    query:     str = typer.Argument(..., help="Query genome FASTA"),
    outdir:    str = typer.Option("funcat_dotplots", help="Output directory"),
):
    """Generate a whole-genome dotplot between two assemblies."""
    run_dotplot(reference, query, outdir)


def _compare_pair(args):
    g1, g2, outdir = args
    pair_name = f"{g1.stem}_vs_{g2.stem}"
    pair_dir  = outdir / pair_name
    pair_dir.mkdir(exist_ok=True)
    run_snp_analysis(g1, g2, pair_dir)
    run_dotplot(g1, g2, pair_dir)
    return pair_name


@app.command()
def compare_folder(
    folder:  str = typer.Argument(..., help="Folder of genome FASTAs"),
    outdir:  str = typer.Option("funcat_comparisons", help="Output directory"),
    threads: int = typer.Option(4, help="Parallel workers"),
):
    """Run SNP + dotplot comparisons for all genome pairs in a folder (parallelised)."""
    folder  = Path(folder)
    genomes = list(folder.glob("*.fasta")) + list(folder.glob("*.fa"))
    if len(genomes) < 2:
        typer.echo("Need at least two genomes.")
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
                typer.echo(f"  ✅ {fut.result()}")
            except Exception as exc:
                g1, g2, _ = futures[fut]
                typer.echo(f"  ❌ {g1.stem}_vs_{g2.stem}: {exc}")
    typer.echo("\n🎉 All comparisons finished\n")


if __name__ == "__main__":
    app()
