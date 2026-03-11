from pathlib import Path
import os
import time
import typer

from .assemble import run_assembly
from .qc import run_qc, discover_telomere_motif
from .cli import analyze_reads, preview_filter
from .compare import run_snp_analysis, run_dotplot

app = typer.Typer()

BANNER = r"""
============================================================

🧬🐉  FungalFlye  🐉🧬
Long-read fungal genome assembly toolkit

Chromosome-scale assemblies from Nanopore reads
with automated QC, telomere validation,
SNP detection, and genome dotplot support

============================================================
"""


def pause(msg="Press Enter to continue..."):
    typer.echo("")
    typer.prompt(msg, default="", show_default=False)
    typer.echo("")


def path_exists(p: str) -> str:
    p2 = os.path.expanduser(p.strip())
    if not Path(p2).exists():
        raise typer.BadParameter(f"Path not found: {p2}")
    return p2


def normalize_gsize(g: str) -> str:
    g = g.strip().lower()
    if g.isdigit():
        val = int(g)
        if val < 1000:
            return f"{val}m"
        return str(val)
    return g


def get_telomere_setup():

    run_telomeres = typer.confirm("Run telomere analysis?", default=True)

    if not run_telomeres:
        return False, None, False

    typer.echo("\nDo you know the telomere motif?")
    typer.echo("  1) Yes")
    typer.echo("  2) Auto discover")

    choice = typer.prompt("Enter 1 or 2", default="2")

    if choice.strip() == "1":
        motif = typer.prompt("Enter telomere motif sequence").strip().upper()
        return True, motif, False

    return True, None, True


@app.command()
def wizard():

    typer.echo(BANNER)
    typer.echo("Welcome to FungalFlye.\n")

    while True:

        typer.echo("\nSelect workflow mode:\n")
        typer.echo("  1) Full pipeline")
        typer.echo("     Reads → Assembly → QC → Telomeres")
        typer.echo("  2) Assembly only")
        typer.echo("     Reads → Assembly")
        typer.echo("  3) QC only")
        typer.echo("     Existing FASTA → QC report")
        typer.echo("  4) Compare genomes")
        typer.echo("     SNP detection + dotplot")
        typer.echo("  5) Exit\n")

        mode = typer.prompt("Enter choice", default="1")

        # ------------------------------------------------
        # EXIT
        # ------------------------------------------------

        if mode == "5":
            return

        # ------------------------------------------------
        # QC ONLY
        # ------------------------------------------------

        if mode == "3":

            fasta = typer.prompt(
                "Path to assembly FASTA",
                value_proc=path_exists
            )

            run_telomeres, tel_motif, auto_tel = get_telomere_setup()

            typer.echo("\n📊 Running QC...\n")

            if run_telomeres and auto_tel:
                tel_motif = discover_telomere_motif(fasta)

            run_qc(
                fasta,
                telomere=tel_motif,
                run_telomeres=run_telomeres
            )

            typer.echo("\n🎉 QC complete.\n")
            return

        # ------------------------------------------------
        # COMPARE GENOMES
        # ------------------------------------------------

     if mode == "4":

    while True:

        reference = typer.prompt(
            "Reference genome",
            value_proc=path_exists
        )

        query = typer.prompt(
            "Query genome",
            value_proc=path_exists
        )

        outdir = typer.prompt(
            "Output folder",
            default="fungalflye_compare"
        )

        if typer.confirm("Run SNP detection?", default=True):
            run_snp_analysis(reference, query, outdir)

        if typer.confirm("Generate dotplot?", default=True):
            run_dotplot(reference, query, outdir)

        again = typer.prompt(
            "\nRun another comparison?\n1 Yes\n2 No",
            default="2"
        )

        if again != "1":
            break

    return

        # ------------------------------------------------
        # ASSEMBLY / FULL
        # ------------------------------------------------

        reads = typer.prompt(
            "Path to raw reads",
            value_proc=path_exists
        )

        gsize = typer.prompt("Genome size (e.g., 40m)", default="40m")
        gsize = normalize_gsize(gsize)

        outdir = typer.prompt("Output folder", default="fungalflye_out")
        threads = typer.prompt("Threads", default=8, type=int)

        Path(outdir).mkdir(exist_ok=True)

        if (Path(outdir) / "final.fasta").exists():
            typer.echo("\n⚠️ Existing assembly detected in output folder.")
            typer.echo("Pipeline will automatically resume from the last completed step.\n")

        typer.echo("\n🧾 Plan:")
        typer.echo(f"Reads: {reads}")
        typer.echo(f"Genome size: {gsize}")
        typer.echo(f"Outdir: {outdir}")
        typer.echo(f"Threads: {threads}")

        pause("Press Enter to analyze reads")

        typer.echo("\n🔎 Step 1 — Read diagnostics\n")

        total_reads, total_bases, read_n50, lengths = analyze_reads(reads, outdir)

        typer.echo(f"Total reads: {total_reads:,}")
        typer.echo(f"Total bases: {total_bases:,}")
        typer.echo(f"Read N50: {read_n50:,} bp")

        suggested_cutoff = int(read_n50 * 0.7)

        typer.echo(f"\nSuggested minimum read length: {suggested_cutoff} bp")

        apply_filter = typer.confirm("Apply read filtering?", default=True)

        min_read_len = 0

        if apply_filter:
            cutoff = typer.prompt("Minimum read length", default=suggested_cutoff)
            kept, removed = preview_filter(lengths, cutoff)

            typer.echo(f"\nReads kept: {kept:,}")
            typer.echo(f"Reads removed: {removed:,}")

            if typer.confirm("Continue?", default=True):
                min_read_len = cutoff

        downsample_cov = typer.prompt(
            "Downsample coverage? (0 = none)",
            default=0,
            type=int
        )

        min_contig_size = typer.prompt(
            "Minimum contig size to keep after pruning (bp)",
            default=20000,
            type=int
        )

        run_telomeres, tel_motif, auto_tel = get_telomere_setup()

        pause("Ready to begin assembly 🚀")

        start_time = time.time()

        typer.echo("\n🧬 Step 2 — Genome assembly\n")

        final_fasta = run_assembly(
            reads=reads,
            genome_size=gsize,
            outdir=outdir,
            threads=threads,
            min_read_len=min_read_len,
            downsample_cov=downsample_cov,
            min_contig_size=min_contig_size,
        )

        if mode == "1":

            pause("Assembly finished. Press Enter to run QC")

            typer.echo("\n📊 Step 3 — Assembly QC\n")

            if run_telomeres and auto_tel:
                tel_motif = discover_telomere_motif(final_fasta)

            run_qc(
                final_fasta,
                telomere=tel_motif,
                run_telomeres=run_telomeres
            )

        elapsed = time.time() - start_time

        typer.echo("\n" + "=" * 60)
        typer.echo("🎉 Pipeline complete.")
        typer.echo("=" * 60)
        typer.echo(f"\nFinal assembly: {final_fasta}\n")
        typer.echo(f"Total runtime: {elapsed/60:.1f} minutes\n")
        typer.echo("Thank you for using FungalFlye 🐉\n")

        return


if __name__ == "__main__":
    app()