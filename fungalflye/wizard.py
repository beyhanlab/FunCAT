from pathlib import Path
import os
import time
import typer

from .assemble import run_assembly, READ_TYPE_CONFIGS, DEFAULT_ENHANCEMENTS, parse_genome_size
from .qc import run_qc, discover_telomere_motif
from .cli import analyze_reads, preview_filter
from .compare import run_snp_analysis
from .dotplot_run import run_dotplot

app = typer.Typer()

BANNER = r"""
============================================================

🧬🐉  FungalFlye  🐉🧬
Long-read fungal genome assembly toolkit
Version 0.3  |  Beyhan Lab  |  J. Craig Venter Institute

FungalFlye takes raw Nanopore or PacBio HiFi reads and
produces chromosome-scale fungal genome assemblies through
a fully automated pipeline. Starting from raw reads, it
performs adaptive Flye assembly, iterative Medaka polishing,
repeat-aware scaffolding to bridge fragmented contigs, and
telomere repeat motifs (auto-discovered or user-supplied).
Every assembly is scored for contig confidence and delivered
with a self-contained HTML report — no manual configuration
or bioinformatics expertise required.

Developed by: Jacob Durazo
Citation: Durazo J, et al. FungalFlye: a purpose-built
  long-read assembly toolkit for fungal genomes. (2025)
  [Manuscript in preparation]

------------------------------------------------------------
  Quick start:  select 1 for the full assembly pipeline
  Docs & code:  github.com/beyhanlab/FungalFlye
------------------------------------------------------------

============================================================
"""

ENHANCEMENT_MENU = {
    "adaptive_params": {
        "label": "Adaptive Flye parameters",
        "desc":  "Analyses your reads and auto-tunes Flye settings",
        "default": True,
        "time":  "< 1 min",
    },
    "iterative_polish": {
        "label": "Iterative Medaka polishing",
        "desc":  "Polishes up to 3 rounds, stops when assembly converges",
        "default": True,
        "time":  "+30–60 min",
    },
    "purge_dups": {
        "label": "Purge Duplicates  [diploid only]",
        "desc":  "Removes haplotig duplicates — dramatically improves diploid assemblies",
        "default": False,
        "time":  "+15 min",
    },
    "scaffolding": {
        "label": "Repeat-aware scaffolding",
        "desc":  "Uses long reads to bridge contig gaps at repeat boundaries",
        "default": False,
        "time":  "+20–40 min",
    },
    "telo_scaffolding": {
        "label": "Telomere-guided scaffolding",
        "desc":  "Attaches small telomeric fragments to chromosome ends",
        "default": False,
        "time":  "+10 min",
    },
    "confidence_scoring": {
        "label": "Contig confidence scoring",
        "desc":  "Flags suspicious contigs (collapsed repeats, contamination)",
        "default": True,
        "time":  "+5 min",
    },
}


def path_exists(p: str) -> str:
    p2 = os.path.expanduser(p.strip())
    if not Path(p2).exists():
        raise typer.BadParameter(f"Path not found: {p2}")
    return p2


def normalize_gsize(g: str) -> str:
    g = g.strip().lower()
    if g.isdigit():
        val = int(g)
        return f"{val}m" if val < 1000 else str(val)
    return g


def get_read_type() -> str:
    typer.echo("\nSelect read type:")
    for i, (key, cfg) in enumerate(READ_TYPE_CONFIGS.items(), 1):
        typer.echo(f"  {i}) {cfg['label']}  [{key}]")
    choice = typer.prompt("Enter choice", default="1")
    keys = list(READ_TYPE_CONFIGS.keys())
    try:
        return keys[int(choice.strip()) - 1]
    except (ValueError, IndexError):
        typer.echo("Invalid — defaulting to nano-hq")
        return "nano-hq"


def get_ploidy() -> str:
    typer.echo("\nSelect ploidy:")
    typer.echo("  1) Haploid  — most fungi (Aspergillus, Neurospora, Fusarium, Histoplasma...)")
    typer.echo("  2) Diploid  — heterozygous species (Candida albicans, Cryptococcus...)")
    choice = typer.prompt("Enter choice", default="1")
    return "diploid" if choice.strip() == "2" else "haploid"


def get_enhancements(ploidy: str) -> dict:
    """
    Interactive feature selection menu.
    Returns dict of {feature_key: bool}.
    """
    keys   = list(ENHANCEMENT_MENU.keys())
    active = {k: ENHANCEMENT_MENU[k]["default"] for k in keys}

    # Purge dups only available for diploid
    if ploidy == "haploid":
        active["purge_dups"] = False

    while True:
        typer.echo("\n" + "=" * 60)
        typer.echo("⚙️   Enhancement modules — choose what to run")
        typer.echo("=" * 60)

        for i, key in enumerate(keys, 1):
            m = ENHANCEMENT_MENU[key]
            status = "ON " if active[key] else "off"

            # Grey out purge_dups for haploid
            if key == "purge_dups" and ploidy == "haploid":
                typer.echo(
                    f"  -) [{status}]  {m['label']}  "
                    f"(not available for haploid)"
                )
            else:
                typer.echo(
                    f"  {i}) [{status}]  {m['label']}  "
                    f"({m['time']})  — {m['desc']}"
                )

        typer.echo("\n  A) Select all")
        typer.echo("  D) Use recommended defaults  [adaptive + iterative polish + confidence]")
        typer.echo("  C) Continue with current selection")

        choice = typer.prompt("\nEnter number to toggle, or A / D / C", default="C").strip().upper()

        if choice == "C":
            break

        if choice == "A":
            for k in keys:
                if not (k == "purge_dups" and ploidy == "haploid"):
                    active[k] = True
            typer.echo("✅ All modules selected")
            continue

        if choice == "D":
            active = {k: ENHANCEMENT_MENU[k]["default"] for k in keys}
            if ploidy == "haploid":
                active["purge_dups"] = False
            typer.echo("✅ Defaults restored")
            continue

        try:
            idx = int(choice) - 1
            key = keys[idx]
            if key == "purge_dups" and ploidy == "haploid":
                typer.echo("⚠️  Purge Duplicates is only available for diploid assemblies")
            else:
                active[key] = not active[key]
                state = "ON" if active[key] else "off"
                typer.echo(f"  → {ENHANCEMENT_MENU[key]['label']} set to {state}")
        except (ValueError, IndexError):
            typer.echo("Invalid choice — enter a number, A, D, or C")

    active_names = [ENHANCEMENT_MENU[k]["label"] for k, v in active.items() if v]
    typer.echo(f"\n✅ Running with: {', '.join(active_names) if active_names else 'no enhancements'}")

    return active


def get_telomere_setup():
    run_telomeres = typer.confirm("Run telomere analysis?", default=True)
    if not run_telomeres:
        return False, None, False
    typer.echo("\nDo you know the telomere motif?")
    typer.echo("  1) Yes — I'll enter it")
    typer.echo("  2) Auto-discover from assembly")
    choice = typer.prompt("Enter 1 or 2", default="2").strip()

    # If user accidentally typed a motif sequence instead of 1 or 2, catch it
    if len(choice) > 2 or (choice not in ("1", "2") and any(c in "ACGT" for c in choice.upper())):
        typer.echo(f"  → Detected motif sequence: {choice.upper()} — using it directly")
        return True, choice.upper(), False

    if choice == "1":
        motif = typer.prompt("Enter telomere motif (e.g. TTAGGG)").strip().upper()
        return True, motif, False
    return True, None, True


def abort():
    typer.echo("\n⚠️  Cancelled — returning to main menu.\n")


@app.command()
def wizard():

    from .qc import check_dependencies
    check_dependencies(assembly_mode=False)

    typer.echo(BANNER)
    typer.echo("Welcome to FungalFlye.\n")

    while True:

        typer.echo("\nSelect workflow mode:\n")
        typer.echo("  1) Full pipeline  (assembly → polish → QC)")
        typer.echo("  2) Assembly only")
        typer.echo("  3) QC only")
        typer.echo("  4) Compare genomes  (SNPs + dotplot)")
        typer.echo("  5) Exit\n")

        mode = typer.prompt("Enter choice", default="1")

        if mode == "5":
            typer.echo("\nGoodbye 🐉\n")
            return

        if mode not in ("1", "2", "3", "4"):
            typer.echo("Invalid choice — please enter 1–5.")
            continue

        # QC ONLY
        if mode == "3":
            fasta = typer.prompt("Path to assembly FASTA", value_proc=path_exists)
            run_telomeres, tel_motif, auto_tel = get_telomere_setup()
            if not typer.confirm("\nReady to run QC — continue?", default=True):
                abort()
                continue
            typer.echo("\n📊 Running QC...\n")
            if run_telomeres and auto_tel:
                tel_motif = discover_telomere_motif(fasta)
            run_qc(fasta, telomere=tel_motif, run_telomeres=run_telomeres)
            typer.echo("\n🎉 QC complete.\n")
            if not typer.confirm("Run another workflow?", default=False):
                return
            continue

        # COMPARE
        if mode == "4":
            while True:
                reference = typer.prompt("Reference genome", value_proc=path_exists)
                query     = typer.prompt("Query genome",     value_proc=path_exists)
                outdir    = typer.prompt("Output folder", default="fungalflye_compare")
                if not typer.confirm("\nReady to compare — continue?", default=True):
                    abort()
                    break
                if typer.confirm("Run SNP detection?", default=True):
                    run_snp_analysis(reference, query, outdir)
                if typer.confirm("Generate dotplot?", default=True):
                    run_dotplot(reference, query, outdir)
                if not typer.confirm("\nRun another comparison?", default=False):
                    break
            if not typer.confirm("Run another workflow?", default=False):
                return
            continue

        # ASSEMBLY / FULL PIPELINE
        reads  = typer.prompt("Path to raw reads", value_proc=path_exists)
        gsize  = normalize_gsize(typer.prompt("Genome size (e.g. 40m, 1.2g)", default="40m"))
        outdir = typer.prompt("Output folder", default="fungalflye_out")
        threads = typer.prompt("Threads", default=8, type=int)

        read_type = get_read_type()
        ploidy    = get_ploidy()

        # Enhancement selection
        enhancements = get_enhancements(ploidy)

        asm_coverage = 60  # handled automatically by adaptive params module

        Path(outdir).mkdir(exist_ok=True)

        if (Path(outdir) / "final.fasta").exists():
            typer.echo("\n⚠️  Existing assembly detected — pipeline will resume.\n")

        typer.echo("\n🧾 Assembly plan:")
        typer.echo(f"  Reads       : {reads}")
        typer.echo(f"  Genome size : {gsize}")
        typer.echo(f"  Read type   : {READ_TYPE_CONFIGS[read_type]['label']}")
        typer.echo(f"  Ploidy      : {ploidy}")
        typer.echo(f"  Threads     : {threads}")
        typer.echo(f"  Outdir      : {outdir}")

        if not typer.confirm("\nReady to analyze reads — continue?", default=True):
            abort()
            continue

        typer.echo("\n🔎 Step 1 — Read diagnostics\n")
        total_reads, total_bases, read_n50, lengths = analyze_reads(reads, outdir)
        typer.echo(f"Total reads : {total_reads:,}")
        typer.echo(f"Total bases : {total_bases:,}")
        typer.echo(f"Read N50    : {read_n50:,} bp")

        suggested_cutoff = int(read_n50 * 0.7)
        typer.echo(f"\nSuggested minimum read length: {suggested_cutoff} bp")

        min_read_len = 0
        if typer.confirm("Apply read length filter?", default=True):
            cutoff = typer.prompt("Minimum read length", default=suggested_cutoff)
            kept, removed = preview_filter(lengths, cutoff)
            genome_bp = parse_genome_size(gsize)
            kept_bases = lengths[lengths >= cutoff].sum()
            est_coverage = kept_bases / genome_bp if genome_bp > 0 else 0
            typer.echo(f"\nReads kept        : {kept:,}")
            typer.echo(f"Reads removed     : {removed:,}")
            typer.echo(f"Bases retained    : {kept_bases:,} bp")
            typer.echo(f"Estimated coverage: {est_coverage:.1f}x")
            if est_coverage < 20:
                typer.echo("  ⚠️  Coverage below 20x — assembly quality may be reduced")
            elif est_coverage > 100:
                typer.echo("  ℹ️  High coverage — consider downsampling to 60–80x for speed")
            else:
                typer.echo("  ✅ Coverage looks good for assembly")
            if not typer.confirm("\nConfirm these filter settings?", default=True):
                abort()
                continue
            min_read_len = cutoff

        downsample_cov  = typer.prompt("Downsample coverage? (0 = none)", default=0, type=int)
        min_contig_size = typer.prompt("Minimum contig size (bp)", default=5000, type=int)

        run_telomeres, tel_motif, auto_tel = get_telomere_setup()

        # Pass telomere motif into enhancements so telo_scaffolding can use it
        if tel_motif:
            enhancements["telo_motif"] = tel_motif

        typer.echo("\n" + "=" * 60)
        typer.echo("Ready to assemble:")
        typer.echo(f"  Read filter   : {min_read_len} bp min" if min_read_len else "  Read filter   : off")
        typer.echo(f"  Downsample    : {downsample_cov}x" if downsample_cov else "  Downsample    : off")
        typer.echo(f"  Contig cutoff : {min_contig_size} bp")
        typer.echo(f"  Ploidy        : {ploidy}")
        active_mods = [
    ENHANCEMENT_MENU[k]["label"]
    for k, v in enhancements.items()
    if v and k in ENHANCEMENT_MENU
]
        typer.echo(f"  Modules       : {', '.join(active_mods) if active_mods else 'none'}")
        typer.echo("=" * 60)

        if not typer.confirm("\nLaunch assembly?", default=True):
            abort()
            continue

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
            read_type=read_type,
            ploidy=ploidy,
            asm_coverage=asm_coverage,
            enhancements=enhancements,
        )

        if mode == "1" and typer.confirm("\nAssembly finished — run QC?", default=True):
            typer.echo("\n📊 Step 3 — Assembly QC\n")
            if run_telomeres and auto_tel:
                tel_motif = discover_telomere_motif(final_fasta)
            run_qc(final_fasta, telomere=tel_motif, run_telomeres=run_telomeres)

        elapsed = time.time() - start_time
        typer.echo("\n" + "=" * 60)
        typer.echo("🎉 Pipeline complete.")
        typer.echo("=" * 60)
        typer.echo(f"\nFinal assembly : {final_fasta}")
        typer.echo(f"Total runtime  : {elapsed / 60:.1f} minutes\n")
        typer.echo("Thank you for using FungalFlye 🐉\n")

        if not typer.confirm("Run another workflow?", default=False):
            return


if __name__ == "__main__":
    app()
