import subprocess
from pathlib import Path


def run(cmd):
    """Run external command safely."""
    print(f"\n[fungalflye] Running:\n{cmd}\n")

    result = subprocess.run(cmd, shell=True)

    if result.returncode != 0:
        raise RuntimeError(f"Command failed:\n{cmd}")


def run_snp_analysis(reference, query, outdir):
    """
    Run nucmer alignment and SNP detection between two genomes.
    Prints total substitution SNPs and saves full SNP table.
    """

    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)

    prefix = outdir / "comparison"

    delta = f"{prefix}.delta"
    snps_file = outdir / "snps.txt"

    print("\n🧬 Running genome alignment (NUCmer)...")

    run(f"nucmer --prefix={prefix} {reference} {query}")

    print("\n🧬 Detecting SNPs...")

    run(f"show-snps -ClrTH {delta} > {snps_file}")

    if not snps_file.exists():
        raise RuntimeError("SNP file was not generated")

    # ------------------------------------------------
    # Count substitution SNPs
    # ------------------------------------------------

    snp_counts = {}
    total_snps = 0

    with open(snps_file) as f:
        for line in f:

            parts = line.strip().split()

            if len(parts) < 3:
                continue

            ref = parts[1]
            qry = parts[2]

            # substitutions only
            if ref != "." and qry != ".":

                total_snps += 1

                change = f"{ref}->{qry}"

                if change not in snp_counts:
                    snp_counts[change] = 0

                snp_counts[change] += 1

    print("\n" + "=" * 60)
    print("🧬 SNP SUMMARY")
    print("=" * 60)

    print(f"\nTotal substitution SNPs: {total_snps}\n")

    for change, count in sorted(snp_counts.items()):
        print(f"{change}: {count}")

    print(f"\n📄 Full SNP table saved to:\n{snps_file}\n")


# ------------------------------------------------
# Dotplot wrapper
# ------------------------------------------------

def run_dotplot(reference, query, outdir):
    """
    Run genome dotplot generation.
    """

    from .dotplot_run import run_dotplot as _run_dotplot

    _run_dotplot(reference, query, outdir)