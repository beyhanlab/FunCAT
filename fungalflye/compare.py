import subprocess
from pathlib import Path


def run(cmd):
    print(f"\n[fungalflye] Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)


def run_snp_analysis(reference, query, outdir):

    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)

    prefix = outdir / "comparison"

    delta = f"{prefix}.delta"
    snps = outdir / "snps.txt"

    run(f"nucmer -p {prefix} {reference} {query}")
    run(f"show-snps -Clr {delta} > {snps}")

    print(f"\nSNP results saved to {snps}")