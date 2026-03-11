import subprocess
from pathlib import Path
import shutil as sh


def run(cmd):
    print(f"\n[fungalflye] Running: {cmd}\n")
    subprocess.run(cmd, shell=True, check=True)


def check_dependencies():

    tools = ["nucmer", "show-snps", "show-coords", "delta-filter", "mummerplot"]

    missing = []

    for t in tools:
        if sh.which(t) is None:
            missing.append(t)

    if missing:
        print("\n❌ Missing required comparison tools:")
        for m in missing:
            print(f"  - {m}")

        print("\nInstall with:")
        print("  conda install -c bioconda mummer\n")

        raise SystemExit(1)


def label_from_path(path):
    return Path(path).stem


def run_snp_analysis(reference, query, outdir):

    check_dependencies()

    reference = Path(reference)
    query = Path(query)
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)

    ref_label = label_from_path(reference)
    qry_label = label_from_path(query)

    prefix = outdir / f"{ref_label}_vs_{qry_label}"

    delta_file = f"{prefix}.delta"
    snps_file = outdir / f"{ref_label}_vs_{qry_label}_snps.txt"
    positions_file = outdir / f"{ref_label}_vs_{qry_label}_snp_positions.txt"
    summary_file = outdir / f"{ref_label}_vs_{qry_label}_snp_summary.tsv"

    print("\n[fungalflye] Running nucmer for SNP analysis")
    run(f"nucmer --prefix={prefix} {reference} {query}")

    print("\n[fungalflye] Calling SNPs")
    run(f"show-snps -Clr {delta_file} > {snps_file}")

    snp_counts = {}
    positions = []

    with open(snps_file) as fh:
        seen_header = False
        for line in fh:
            line = line.strip()
            if not seen_header:
                if line.startswith("===="):
                    seen_header = True
                continue

            parts = line.split()
            if len(parts) < 3:
                continue

            position = parts[0]
            ref_base = parts[1]
            qry_base = parts[2]

            if ref_base != "." and qry_base != ".":
                positions.append((int(position), ref_base, qry_base))
                change = f"{ref_base}->{qry_base}"
                snp_counts[change] = snp_counts.get(change, 0) + 1

    positions.sort()

    with open(positions_file, "w") as out:
        out.write("Position\tReference_Base\tQuery_Base\n")
        for position, ref_base, qry_base in positions:
            out.write(f"{position}\t{ref_base}\t{qry_base}\n")

    with open(summary_file, "w") as out:
        out.write("Change\tCount\n")
        for change, count in sorted(snp_counts.items()):
            out.write(f"{change}\t{count}\n")

    total_snps = sum(snp_counts.values())

    print("\n" + "=" * 60)
    print("🧬 SNP Analysis Complete")
    print("=" * 60)
    print(f"Reference: {reference}")
    print(f"Query:     {query}")
    print(f"Total SNPs: {total_snps}")
    print(f"Positions: {positions_file}")
    print(f"Summary:   {summary_file}")
    print("=" * 60 + "\n")

    return {
        "delta": str(delta_file),
        "snps": str(snps_file),
        "positions": str(positions_file),
        "summary": str(summary_file),
    }


def run_dotplot(reference, query, outdir):

    check_dependencies()

    reference = Path(reference)
    query = Path(query)
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)

    ref_label = label_from_path(reference)
    qry_label = label_from_path(query)

    prefix = outdir / f"{ref_label}_vs_{qry_label}"

    delta_file = f"{prefix}.delta"
    filtered_delta = f"{prefix}.1delta"

    print("\n[fungalflye] Running nucmer for dotplot")
    run(f"nucmer --prefix={prefix} {reference} {query}")

    print("\n[fungalflye] Filtering alignments for clean dotplot")
    run(f"delta-filter -1 {delta_file} > {filtered_delta}")

    print("\n[fungalflye] Generating dotplot")
    run(f"mummerplot --png --layout -p {prefix} {filtered_delta}")

    png_file = f"{prefix}.png"

    print("\n" + "=" * 60)
    print("🧬 Dotplot Complete")
    print("=" * 60)
    print(f"Dotplot: {png_file}")
    print("=" * 60 + "\n")

    return {"dotplot_png": png_file}