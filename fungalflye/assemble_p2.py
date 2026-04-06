
def prune_small_contigs(fasta, out_fasta, min_size=5000):
    print(f"\n[fungalflye] Removing contigs < {min_size} bp\n")
    fasta, out_fasta = Path(fasta), Path(out_fasta)
    kept, removed_count = [], 0
    for r in SeqIO.parse(str(fasta), "fasta"):
        if len(r.seq) < min_size:
            removed_count += 1
        else:
            kept.append(r)
    SeqIO.write(kept, str(out_fasta), "fasta")
    print(f"[fungalflye] Removed {removed_count} small contigs")
    return removed_count


def _warn_if_at_rich(fasta, sample=3):
    try:
        total, gc = 0, 0
        for i, r in enumerate(SeqIO.parse(str(fasta), "fasta")):
            if i >= sample:
                break
            s = str(r.seq).upper()
            total += len(s)
            gc += s.count("G") + s.count("C")
        if total > 0:
            pct = 100 * gc / total
            if pct < 35:
                print(f"\n⚠️  Low GC content ({pct:.1f}%) — AT-rich genome. Consider --asm-coverage 80.\n")
    except Exception:
        pass


def _separate_mito(polished_fasta, assembly_info, outdir):
    mito_ids = set()
    try:
        with open(assembly_info) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 4:
                    continue
                if parts[3].strip() == "Y" and 20_000 <= int(parts[1]) <= 100_000:
                    mito_ids.add(parts[0])
    except Exception:
        return
    if not mito_ids:
        return
    nuclear, mito = [], []
    for r in SeqIO.parse(str(polished_fasta), "fasta"):
        (mito if r.id in mito_ids else nuclear).append(r)
    if mito:
        SeqIO.write(nuclear, str(outdir / "nuclear.fasta"), "fasta")
        SeqIO.write(mito,    str(outdir / "mitochondrial.fasta"), "fasta")
        print(f"\n[fungalflye] Separated {len(mito)} mitochondrial contig(s)")
