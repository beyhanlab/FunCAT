"""
fungalflye.scaffold
~~~~~~~~~~~~~~~~~~~
Module 5 — Repeat-aware scaffolding

Uses long reads to bridge contig gaps at repeat boundaries.
Strategy:
  1. Map all reads to the assembly with minimap2
  2. Find reads that map to the END of one contig and the START of another
     (i.e. they span a gap between two contigs)
  3. Count supporting reads for each possible contig-pair join
  4. Accept joins with >= min_support reads and build scaffolds
  5. Write scaffolded FASTA with gaps filled by Ns

This is pure Python + minimap2 — no extra tools required.
"""

import subprocess
import shutil
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


GAP_SEQUENCE = "N" * 100   # 100 Ns between joined contigs


def run(cmd):
    print(f"\n[fungalflye] Running: {cmd}\n")
    subprocess.run(cmd, shell=True, check=True)


def _parse_paf(paf_path, contig_lengths, end_window=2000):
    """
    Parse PAF alignments and find reads that anchor to the
    end-region of one contig and the start-region of another.

    end_window: how many bp from each end counts as a "terminal" region.

    Returns: dict of (contig_a, side_a, contig_b, side_b) -> support_count
    """

    # For each read, collect all contig terminal hits
    read_hits = defaultdict(list)   # read_id -> list of (contig, side)

    with open(paf_path) as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 12:
                continue

            read_id = cols[0]
            t_name  = cols[5]
            t_start = int(cols[7])
            t_end   = int(cols[8])
            t_len   = int(cols[6])
            mapq    = int(cols[11])

            if mapq < 10:
                continue

            t_len_actual = contig_lengths.get(t_name, t_len)

            at_start = t_start <= end_window
            at_end   = t_end   >= t_len_actual - end_window

            if at_start and not at_end:
                read_hits[read_id].append((t_name, "start"))
            elif at_end and not at_start:
                read_hits[read_id].append((t_name, "end"))
            # reads entirely internal or spanning whole contig — skip

    # Count contig-pair joins supported by reads
    join_support = defaultdict(int)

    for read_id, hits in read_hits.items():
        if len(hits) < 2:
            continue

        # Deduplicate hits for this read
        unique_hits = list(set(hits))

        for i in range(len(unique_hits)):
            for j in range(i + 1, len(unique_hits)):
                a_contig, a_side = unique_hits[i]
                b_contig, b_side = unique_hits[j]

                if a_contig == b_contig:
                    continue

                # Canonical key: sort so (A,B) and (B,A) are the same
                key = tuple(sorted([
                    (a_contig, a_side),
                    (b_contig, b_side)
                ]))
                join_support[key] += 1

    return join_support


def _build_scaffolds(records_dict, joins, min_support):
    """
    Given contig records and supported joins, greedily build scaffolds.

    Strategy: treat as a path-building problem.
    Each contig has two ends (start / end). A join links two ends.
    We walk chains of joins to build scaffolds.

    Returns list of SeqRecord (scaffolds + unjoined contigs).
    """

    # Filter to high-confidence joins only
    confident = {
        k: v for k, v in joins.items() if v >= min_support
    }

    if not confident:
        print(f"[fungalflye] No joins with >= {min_support} supporting reads found")
        return list(records_dict.values())

    # Build adjacency: end -> (other_contig, other_side, support)
    adjacency = defaultdict(list)  # (contig, side) -> [(contig, side, support)]

    for (a, a_side), (b, b_side) in confident:
        support = confident[((a, a_side), (b, b_side))]
        adjacency[(a, a_side)].append((b, b_side, support))
        adjacency[(b, b_side)].append((a, a_side, support))

    used = set()
    scaffolds = []
    scaffold_num = 1

    def get_seq(contig, orientation):
        """Return sequence, reverse-complemented if orientation is '-'."""
        seq = records_dict[contig].seq
        return seq if orientation == "+" else seq.reverse_complement()

    # Greedy chain building: start from each unused contig
    for start_contig in list(records_dict.keys()):
        if start_contig in used:
            continue

        # Try to extend from the END of this contig
        chain = [(start_contig, "+")]   # (contig_id, orientation)
        used.add(start_contig)

        # Walk forward (from end of last contig)
        while True:
            last_contig, last_orient = chain[-1]

            # The "exit" side depends on orientation
            exit_side = "end" if last_orient == "+" else "start"
            key = (last_contig, exit_side)

            candidates = [
                c for c in adjacency.get(key, [])
                if c[0] not in used
            ]

            if not candidates:
                break

            # Pick highest-support join
            candidates.sort(key=lambda x: x[2], reverse=True)
            next_contig, next_side, support = candidates[0]

            # Determine orientation of next contig
            # If we're joining end->start, orientation is "+"
            # If we're joining end->end, we need to flip
            next_orient = "+" if next_side == "start" else "-"

            chain.append((next_contig, next_orient))
            used.add(next_contig)

        if len(chain) == 1:
            # No joins — keep as individual contig
            scaffolds.append(records_dict[start_contig])
        else:
            # Build scaffold sequence
            seq_parts = []
            for i, (ctg, orient) in enumerate(chain):
                seq_parts.append(str(get_seq(ctg, orient)))
                if i < len(chain) - 1:
                    seq_parts.append(GAP_SEQUENCE)

            scaffold_seq = "".join(seq_parts)
            scaffold_id = f"scaffold_{scaffold_num:04d}"
            scaffold_num += 1

            record = SeqRecord(
                Seq(scaffold_seq),
                id=scaffold_id,
                description=(
                    f"joined {len(chain)} contigs: "
                    + " + ".join(c for c, _ in chain)
                )
            )
            scaffolds.append(record)
            print(
                f"  Scaffold {scaffold_id}: joined "
                + " → ".join(
                    f"{c}({'fwd' if o == '+' else 'rev'})"
                    for c, o in chain
                )
            )

    return scaffolds


def run_scaffold(
    assembly,
    reads,
    outdir,
    threads,
    minimap2_preset,
    min_support=3,
    end_window=2000,
):
    """
    Main entry point for repeat-aware scaffolding.
    Returns path to scaffolded FASTA.
    """

    print("\n" + "=" * 60)
    print("🔬 Module 5 — Repeat-aware scaffolding")
    print(f"   Min supporting reads : {min_support}")
    print(f"   Terminal window      : {end_window} bp from each contig end")
    print("=" * 60)

    scaffold_dir = Path(outdir) / "scaffolding"
    scaffold_dir.mkdir(exist_ok=True)

    scaffolded_fasta = scaffold_dir / "scaffolded.fasta"

    if scaffolded_fasta.exists():
        print("[fungalflye] Existing scaffolded assembly detected — skipping")
        return scaffolded_fasta

    assembly = Path(assembly)

    # Load assembly
    records_dict = {
        r.id: r for r in SeqIO.parse(str(assembly), "fasta")
    }
    contig_lengths = {r.id: len(r.seq) for r in records_dict.values()}

    n_input = len(records_dict)
    print(f"\n[fungalflye] Input: {n_input} contigs")

    # Map reads
    paf = scaffold_dir / "reads.paf"
    if not paf.exists():
        print("[fungalflye] Mapping reads to assembly...")
        run(
            f"minimap2 -x {minimap2_preset} -t {threads} "
            f"--secondary=no {assembly} {reads} > {paf}"
        )

    # Parse bridges
    print("[fungalflye] Identifying contig bridges...")
    join_support = _parse_paf(paf, contig_lengths, end_window=end_window)

    # Report candidates
    confident_joins = {k: v for k, v in join_support.items() if v >= min_support}
    print(f"[fungalflye] Found {len(confident_joins)} high-confidence joins "
          f"(>= {min_support} supporting reads)")

    if confident_joins:
        print("\n  Top joins by read support:")
        top = sorted(confident_joins.items(), key=lambda x: x[1], reverse=True)[:10]
        for (a, a_side), (b, b_side) in [k for k, _ in top]:
            support = confident_joins[((a, a_side), (b, b_side))]
            print(f"    {a} ({a_side}) ↔ {b} ({b_side})  [{support} reads]")

    # Build scaffolds
    print("\n[fungalflye] Building scaffolds...")
    scaffolds = _build_scaffolds(records_dict, join_support, min_support)

    # Write output
    SeqIO.write(scaffolds, str(scaffolded_fasta), "fasta")

    n_out = len(scaffolds)
    joined = n_input - n_out
    print(f"\n✅ Scaffolding complete")
    print(f"   Input contigs   : {n_input}")
    print(f"   Output scaffolds: {n_out}  ({max(0, joined)} contigs joined)")
    print(f"   Output file     : {scaffolded_fasta}")

    return scaffolded_fasta


# ================================================
# Telomere-guided scaffolding
# ================================================

def _scan_telo_signal(seq, motif, window=3000, min_repeats=3):
    """Return (start_has_telo, end_has_telo) for a sequence."""
    rc = revcomp(motif)
    m  = len(motif)

    def count_tandem(s):
        best = 0
        for i in range(len(s) - m + 1):
            run, j = 0, i
            while j + m <= len(s):
                unit = s[j:j+m]
                if sum(x != y for x, y in zip(unit, motif)) <= 1 or \
                   sum(x != y for x, y in zip(unit, rc))    <= 1:
                    run += 1; j += m
                else:
                    break
            best = max(best, run)
        return best

    seq = seq.upper()
    start_rep = count_tandem(seq[:window])
    end_rep   = count_tandem(seq[-window:])
    return start_rep >= min_repeats, end_rep >= min_repeats


def run_telomere_scaffolding(
    assembly,
    reads,
    outdir,
    threads,
    minimap2_preset,
    telomere_motif="TTAGGG",
    min_support=2,
    end_window=3000,
    telo_window=3000,
    min_telo_repeats=3,
):
    """
    Telomere-guided scaffolding.

    Strategy:
      1. Scan all contigs for TTAGGG signal at their ends
      2. Identify small telomeric fragments (have telo signal, likely chromosome caps)
      3. Identify large contigs with uncapped ends (no telo signal on one or both ends)
      4. Use read bridging to attach telomeric fragments to uncapped chromosome ends
      5. Return path to telomere-extended assembly

    This specifically targets the biology where subtelomeric fragments
    get separated from chromosome bodies during assembly.
    """

    print("\n" + "=" * 60)
    print("🔬 Telomere-guided scaffolding")
    print(f"   Motif: {telomere_motif}  |  Min support: {min_support} reads")
    print("=" * 60)

    telo_dir = Path(outdir)
    telo_dir.mkdir(parents=True, exist_ok=True)

    out_fasta = telo_dir / "telo_scaffolded.fasta"

    if out_fasta.exists():
        print("[fungalflye] Existing telomere scaffold detected — skipping")
        return out_fasta

    assembly = Path(assembly)

    # Load all contigs
    records = {r.id: r for r in SeqIO.parse(str(assembly), "fasta")}
    lengths = {r.id: len(r.seq) for r in records.values()}

    # Step 1 — scan telomere signal on all contig ends
    print("\n[fungalflye] Scanning contig ends for telomere signal...")
    telo_status = {}
    for cid, rec in records.items():
        seq = str(rec.seq).upper()
        has_start, has_end = _scan_telo_signal(
            seq, telomere_motif, telo_window, min_telo_repeats
        )
        telo_status[cid] = {"start": has_start, "end": has_end}

    # Step 2 — classify contigs
    median_len = sorted(lengths.values())[len(lengths) // 2]
    size_threshold = max(median_len * 0.3, 50_000)

    telo_fragments = []   # small contigs with telo signal
    large_contigs  = []   # large contigs (chromosome bodies)

    for cid, ts in telo_status.items():
        has_any_telo = ts["start"] or ts["end"]
        is_small     = lengths[cid] < size_threshold

        if has_any_telo and is_small:
            telo_fragments.append(cid)
        if lengths[cid] >= size_threshold:
            large_contigs.append(cid)

    print(f"[fungalflye] Large chromosome contigs : {len(large_contigs)}")
    print(f"[fungalflye] Small telomeric fragments : {len(telo_fragments)}")

    if not telo_fragments:
        print("[fungalflye] No telomeric fragments found — skipping telo scaffolding")
        import shutil as _sh
        _sh.copy(str(assembly), str(out_fasta))
        return out_fasta

    # Step 3 — map reads to find bridges
    paf = telo_dir / "reads.paf"
    if not paf.exists():
        print("\n[fungalflye] Mapping reads for telomere bridge detection...")
        run(f"minimap2 -x {minimap2_preset} -t {threads} "
            f"--secondary=no {assembly} {reads} > {paf}")

    # Parse PAF for bridges between telo fragments and large contigs
    print("[fungalflye] Identifying telomere bridges...")

    read_hits = {}
    with open(paf) as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 12: continue
            read_id  = cols[0]
            t_name   = cols[5]
            t_start  = int(cols[7])
            t_end    = int(cols[8])
            t_len    = int(cols[6])
            mapq     = int(cols[11])
            if mapq < 10: continue

            at_start = t_start <= end_window
            at_end   = t_end   >= t_len - end_window

            if at_start and not at_end:
                read_hits.setdefault(read_id, []).append((t_name, "start"))
            elif at_end and not at_start:
                read_hits.setdefault(read_id, []).append((t_name, "end"))

    # Count bridges: telo_fragment_end -> large_contig_end
    bridge_support = {}
    for read_id, hits in read_hits.items():
        if len(hits) < 2: continue
        unique_hits = list(set(hits))
        ctg_names   = [h[0] for h in unique_hits]

        has_telo_frag  = any(c in telo_fragments for c in ctg_names)
        has_large_ctg  = any(c in large_contigs  for c in ctg_names)

        if not (has_telo_frag and has_large_ctg): continue

        for i in range(len(unique_hits)):
            for j in range(i + 1, len(unique_hits)):
                a, a_side = unique_hits[i]
                b, b_side = unique_hits[j]
                if a == b: continue
                if not ((a in telo_fragments and b in large_contigs) or
                        (b in telo_fragments and a in large_contigs)):
                    continue
                key = tuple(sorted([(a, a_side), (b, b_side)]))
                bridge_support[key] = bridge_support.get(key, 0) + 1

    confident = {k: v for k, v in bridge_support.items() if v >= min_support}

    if not confident:
        print(f"[fungalflye] No telomere bridges with >= {min_support} reads — skipping")
        import shutil as _sh
        _sh.copy(str(assembly), str(out_fasta))
        return out_fasta

    print(f"\n[fungalflye] Found {len(confident)} telomere bridges:")
    for (a, a_side), (b, b_side) in sorted(confident, key=lambda k: confident[k], reverse=True):
        support = confident[((a, a_side), (b, b_side))]
        frag    = a if a in telo_fragments else b
        large   = b if a in telo_fragments else a
        f_side  = a_side if a in telo_fragments else b_side
        l_side  = b_side if a in telo_fragments else a_side
        print(f"  {frag}({f_side}) → {large}({l_side})  [{support} reads]")

    # Step 4 — attach telomeric fragments to chromosome ends
    GAP = "N" * 100
    used_frags   = set()
    used_sides   = {}   # large_contig -> set of sides already extended
    new_records  = {}

    # Copy all large contigs as starting point
    for cid in large_contigs:
        new_records[cid] = str(records[cid].seq)
    used_sides = {cid: set() for cid in large_contigs}

    joins_made = 0

    for (a, a_side), (b, b_side) in sorted(confident.items(),
                                             key=lambda x: x[1], reverse=True):
        frag_id  = a if a in telo_fragments else b
        large_id = b if a in telo_fragments else a
        frag_side  = a_side if a in telo_fragments else b_side
        large_side = b_side if a in telo_fragments else a_side

        if frag_id in used_frags: continue
        if large_id not in new_records: continue
        if large_side in used_sides.get(large_id, set()): continue

        frag_seq  = str(records[frag_id].seq).upper()
        large_seq = new_records[large_id]

        # Orient the fragment: telo end should face outward
        frag_has_start, frag_has_end = _scan_telo_signal(
            frag_seq, telomere_motif, telo_window, min_telo_repeats
        )

        if large_side == "end":
            # Attach to end of large contig
            # Fragment should have telo on its END (facing out)
            if frag_has_end:
                new_records[large_id] = large_seq + GAP + frag_seq
            else:
                # Reverse complement so telo end faces out
                rc_seq = str(records[frag_id].seq.reverse_complement())
                new_records[large_id] = large_seq + GAP + rc_seq
        else:
            # Attach to start of large contig
            # Fragment should have telo on its START (facing out)
            if frag_has_start:
                new_records[large_id] = frag_seq + GAP + large_seq
            else:
                rc_seq = str(records[frag_id].seq.reverse_complement())
                new_records[large_id] = rc_seq + GAP + large_seq

        used_frags.add(frag_id)
        used_sides[large_id].add(large_side)
        joins_made += 1
        print(f"  Attached {frag_id} to {large_side} of {large_id}")

    # Step 5 — write output
    # Include large contigs (extended), unattached telo fragments, and other small contigs
    out_records = []
    for cid, seq in new_records.items():
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        out_records.append(SeqRecord(Seq(seq), id=cid, description=""))

    for cid, rec in records.items():
        if cid in large_contigs: continue       # already in new_records
        if cid in used_frags:    continue       # attached to a chromosome
        out_records.append(rec)                 # keep remaining small contigs

    SeqIO.write(out_records, str(out_fasta), "fasta")

    n_in  = len(records)
    n_out = len(out_records)
    print(f"\n✅ Telomere scaffolding complete")
    print(f"   Telomeric fragments attached : {joins_made}")
    print(f"   Input contigs  : {n_in}")
    print(f"   Output contigs : {n_out}  ({n_in - n_out} fewer)")
    print(f"   Output file    : {out_fasta}")

    return out_fasta
