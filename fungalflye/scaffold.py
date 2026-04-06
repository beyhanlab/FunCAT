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
