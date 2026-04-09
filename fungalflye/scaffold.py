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
    print(f"\n[funcat] Running: {cmd}\n")
    subprocess.run(cmd, shell=True, check=True)


def _parse_paf(paf_path, contig_lengths, end_window=500, min_mapq=30,
               min_read_span_fraction=0.8, flagged_contigs=None):
    """
    Parse PAF alignments and find reads that anchor to the
    end-region of one contig and the start-region of another.

    Parameters
    ----------
    end_window : int
        How many bp from each end counts as a terminal region.
        Smaller = more specific, fewer false joins. Default 500bp.
    min_mapq : int
        Minimum mapping quality. 30+ rejects multi-mapping repeat reads.
    min_read_span_fraction : float
        Read alignment must span at least this fraction of end_window
        to count. Prevents short clips from being counted as bridges.
    flagged_contigs : set or None
        Contigs flagged as collapsed repeats by confidence scoring.
        These are excluded from bridging to prevent repeat-driven joins.

    Returns: dict of (contig_a, side_a, contig_b, side_b) -> support_count
    """

    if flagged_contigs is None:
        flagged_contigs = set()

    # Minimum alignment length needed to span end_window meaningfully
    min_aln_len = int(end_window * min_read_span_fraction)

    # How close to the absolute tip a read must reach to count as a bridge.
    # A read that anchors within end_window but never reaches the tip is
    # likely an internal read — not a true gap-bridging read.
    tip_distance = max(200, end_window // 10)

    # For each read, collect all contig terminal hits
    read_hits = defaultdict(list)   # read_id -> list of (contig, side)

    with open(paf_path) as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 12:
                continue

            read_id  = cols[0]
            t_name   = cols[5]
            t_start  = int(cols[7])
            t_end    = int(cols[8])
            t_len    = int(cols[6])
            mapq     = int(cols[11])
            aln_len  = t_end - t_start

            # Reject low-quality and short alignments
            if mapq < min_mapq:
                continue
            if aln_len < min_aln_len:
                continue

            # Skip collapsed repeat contigs entirely
            if t_name in flagged_contigs:
                continue

            t_len_actual = contig_lengths.get(t_name, t_len)

            # Contig must be large enough for end_window to make sense
            if t_len_actual < end_window * 3:
                continue

            # The read must ANCHOR in the terminal window (within end_window)
            # AND EXIT the contig — reach within tip_distance of the very end.
            # This rejects long internal reads that happen to sit near an end
            # but are not actually bridging out of the contig.
            in_start_window = t_start <= end_window
            exits_start     = t_start <= tip_distance

            in_end_window   = t_end >= t_len_actual - end_window
            exits_end       = t_end >= t_len_actual - tip_distance

            if in_start_window and exits_start and not in_end_window:
                read_hits[read_id].append((t_name, "start"))
            elif in_end_window and exits_end and not in_start_window:
                read_hits[read_id].append((t_name, "end"))
            # reads entirely internal or spanning whole contig — skip

    # Count contig-pair joins supported by reads
    # Enforce max ONE join partner per contig end (prevents fan-out from repeats)
    join_support = defaultdict(int)

    for read_id, hits in read_hits.items():
        if len(hits) < 2:
            continue

        # Deduplicate hits for this read
        unique_hits = list(set(hits))

        # Safety: ignore reads hitting more than 4 distinct contig ends
        # (those are almost certainly repeat reads spanning many loci)
        if len(unique_hits) > 4:
            continue

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

    # Enforce max ONE join per contig end — keep only the highest-support
    # partner for each (contig, side). This prevents a single contig end
    # from being joined to multiple partners (which signals a repeat, not a gap).
    best_per_end = {}   # (contig, side) -> best_key
    for key, support in sorted(confident.items(), key=lambda x: x[1], reverse=True):
        (a, a_side), (b, b_side) = key
        end_a = (a, a_side)
        end_b = (b, b_side)
        if end_a not in best_per_end and end_b not in best_per_end:
            best_per_end[end_a] = key
            best_per_end[end_b] = key

    confident = {v: confident[v] for v in best_per_end.values()}

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
    min_support=None,
    end_window=2000,
    confidence_tsv=None,
):
    """
    Main entry point for repeat-aware scaffolding.

    min_support: If None, auto-scales with estimated coverage (max(5, cov//10)).
                 Pass an explicit int to override.
    end_window:  Terminal window size in bp. Default 500 (tighter than before).
    confidence_tsv: Path to contig_confidence.tsv — flagged contigs are excluded.

    Returns path to scaffolded FASTA.
    """

    # Load flagged (collapsed repeat) contigs from confidence scoring
    flagged_contigs = set()
    if confidence_tsv and Path(confidence_tsv).exists():
        with open(confidence_tsv) as f:
            headers = None
            for line in f:
                parts = line.strip().split("\t")
                if headers is None:
                    headers = parts
                    continue
                row = dict(zip(headers, parts))
                if row.get("label") == "FLAG":
                    flagged_contigs.add(row["contig"])
        if flagged_contigs:
            print(f"[fungalflye] Excluding {len(flagged_contigs)} FLAG contigs "
                  f"from scaffolding: {', '.join(sorted(flagged_contigs))}")

    # Estimate coverage from reads to auto-scale min_support
    if min_support is None:
        try:
            import subprocess as _sp
            result = _sp.run(
                f"seqkit fx2tab -n -l {reads} 2>/dev/null | awk '{{sum+=$2}} END{{print sum}}'",
                shell=True, capture_output=True, text=True
            )
            total_bases = int(result.stdout.strip() or 0)
            assembly_size = sum(len(r.seq) for r in SeqIO.parse(str(assembly), "fasta"))
            if assembly_size > 0 and total_bases > 0:
                est_cov = total_bases / assembly_size
                min_support = max(5, int(est_cov // 10))
                print(f"[fungalflye] Estimated coverage: {est_cov:.1f}x → "
                      f"auto min_support = {min_support}")
            else:
                min_support = 5
        except Exception:
            min_support = 5
            print(f"[fungalflye] Could not estimate coverage — using min_support={min_support}")

    print("\n" + "=" * 60)
    print("🔬 Module 5 — Repeat-aware scaffolding")
    print(f"   Min supporting reads : {min_support}  (coverage-scaled)")
    print(f"   Terminal window      : {end_window} bp from each contig end")
    print(f"   Min MAPQ             : 30  (rejects multi-mapping reads)")
    print(f"   Flagged contigs skip : {len(flagged_contigs)}")
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
    join_support = _parse_paf(
        paf, contig_lengths,
        end_window=end_window,
        min_mapq=30,
        min_read_span_fraction=0.8,
        flagged_contigs=flagged_contigs,
    )

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

    # ── Small contig rescue pass ──────────────────────────────────────────
    # After the main scaffolding, check if any small GOOD contigs (<2Mb)
    # remain unjoined. These are likely subtelomeric fragments that had too
    # few bridging reads to pass the coverage-scaled min_support threshold.
    # Re-examine the PAF with a much lower threshold (min 3 reads) restricted
    # ONLY to joins between a small unjoined contig and a large contig (>=2Mb).
    # This prevents false joins between two small fragments while still
    # rescuing genuine chromosome-tip attachments.

    RESCUE_MAX_SMALL   = 1_000_000   # contigs smaller than 1Mb are rescue candidates
    RESCUE_MIN_LARGE   = 2_000_000   # must attach to a contig at least 2Mb (real chromosome body)
    RESCUE_MIN_SUPPORT = 3           # lower threshold for rescue pass only

    # Find which contigs are still unjoined after main scaffolding
    joined_ids = set()
    for s in scaffolds:
        if s.description and "joined" in s.description:
            for part in s.description.split(": ")[-1].split(" + "):
                joined_ids.add(part.strip())

    small_unjoined = {
        cid for cid, length in contig_lengths.items()
        if length < RESCUE_MAX_SMALL
        and cid not in flagged_contigs
        and cid not in joined_ids
    }
    large_contigs = {
        cid for cid, length in contig_lengths.items()
        if length >= RESCUE_MIN_LARGE
        and cid not in flagged_contigs
    }

    if small_unjoined and large_contigs:
        print(f"\n[fungalflye] Small contig rescue pass:")
        print(f"   Candidates : {len(small_unjoined)} small unjoined contigs")
        print(f"   Targets    : {len(large_contigs)} large chromosome-body contigs")

        # Re-parse PAF with rescue thresholds
        rescue_support = _parse_paf(
            paf, contig_lengths,
            end_window=end_window,
            min_mapq=30,
            min_read_span_fraction=0.5,   # more lenient — subtelomeric reads are harder
            flagged_contigs=flagged_contigs,
        )

        # Filter to only small↔large joins with >= RESCUE_MIN_SUPPORT
        rescue_joins = {}
        for key, support in rescue_support.items():
            (a, a_side), (b, b_side) = key
            if support < RESCUE_MIN_SUPPORT:
                continue
            small_to_large = (a in small_unjoined and b in large_contigs)
            large_to_small = (b in small_unjoined and a in large_contigs)
            if small_to_large or large_to_small:
                rescue_joins[key] = support

        if rescue_joins:
            print(f"   Found {len(rescue_joins)} rescue join(s):")
            for key, support in sorted(rescue_joins.items(), key=lambda x: x[1], reverse=True):
                (a, a_side), (b, b_side) = key
                print(f"     {a}({a_side}) ↔ {b}({b_side})  [{support} reads]")

            # Rebuild scaffolds including rescue joins
            combined_support = dict(join_support)
            combined_support.update(rescue_joins)
            scaffolds = _build_scaffolds(records_dict, combined_support, RESCUE_MIN_SUPPORT)
            print(f"   After rescue: {len(scaffolds)} contigs/scaffolds")
        else:
            print(f"   No rescue joins found — small contigs remain as fragments")
    # ── end rescue pass ───────────────────────────────────────────────────

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

def revcomp(s):
    table = str.maketrans("ACGT", "TGCA")
    return s.upper().translate(table)[::-1]


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
    min_support=5,
    end_window=500,
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
    # Classify contigs by telomere status, not size
    # telo_fragments = contigs with at least one telomeric end (chromosome caps)
    # large_contigs  = all contigs (we try to join any uncapped end)
    telo_fragments = []
    # large_contigs = contigs WITHOUT telomeric ends (attachment targets)
    # telo_fragments = contigs WITH at least one telomeric end (to be attached)
    large_contigs  = []

    for cid, ts in telo_status.items():
        if ts["start"] or ts["end"]:
            telo_fragments.append(cid)
        else:
            large_contigs.append(cid)

    print(f"[fungalflye] Contigs with telomeric ends : {len(telo_fragments)}")
    print(f"[fungalflye] Total contigs               : {len(large_contigs)}")

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

    min_aln_len = int(end_window * 0.8)   # read must span 80% of end_window
    tip_distance = max(200, end_window // 10)  # read must exit the contig tip

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
            aln_len  = t_end - t_start

            if mapq < 30: continue
            if aln_len < min_aln_len: continue

            in_start_window = t_start <= end_window
            exits_start     = t_start <= tip_distance
            in_end_window   = t_end >= t_len - end_window
            exits_end       = t_end >= t_len - tip_distance

            if in_start_window and exits_start and not in_end_window:
                read_hits.setdefault(read_id, []).append((t_name, "start"))
            elif in_end_window and exits_end and not in_start_window:
                read_hits.setdefault(read_id, []).append((t_name, "end"))

    # Count bridges: telo_fragment_end -> large_contig_end
    bridge_support = {}
    for read_id, hits in read_hits.items():
        if len(hits) < 2: continue
        unique_hits = list(set(hits))

        # Reject reads hitting too many distinct contig ends — repeat reads
        if len(unique_hits) > 4: continue

        ctg_names   = [h[0] for h in unique_hits]

        has_telo_frag  = any(c in telo_fragments for c in ctg_names)
        has_other      = len(set(ctg_names)) >= 2

        if not (has_telo_frag and has_other): continue

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
    for key in sorted(confident.keys(), key=lambda k: confident[k], reverse=True):
        (a, a_side), (b, b_side) = key
        support = confident[key]
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

    for key in sorted(confident.keys(), key=lambda k: confident[k], reverse=True):
        (a, a_side), (b, b_side) = key
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
    out_records = []

    # Write extended large contigs
    for cid, seq in new_records.items():
        out_records.append(SeqRecord(Seq(seq), id=cid, description=""))

    # Write any contigs that were NOT large contigs and NOT used as fragments
    for cid, rec in records.items():
        if cid in large_contigs:
            continue       # already written above
        if cid in used_frags:
            continue       # successfully attached — don't write separately
        out_records.append(rec)   # unattached small contigs kept as-is

    SeqIO.write(out_records, str(out_fasta), "fasta")

    n_in  = len(records)
    n_out = len(out_records)
    print(f"\n✅ Telomere scaffolding complete")
    print(f"   Telomeric fragments attached : {joins_made}")
    print(f"   Input contigs  : {n_in}")
    print(f"   Output contigs : {n_out}  ({n_in - n_out} fewer)")
    print(f"   Output file    : {out_fasta}")

    return out_fasta
