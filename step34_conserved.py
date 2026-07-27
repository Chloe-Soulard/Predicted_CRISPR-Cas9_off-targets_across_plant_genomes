"""
Steps 3-4 — Select one conserved, fully exonic 100 bp window per gene.

Step 3  Star alignment: every species CDS is aligned to the Arabidopsis reference
        CDS (BioPython PairwiseAligner, global, NUC.4.4, gap-open -10,
        gap-extend -0.5), giving for each reference position the aligned base and
        index in each species.

Step 4  Window choice. A 100 bp window is slid along the reference CDS and scored
        by the number of positions matching the reference, summed over species.
        Rather than taking the single most-conserved window and then shifting it
        per species, the window chosen is the highest-scoring one that lies
        entirely within a single exon in *every* species, so no species needs a
        shift across an exon junction. If no window is single-exon everywhere
        (AtUBQ2, whose polyubiquitin repeats put junctions throughout the CDS),
        the window maximising the number of single-exon species wins, ties broken
        by conservation.

        Each species then contributes exactly `window` bp of single-exon CDS
        homologous to that reference window. Species whose exon still cannot
        accommodate it are flagged rather than silently truncated.

Output:
    results/step34/{gene_slug}.json                  windows + provenance per species
    results/step34/short_exon_flags.json             pairs needing manual review
    results/alignments/{gene_slug}_alignment.fasta   star alignment, reference-framed
    results/alignments/{gene_slug}_conservation.txt  per-position identity over the window

This step is purely local: no network access.

Run from pipeline/:
    python step34_conserved.py --gene AtActin
    python step34_conserved.py --all
"""

from __future__ import annotations

import argparse
from collections import Counter

from Bio.Align import PairwiseAligner, substitution_matrices

from paths import (
    ALIGN_DIR, AT_CDS_FILE, STEP2_DIR, STEP34_DIR, WINDOW_BP, gene_slug,
    load_json, resolve_genes, save_json,
)

FASTA_LINE_WIDTH = 60


# ── Star alignment ────────────────────────────────────────────────────────────

def make_aligner() -> PairwiseAligner:
    """Global nucleotide aligner, NUC.4.4 scoring."""
    aligner = PairwiseAligner(mode="global")
    aligner.substitution_matrix = substitution_matrices.load("NUC.4.4")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    return aligner


def align_to_reference(reference: str, query: str
                       ) -> tuple[dict[int, str], dict[int, int]]:
    """Align one species CDS to the reference CDS.

    Returns (base_map, pos_map) keyed by 0-based reference position:
        base_map[ref_pos]  aligned query base, or '-' where the query has a gap
        pos_map[ref_pos]   index of that base in `query` (absent at gaps)
    """
    alignment = make_aligner().align(reference, query)[0]

    base_map: dict[int, str] = {}
    pos_map: dict[int, int] = {}
    for (ref_start, ref_end), (query_start, _) in zip(*alignment.aligned):
        ref_start, query_start = int(ref_start), int(query_start)
        for offset in range(int(ref_end) - ref_start):
            base_map[ref_start + offset] = query[query_start + offset]
            pos_map[ref_start + offset] = query_start + offset

    for pos in range(len(reference)):
        base_map.setdefault(pos, "-")          # uncovered reference position = gap
    return base_map, pos_map


def build_column_maps(reference: str, species_cds: dict[str, str]
                      ) -> tuple[dict[str, dict[int, str]], dict[str, dict[int, int]]]:
    base_maps, pos_maps = {}, {}
    for species, cds in species_cds.items():
        base_maps[species], pos_maps[species] = align_to_reference(reference, cds)
    return base_maps, pos_maps


# ── Window selection ──────────────────────────────────────────────────────────

def match_counts(reference: str, base_maps: dict[str, dict[int, str]]) -> list[int]:
    """Per reference position, how many species carry the reference base."""
    return [sum(1 for m in base_maps.values() if m.get(pos, "-") == base)
            for pos, base in enumerate(reference)]


def best_exonic_window(reference: str, base_maps: dict[str, dict[int, str]],
                       species_cds: dict[str, str], junctions: dict[str, list[int]],
                       pos_maps: dict[str, dict[int, int]], window: int = WINDOW_BP):
    """Choose the window: fully exonic first, conservation second.

    Returns (start, window, score, n_exonic, fully_exonic).
    """
    n_species = len(base_maps)
    counts = match_counts(reference, base_maps)

    # Prefix sums give the conservation score of every window in one pass.
    prefix = [0]
    for count in counts:
        prefix.append(prefix[-1] + count)

    candidates = []                                  # (start, score, n_exonic)
    for start in range(len(reference) - window + 1):
        n_exonic = sum(
            1 for species in species_cds
            if extract_exonic_window(species_cds[species], junctions[species],
                                     pos_maps[species], start, start + window,
                                     window).get("source") == "exonic"
        )
        candidates.append((start, prefix[start + window] - prefix[start], n_exonic))

    fully_exonic = [c for c in candidates if c[2] == n_species]
    if fully_exonic:
        start, score, n_exonic = max(fully_exonic, key=lambda c: c[1])
        return start, window, score, n_exonic, True

    start, score, n_exonic = max(candidates, key=lambda c: (c[2], c[1]))
    return start, window, score, n_exonic, False


# ── Per-species window extraction ─────────────────────────────────────────────

def exon_bounds_of(pos: int, junctions: list[int], cds_len: int) -> tuple[int, int, int]:
    """Half-open bounds [start, end) and index of the exon containing `pos`."""
    starts = [0] + list(junctions)
    ends = list(junctions) + [cds_len]
    for index, (start, end) in enumerate(zip(starts, ends)):
        if start <= pos < end:
            return start, end, index
    return starts[-1], cds_len, len(starts) - 1        # past the end: last exon


def extract_exonic_window(species_cds: str, junctions: list[int],
                          pos_map: dict[int, int], ref_start: int, ref_end: int,
                          window: int = WINDOW_BP) -> dict:
    """Take exactly `window` bp of single-exon CDS homologous to the reference window.

        exonic          the homologous position already sits inside one exon
        exonic_shifted  it straddled a junction, so the window slid within the
                        nearest exon that can hold it (`shift_bp` records how far)
        short_exon      no covering exon is long enough — flagged for review
        no_mapping      the reference window does not align to this species at all
    """
    cds_len = len(species_cds)

    # Map the reference window start onto this species, searching outward for the
    # nearest aligned reference position when the exact one falls in a gap.
    start = None
    for distance in range(ref_end - ref_start + 50):
        for candidate in (ref_start + distance, ref_start - distance):
            if candidate in pos_map:
                start = pos_map[candidate] - (candidate - ref_start)
                break
        if start is not None:
            break
    if start is None:
        return {"flag": "no_mapping",
                "reason": "reference window does not align to this species"}
    start = max(0, min(start, max(0, cds_len - 1)))

    exon_start, exon_end, exon_index = exon_bounds_of(start, junctions, cds_len)

    if start + window <= exon_end and start + window <= cds_len:
        return {"sequence": species_cds[start:start + window], "source": "exonic",
                "exon_index": exon_index, "species_cds_start": start,
                "species_cds_end": start + window, "shift_bp": 0, "length": window}

    # The window crosses a junction or runs off the CDS. Consider sliding left
    # within the current exon, or left-aligning in the exon the window ends in,
    # and take whichever moves least.
    options = []                                       # (start, exon_index, distance)
    if exon_end - exon_start >= window and exon_end <= cds_len:
        shifted = min(max(start, exon_start), exon_end - window)
        options.append((shifted, exon_index, abs(shifted - start)))

    next_start, next_end, next_index = exon_bounds_of(
        min(start + window - 1, cds_len - 1), junctions, cds_len)
    if ((next_start, next_end) != (exon_start, exon_end)
            and next_end - next_start >= window and next_end <= cds_len):
        options.append((next_start, next_index, abs(next_start - start)))

    if not options:
        return {"flag": "short_exon",
                "reason": f"no single exon >= {window} bp covers the conserved window",
                "species_cds_start": start}

    chosen, exon_index, _ = min(options, key=lambda o: o[2])
    return {"sequence": species_cds[chosen:chosen + window], "source": "exonic_shifted",
            "exon_index": exon_index, "species_cds_start": chosen,
            "species_cds_end": chosen + window, "shift_bp": chosen - start,
            "length": window}


# ── Alignment outputs ─────────────────────────────────────────────────────────

def write_fasta(lines, header: str, sequence: str) -> None:
    lines.append(header)
    for i in range(0, len(sequence), FASTA_LINE_WIDTH):
        lines.append(sequence[i:i + FASTA_LINE_WIDTH])


def write_alignment_fasta(gene_name: str, reference: str, species_cds: dict[str, str],
                          base_maps: dict[str, dict[int, str]],
                          start: int, end: int) -> None:
    """Write the star alignment in reference coordinates (one column per reference base)."""
    out_file = ALIGN_DIR / f"{gene_slug(gene_name)}_alignment.fasta"
    n_ref = len(reference)

    lines: list[str] = []
    write_fasta(lines, f">Arabidopsis_thaliana|AT_reference|{n_ref}bp", reference)
    lines.append(f"# Conserved window: AT CDS positions {start}-{end} ({end - start} bp)")
    for species, cds in species_cds.items():
        gapped = "".join(base_maps[species].get(pos, "-") for pos in range(n_ref))
        write_fasta(lines, f">{species.replace(' ', '_')}|{len(cds)}bp", gapped)

    out_file.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"  wrote {out_file}", flush=True)


def write_conservation_table(gene_name: str, reference: str,
                             base_maps: dict[str, dict[int, str]],
                             start: int, end: int, score: int) -> None:
    """Per-position identity across the selected window."""
    out_file = ALIGN_DIR / f"{gene_slug(gene_name)}_conservation.txt"
    n_species = len(base_maps)
    max_score = (end - start) * n_species

    lines = [
        f"Gene: {gene_name}",
        f"Conserved window: AT CDS positions {start}–{end - 1} ({end - start} bp)",
        f"Conservation score: {score} / {max_score} ({100 * score / max_score:.1f}%)",
        f"Species aligned: {n_species}",
        "",
        f"{'Pos':>5}  {'AT':>2}  {'Count':>5}  {'Frac':>5}  {'Consensus':>9}",
        "-" * 40,
    ]
    for pos in range(start, end):
        reference_base = reference[pos]
        bases = [base_maps[sp].get(pos, "-") for sp in base_maps]
        count = sum(1 for base in bases if base == reference_base)
        fraction = count / n_species
        consensus = Counter(b for b in bases if b != "-").most_common(1)
        marker = " ***" if fraction == 1.0 else ""
        lines.append(f"{pos:>5}  {reference_base:>2}  {count:>5}  {fraction:>5.2f}  "
                     f"{(consensus[0][0] if consensus else '?'):>9}{marker}")

    out_file.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"  wrote {out_file}", flush=True)


# ── Per-gene driver ───────────────────────────────────────────────────────────

def run_gene(gene_name: str, window: int = WINDOW_BP) -> list[str]:
    """Select the window for one gene. Returns the species flagged short-exon."""
    slug = gene_slug(gene_name)
    step2_file = STEP2_DIR / f"{slug}.json"
    if not step2_file.exists():
        print(f"[{gene_name}] no step2 file — run step2_fetch_genomic.py first")
        return []

    at_cds_entry = load_json(AT_CDS_FILE).get(gene_name, {})
    if "cds" not in at_cds_entry:
        print(f"[{gene_name}] no reference CDS in {AT_CDS_FILE.name} — "
              f"run fetch_at_cds.py first")
        return []
    reference = at_cds_entry["cds"]

    step2 = load_json(step2_file)
    species_cds = {sp: entry["cds"] for sp, entry in step2.items() if "cds" in entry}

    print(f"\n{'=' * 68}\nGene: {gene_name} | reference CDS {len(reference)} bp"
          f"\n{'=' * 68}", flush=True)
    print(f"  usable CDS: {len(species_cds)}/{len(step2)} species", flush=True)
    if len(species_cds) < 2:
        print("  ERROR: fewer than 2 species with a CDS — cannot align", flush=True)
        return []

    print("  building star alignment ...", flush=True)
    base_maps, pos_maps = build_column_maps(reference, species_cds)

    print(f"  scanning for the best {window} bp window (fully exonic preferred) ...",
          flush=True)
    junctions = {sp: step2[sp].get("cds_junctions", []) for sp in species_cds}
    start, window, score, n_exonic_at_window, fully_exonic = best_exonic_window(
        reference, base_maps, species_cds, junctions, pos_maps, window)
    end = start + window
    max_score = window * len(species_cds)

    print(f"  window: AT CDS [{start}, {end}) "
          f"score={score}/{max_score} ({100 * score / max_score:.1f}%) | "
          f"{n_exonic_at_window}/{len(species_cds)} species single-exon "
          f"({'fully exonic' if fully_exonic else 'best achievable'})", flush=True)

    write_alignment_fasta(gene_name, reference, species_cds, base_maps, start, end)
    write_conservation_table(gene_name, reference, base_maps, start, end, score)

    print(f"  extracting {window} bp per species ...", flush=True)
    cds_sequences: dict[str, str] = {}
    exonic_windows: dict[str, dict] = {}
    flagged: list[str] = []
    n_exonic = n_shifted = 0

    for species in sorted(species_cds):
        # The gap-stripped alignment slice, kept for reference alongside the
        # single-exon window that is actually submitted to CRISPOR.
        cds_sequences[species] = "".join(
            base_maps[species].get(pos, "-") for pos in range(start, end)
        ).replace("-", "")

        info = extract_exonic_window(species_cds[species], junctions[species],
                                     pos_maps[species], start, end, window)
        # Species with no real gene model are treated as single-exon; mark them so
        # the window can be checked before it is trusted.
        if step2[species].get("exon_source") in ("assumed_single_exon", "error"):
            info["exon_unverified"] = True
        exonic_windows[species] = info

        if "flag" in info:
            flagged.append(species)
            print(f"    [{species}] FLAG {info['flag']}: {info.get('reason', '')}",
                  flush=True)
        elif info["source"] == "exonic":
            n_exonic += 1
        else:
            n_shifted += 1
            print(f"    [{species}] shifted {info['shift_bp']:+d} bp to stay within "
                  f"one exon", flush=True)

    print(f"  windows: {n_exonic} exonic, {n_shifted} shifted, {len(flagged)} flagged",
          flush=True)

    out_file = STEP34_DIR / f"{slug}.json"
    save_json(out_file, {
        "gene": gene_name,
        "cds_window_start": start,
        "cds_window_end": end,
        "cds_window_length": window,
        "conservation_score": score,
        "conservation_max_possible": max_score,
        "conservation_pct": round(100 * score / max_score, 1),
        "n_species_aligned": len(species_cds),
        "window_fully_exonic": fully_exonic,
        "n_exonic_at_window": n_exonic_at_window,
        "cds_sequences": cds_sequences,
        "exonic_windows": exonic_windows,
    })
    print(f"  saved {out_file}", flush=True)
    return flagged


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--gene", help="single gene (e.g. AtActin)")
    group.add_argument("--all", action="store_true", help="run every gene in the panel")
    parser.add_argument("--window", type=int, default=WINDOW_BP,
                        help=f"window length in bp (default {WINDOW_BP})")
    args = parser.parse_args()

    STEP34_DIR.mkdir(parents=True, exist_ok=True)
    ALIGN_DIR.mkdir(parents=True, exist_ok=True)

    flags: dict[str, list[str]] = {}
    for gene_name in resolve_genes(args.gene):
        flagged = run_gene(gene_name, window=args.window)
        if flagged:
            flags[gene_name] = flagged

    flags_file = STEP34_DIR / "short_exon_flags.json"
    save_json(flags_file, flags)

    print("\nAll done.")
    if flags:
        total = sum(len(v) for v in flags.values())
        print(f"  {total} short-exon flag(s) needing manual review -> {flags_file}")
        for gene_name, species_list in flags.items():
            for species in species_list:
                print(f"    {gene_name} / {species}")
    else:
        print("  No short-exon flags.")


if __name__ == "__main__":
    main()
