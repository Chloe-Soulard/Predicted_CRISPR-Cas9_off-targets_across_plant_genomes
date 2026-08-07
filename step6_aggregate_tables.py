"""
Step 6a — Aggregate the CRISPOR results into the analysis tables.

Reads results/step5/*.json and joins it to the species metadata, producing the
tidy tables the statistics and figures are built from. One window = one gene x
species conserved region; one guide = one SpCas9/NGG protospacer inside a window.

Outputs (all four from a single pass over the same records):
    results/offtarget_aggregate.csv             per window, raw counts
    results/genic_offtarget_breakdown.csv       per window, genic vs intergenic split
    results/R_export/offtargets_perwindow.csv   per window, tidy, for R
    results/R_export/offtargets_perguide.csv    per guide, tidy, for R

Two conventions carry through every table, and both matter for interpretation:

  * Off-target load is reported **per guide**. The number of guides in a 100 bp
    window varies with the sequence, so raw window totals are not comparable
    across species; dividing by the number of guides makes them so.

  * Only guides **present in the target genome** are counted. CRISPOR gives a
    guide absent from the assembly a specificity of -1 and no off-targets;
    including them would dilute the means with structural zeros.

Windows CRISPOR could not map to their genome are excluded entirely and reported
on the console — they measure an assembly mismatch, not an off-target rate.

Run from pipeline/:
    python step6_aggregate_tables.py
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from paths import (
    DEFAULT_PAM, R_EXPORT_DIR, RESULTS_DIR, gene_names, gene_slug,
    load_json, load_species_metadata, pam_suffix, ploidy_group, step5_dir,
)

MISMATCH_LEVELS = range(5)          # CRISPOR reports off-targets up to 4 mismatches

# Decimal places: 2 for the human-readable tables, 4-5 for the tables R consumes.
DP_REPORT   = 2
DP_ANALYSIS = 4
DP_FRACTION = 5


def per_guide(total: int, n_guides: int, decimals: int) -> float:
    return round(total / n_guides, decimals)


def fraction(part: int, whole: int, decimals: int) -> float:
    """Genic fraction. A window with no off-targets at all reports 0."""
    return round(part / whole, decimals) if whole else 0.0


def tally(guides: list[dict]) -> dict[str, int]:
    """Sum off-target counts over the guides of one window."""
    totals = {
        "n_guides":    len(guides),
        "total":       sum(g["offtarget_count"] for g in guides),
        "genic":       sum(g["ot_genic"] for g in guides),
        "intergenic":  sum(g["ot_intergenic"] for g in guides),
        "unannotated": sum(g["ot_unannotated"] for g in guides),
    }
    for level in MISMATCH_LEVELS:
        totals[f"mm{level}"] = sum(g["ot_by_mismatch"].get(str(level), 0) for g in guides)
    return totals


def collect_windows(metadata: dict[str, dict],
                    directory) -> tuple[list[dict], list[str]]:
    """Load every analysable window, in panel gene order then species order.

    Returns (windows, excluded, repeats) where each window carries its species
    metadata, its usable guides and their tallies; `excluded` describes what was
    dropped and why, and `repeats` the guides CRISPOR could not enumerate.
    """
    windows, excluded, repeats = [], [], []

    for gene in gene_names():
        path = directory / f"{gene_slug(gene)}.json"
        if not path.exists():
            excluded.append(f"{gene}: no step5 results file")
            continue

        results = load_json(path)
        for species in sorted(results):
            record = results[species]

            if record.get("status") != "ok":
                excluded.append(f"{gene} / {species}: status={record.get('status')}")
                continue
            if not record.get("found_in_genome"):
                # A window with no guides did map; it simply carries no PAM site.
                # Reporting that as "not found" would read as an assembly mismatch.
                reason = ("no PAM site in the window"
                          if record.get("found_state") == "no_guides"
                          else f"window not found in {record.get('genome')}")
                excluded.append(f"{gene} / {species}: {reason}")
                continue
            if species not in metadata:
                excluded.append(f"{gene} / {species}: no genome size / ploidy on file")
                continue

            # Guides CRISPOR refused to enumerate ("repeated region, too unspecific")
            # carry a 0 that means "unknown", not "none". Averaging them in would
            # credit the most repetitive guides as the most specific, and they fall
            # disproportionately in the largest genomes — the very axis under test.
            in_genome = [g for g in record.get("guides", []) if g.get("in_genome")]
            guides = [g for g in in_genome if not g.get("repeat_unspecific")]
            n_repeat = len(in_genome) - len(guides)
            if n_repeat:
                repeats.append(f"{gene} / {species}: {n_repeat} repeat guide(s) dropped")
            if not guides:
                excluded.append(f"{gene} / {species}: no guides present in the genome"
                                + (f" ({n_repeat} repeat-only)" if n_repeat else ""))
                continue

            windows.append({
                "gene": gene,
                "species": species,
                "genome_size_mb": metadata[species]["genome_size_mb"],
                "ploidy": metadata[species]["ploidy"],
                "ploidy_group": ploidy_group(metadata[species]["ploidy"]),
                "record": record,
                "guides": guides,
                "totals": tally(guides),
            })

    return windows, excluded, repeats


# ── Table builders ────────────────────────────────────────────────────────────

def aggregate_rows(windows: list[dict]) -> list[dict]:
    """results/offtarget_aggregate.csv — per window, raw counts."""
    rows = []
    for w in windows:
        t, n = w["totals"], w["totals"]["n_guides"]
        rows.append({
            "gene": w["gene"], "species": w["species"],
            "genome_size_mb": w["genome_size_mb"], "ploidy": w["ploidy"],
            "found_in_genome": w["record"].get("found_in_genome"),
            "needs_manual_recovery": w["record"].get("needs_genomic_dna"),
            "n_guides": n, "ot_total": t["total"],
            "ot_per_guide": per_guide(t["total"], n, DP_REPORT),
            **{f"ot_mm{level}": t[f"mm{level}"] for level in MISMATCH_LEVELS},
            "ot_genic": t["genic"], "ot_intergenic": t["intergenic"],
        })
    return rows


def genic_rows(windows: list[dict]) -> list[dict]:
    """results/genic_offtarget_breakdown.csv — per window, locus-class split."""
    rows = []
    for w in windows:
        t, n = w["totals"], w["totals"]["n_guides"]
        rows.append({
            "gene": w["gene"], "species": w["species"],
            "genome_size_mb": w["genome_size_mb"], "ploidy": w["ploidy"],
            "n_guides": n,
            "total_ot": t["total"], "genic_ot": t["genic"],
            "intergenic_ot": t["intergenic"], "unannotated_ot": t["unannotated"],
            "fraction_genic": fraction(t["genic"], t["total"], DP_ANALYSIS),
            "total_ot_per_guide": per_guide(t["total"], n, DP_REPORT),
            "genic_ot_per_guide": per_guide(t["genic"], n, DP_REPORT),
            "intergenic_ot_per_guide": per_guide(t["intergenic"], n, DP_REPORT),
        })
    return rows


def per_window_rows(windows: list[dict]) -> list[dict]:
    """R_export/offtargets_perwindow.csv — the main analysis table."""
    rows = []
    for w in windows:
        t, n = w["totals"], w["totals"]["n_guides"]
        counts = [g["offtarget_count"] for g in w["guides"]]
        rows.append({
            "gene": w["gene"], "species": w["species"],
            "genome_size_mb": w["genome_size_mb"], "ploidy": w["ploidy"],
            "ploidy_group": w["ploidy_group"],
            "n_guides": n, "total_ot": t["total"],
            "total_ot_per_guide": per_guide(t["total"], n, DP_ANALYSIS),
            "genic_ot": t["genic"], "intergenic_ot": t["intergenic"],
            "unannotated_ot": t["unannotated"],
            "genic_ot_per_guide": per_guide(t["genic"], n, DP_ANALYSIS),
            "intergenic_ot_per_guide": per_guide(t["intergenic"], n, DP_ANALYSIS),
            "fraction_genic": fraction(t["genic"], t["total"], DP_FRACTION),
            # Best and worst guide in the window: the spread a designer would face.
            "min_ot": min(counts), "max_ot": max(counts),
            **{f"mm{level}_per_guide": per_guide(t[f"mm{level}"], n, DP_ANALYSIS)
               for level in MISMATCH_LEVELS},
        })
    return rows


def per_guide_rows(windows: list[dict]) -> list[dict]:
    """R_export/offtargets_perguide.csv — one row per gRNA."""
    rows = []
    for w in windows:
        for guide in w["guides"]:
            total = guide["offtarget_count"]
            rows.append({
                "gene": w["gene"], "species": w["species"],
                "guide_id": guide["guide_id"], "target_seq": guide["target_seq"],
                "genome_size_mb": w["genome_size_mb"], "ploidy": w["ploidy"],
                "ploidy_group": w["ploidy_group"],
                "mit_spec": guide["mit_spec"], "cfd_spec": guide["cfd_spec"],
                "total_ot": total, "genic_ot": guide["ot_genic"],
                "intergenic_ot": guide["ot_intergenic"],
                "unannotated_ot": guide["ot_unannotated"],
                "fraction_genic": fraction(guide["ot_genic"], total, DP_FRACTION),
                **{f"mm{level}": guide["ot_by_mismatch"].get(str(level), 0)
                   for level in MISMATCH_LEVELS},
            })
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"  {path.relative_to(RESULTS_DIR.parent)}  ({len(rows)} rows)")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("--pam", default=DEFAULT_PAM,
                        help=f"nuclease to aggregate (default {DEFAULT_PAM}; TTTV = "
                             f"Cas12a). Reads results/step5<suffix> and writes the "
                             f"tables with the same suffix, so nucleases never mix.")
    args = parser.parse_args()

    pam = args.pam.strip().upper()
    suffix = pam_suffix(pam)
    directory = step5_dir(pam)
    if not directory.exists():
        raise SystemExit(f"{directory} does not exist — run "
                         f"step5_crispor.py --pam {pam} first.")

    metadata = load_species_metadata()
    windows, excluded, repeats = collect_windows(metadata, directory)
    if not windows:
        raise SystemExit(f"No analysable windows in {directory} — "
                         f"run step5_crispor.py --pam {pam} first.")

    genes = sorted({w["gene"] for w in windows})
    species = sorted({w["species"] for w in windows})
    n_guides = sum(w["totals"]["n_guides"] for w in windows)
    print(f"{len(windows)} windows ({len(genes)} genes x {len(species)} species), "
          f"{n_guides} guides in genome")

    if repeats:
        n_repeat = sum(int(r.split(": ")[1].split()[0]) for r in repeats)
        print(f"\n{n_repeat} guide(s) across {len(repeats)} window(s) dropped as "
              f"repeats — CRISPOR reports 'repeated region, too unspecific' and "
              f"gives 0 off-targets, which means unknown, not none:")
        for reason in repeats:
            print(f"  {reason}")

    if excluded:
        print(f"\n{len(excluded)} window(s) excluded:")
        for reason in excluded:
            print(f"  {reason}")

    print("\nWriting tables:")
    write_csv(RESULTS_DIR / f"offtarget_aggregate{suffix}.csv", aggregate_rows(windows))
    write_csv(RESULTS_DIR / f"genic_offtarget_breakdown{suffix}.csv", genic_rows(windows))
    write_csv(R_EXPORT_DIR / f"offtargets_perwindow{suffix}.csv", per_window_rows(windows))
    write_csv(R_EXPORT_DIR / f"offtargets_perguide{suffix}.csv", per_guide_rows(windows))


if __name__ == "__main__":
    main()
