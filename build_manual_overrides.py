"""
Step 1d — Convert the filled Excel form into config/manual_cds_overrides.json.

Only rows with a CDS are written; blank rows are reported and ignored. Step 2
reads this file and uses a curated entry in place of an NCBI lookup.

Every entry whose row also carries a protein sequence is validated by checking
that its CDS translates to exactly that protein — the gate that catches a
mis-pasted, truncated or frame-shifted sequence before it reaches the analysis.
A failure aborts the write unless --force is given.

Output schema (gene -> species -> entry):
    {"AtActin": {"Rubus occidentalis": {
        "protein_accession": "manual",          # source id, or 'manual' if none
        "protein": "MADGE...",                  # validated against the CDS
        "cds": "ATG...",
        "cds_length": 1338,
        "manual_genomic_sequence": "ATG...",    # only if a genomic sequence was given
        "manual_modif_reason": "...",           # only if a reason was given
        "source": "manual"}}}

The conservation analysis needs only the CDS. No NCBI accession or strand is
required — these species are precisely the ones NCBI does not cover.

Run from pipeline/:
    python build_manual_overrides.py
"""

from __future__ import annotations

import argparse
import re

import openpyxl
from Bio.Seq import Seq

from paths import MANUAL_OVERRIDES_FILE, ROOT, save_json

IN_FILE = ROOT / "manual_overrides_TEMPLATE.xlsx"
SHEET   = "Manual overrides"
FIRST_DATA_ROW = 3

# Column order must match make_manual_template.COLUMNS
GENE, SPECIES, GENOME, STATUS, ACCESSION, PROTEIN, CDS, GENOMIC_SEQ, REASON = range(9)

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWYBZJUOX*")


def clean_sequence(value) -> str:
    """Strip whitespace and upper-case a pasted sequence."""
    return re.sub(r"\s+", "", str(value)).upper() if value else ""


def is_protein_sequence(text: str) -> bool:
    """Distinguish a pasted protein from an accession or a free-text note."""
    return len(text) > 20 and set(text) <= AMINO_ACIDS


def translation_error(cds: str, protein: str) -> str | None:
    """None if the CDS translates to `protein`, else a description of the mismatch."""
    if len(cds) % 3:
        return f"CDS length {len(cds)} is not a multiple of 3"
    translated = str(Seq(cds).translate(to_stop=True))
    expected = protein.rstrip("*")
    if translated == expected:
        return None
    return (f"CDS translates to {len(translated)} aa, protein is {len(expected)} aa"
            if len(translated) != len(expected) else
            "CDS translation differs from the pasted protein")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("--force", action="store_true",
                        help="write the file even if a CDS fails to translate to its protein")
    args = parser.parse_args()

    if not IN_FILE.exists():
        raise SystemExit(f"{IN_FILE} not found — run make_manual_template.py first, "
                         f"then fill it in.")

    sheet = openpyxl.load_workbook(IN_FILE)[SHEET]

    overrides: dict[str, dict[str, dict]] = {}
    n_entries = n_validated = 0
    skipped: list[tuple[str, str]] = []
    failures: list[tuple[str, str, str]] = []

    for row in sheet.iter_rows(min_row=FIRST_DATA_ROW, values_only=True):
        if not row or not row[GENE]:
            continue
        gene    = str(row[GENE]).strip()
        species = str(row[SPECIES]).strip()
        cds     = clean_sequence(row[CDS])

        if not cds:
            skipped.append((gene, species))
            continue

        entry = {
            "protein_accession": str(row[ACCESSION] or "manual").strip(),
            "cds": cds,
            "cds_length": len(cds),
            "source": "manual",
        }

        protein = clean_sequence(row[PROTEIN])
        if is_protein_sequence(protein):
            entry["protein"] = protein
            error = translation_error(cds, protein)
            if error:
                failures.append((gene, species, error))
            else:
                n_validated += 1

        genomic = clean_sequence(row[GENOMIC_SEQ])
        if genomic:
            entry["manual_genomic_sequence"] = genomic
        reason = str(row[REASON] or "").strip()
        if reason:
            entry["manual_modif_reason"] = reason

        overrides.setdefault(gene, {})[species] = entry
        n_entries += 1

    print(f"{n_entries} curated entries across {len(overrides)} genes")
    print(f"  {n_validated} validated (CDS translates to the pasted protein)")
    if skipped:
        print(f"  {len(skipped)} rows skipped (no CDS):")
        for gene, species in skipped:
            print(f"    {gene} / {species}")
    if failures:
        print(f"  {len(failures)} FAILED validation:")
        for gene, species, error in failures:
            print(f"    {gene} / {species}: {error}")
        if not args.force:
            raise SystemExit("\nNothing written. Fix the rows above, or pass --force.")

    save_json(MANUAL_OVERRIDES_FILE, overrides)
    print(f"\nWrote {MANUAL_OVERRIDES_FILE}")


if __name__ == "__main__":
    main()
