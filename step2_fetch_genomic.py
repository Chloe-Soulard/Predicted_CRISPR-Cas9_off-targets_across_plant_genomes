"""
Step 2 — Reconstruct the spliced CDS of each species' ortholog.

For every species, fetches the GenPept record of its top step 1 hit, parses the
`coded_by` qualifier for the source accession, coding intervals and strand, then
fetches each interval and concatenates them into the spliced CDS. Minus-strand
transcripts are concatenated in reverse genomic order so the CDS reads 5'->3'.

Curated species (config/manual_cds_overrides.json, written by
build_manual_overrides.py) are taken from that file and NCBI is not consulted.

Output: results/step2/{gene_slug}.json
    {"Arabidopsis thaliana": {"protein_accession": "NP_001031504",
                              "cds": "ATG...", "cds_length": 1134,
                              "genomic_accession": "NM_001036427.3",
                              "exons": [[59, 1192]], "strand": "+", "n_exons": 1},
     "Rubus occidentalis":   {... "source": "manual"},
     "SomeSpecies":          {"error": "no_hits", "reason": "..."},
     ...}

Note that `genomic_accession` is whatever `coded_by` points at, which for RefSeq
proteins is usually an mRNA rather than a chromosome. Real coding-exon boundaries
are added by step2b_exon_structure.py.

Run from pipeline/:
    python step2_fetch_genomic.py --gene ACT1
    python step2_fetch_genomic.py --all
    python step2_fetch_genomic.py --all --resume     # skip species already reconstructed
"""

from __future__ import annotations

import argparse

import ncbi
from paths import (
    MANUAL_OVERRIDES_FILE, STEP1_DIR, STEP2_DIR, gene_slug, load_json,
    resolve_genes, save_json,
)


def process_species(step1_entry: dict) -> dict:
    """Reconstruct one species' CDS, or return an error record."""
    hits = step1_entry.get("hits", [])
    if not hits:
        return {"error": "no_hits", "reason": "step 1 found no ortholog for this species"}

    accession = hits[0]["accession"]
    try:
        return ncbi.fetch_cds_from_protein(accession)
    except Exception as exc:
        return {"error": "fetch_failed", "reason": str(exc), "protein_accession": accession}


def run_gene(gene_name: str, resume: bool = False) -> None:
    slug = gene_slug(gene_name)
    step1_file = STEP1_DIR / f"{slug}.json"
    out_file   = STEP2_DIR / f"{slug}.json"

    if not step1_file.exists():
        print(f"[{gene_name}] no step1 file at {step1_file} — run step1_blast.py first")
        return

    step1     = load_json(step1_file)
    overrides = load_json(MANUAL_OVERRIDES_FILE).get(gene_name, {})
    results   = load_json(out_file)

    print(f"\n{'=' * 68}\nGene: {gene_name} | CDS for {len(step1)} species\n{'=' * 68}",
          flush=True)

    for species, step1_entry in step1.items():
        if species in overrides:
            results[species] = overrides[species]
            print(f"  [{species}] curated entry", flush=True)
            continue

        if resume and "cds" in results.get(species, {}):
            continue

        print(f"  fetching [{species}] ...", flush=True)
        result = process_species(step1_entry)
        results[species] = result

        if "error" in result:
            print(f"    ERROR: {result['reason']}", flush=True)
        else:
            print(f"    OK: {result['cds_length']} bp, {result['n_exons']} interval(s), "
                  f"strand={result['strand']}, {result['genomic_accession']}", flush=True)

        save_json(out_file, results)      # checkpoint after every species

    # Curated entries are copied in memory without a write, so flush once at the end
    # in case the last species processed was a curated one.
    save_json(out_file, results)

    n_ok = sum(1 for entry in results.values() if "cds" in entry)
    failures = {sp: entry["reason"] for sp, entry in results.items() if "error" in entry}
    print(f"\n[{gene_name}] {n_ok} CDS, {len(failures)} error(s)", flush=True)
    for species, reason in failures.items():
        print(f"  FAILED [{species}]: {reason}", flush=True)
    if failures:
        print("  Curate these by hand: make_manual_template.py -> build_manual_overrides.py",
              flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--gene", help="single gene (e.g. ACT1)")
    group.add_argument("--all", action="store_true", help="run every gene in the panel")
    parser.add_argument("--resume", action="store_true",
                        help="skip species that already have a CDS")
    args = parser.parse_args()

    ncbi.configure()
    STEP2_DIR.mkdir(parents=True, exist_ok=True)

    for gene_name in resolve_genes(args.gene):
        run_gene(gene_name, resume=args.resume)

    print("\nAll done.")


if __name__ == "__main__":
    main()
