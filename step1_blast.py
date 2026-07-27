"""
Step 1 — Find the best ortholog protein per species, for each bait gene.

One broad blastp per gene against refseq_protein restricted to Viridiplantae
(hitlist_size=500), from which the best hit (highest bit score) is kept for each
panel species. This costs 10 BLAST queries instead of 10 x 47 individual ones.

Species with no hit in the broad search are left empty here and picked up by
`step1_fallback.py`, which queries them one at a time.

Output: results/step1/{gene_slug}.json
    {"Zea mays": {"db": "refseq_protein_broad", "query": "Viridiplantae",
                  "hits": [{"accession": ..., "identity_pct": ..., "bits": ...}]},
     "Rubus occidentalis": {"db": "none", "hits": []},
     ...}

Run from pipeline/:
    python step1_blast.py --gene AtActin
    python step1_blast.py --all
    python step1_blast.py --all --resume     # skip genes already complete
"""

from __future__ import annotations

import argparse

from Bio.Blast import NCBIWWW, NCBIXML

import ncbi
from blast import extract_organism, hit_record, match_species
from paths import (
    STEP1_DIR, gene_slug, load_genes, load_json, load_species_genomes,
    resolve_genes, save_json,
)

HITLIST_SIZE = 500
EVALUE       = 1e-5


def best_hits_per_species(blast_record, panel: list[str]) -> dict[str, dict]:
    """Best hit (highest bit score) per panel species from one BLAST record."""
    best: dict[str, dict] = {}
    for alignment in blast_record.alignments:
        species = match_species(extract_organism(alignment.title), panel)
        if species is None:
            continue
        hit = hit_record(alignment)
        if species not in best or hit["bits"] > best[species]["bits"]:
            best[species] = hit
    return best


def blast_viridiplantae(gene_name: str, protein_seq: str,
                        panel: list[str]) -> dict[str, dict]:
    """One blastp across all green plants; return the best hit per panel species."""
    print(f"  blastp {gene_name} vs Viridiplantae refseq_protein "
          f"(hitlist={HITLIST_SIZE}) ...", flush=True)
    handle = NCBIWWW.qblast(
        "blastp", "refseq_protein", protein_seq,
        entrez_query="Viridiplantae[organism]",
        hitlist_size=HITLIST_SIZE, expect=EVALUE,
    )
    records = list(NCBIXML.parse(handle))
    if not records or not records[0].alignments:
        print("  WARNING: broad BLAST returned 0 hits", flush=True)
        return {}
    return best_hits_per_species(records[0], panel)


def run_gene(gene_name: str, protein_seq: str, panel: list[str],
             resume: bool = False) -> None:
    out_file = STEP1_DIR / f"{gene_slug(gene_name)}.json"
    results = load_json(out_file)

    if resume and {sp for sp, v in results.items() if v.get("hits")} == set(panel):
        print(f"[{gene_name}] all {len(panel)} species already have a hit, skipping")
        return

    print(f"\n{'=' * 68}\nGene: {gene_name}  ({len(protein_seq)} aa)\n{'=' * 68}",
          flush=True)

    best = blast_viridiplantae(gene_name, protein_seq, panel)
    print(f"  broad search: hits for {len(best)}/{len(panel)} species", flush=True)

    for species in panel:
        if species in best:
            results[species] = {"db": "refseq_protein_broad",
                                "query": "Viridiplantae",
                                "hits": [best[species]]}
        elif not results.get(species, {}).get("hits"):
            results[species] = {"db": "none", "hits": []}

    save_json(out_file, results)

    missing = [sp for sp in panel if not results[sp].get("hits")]
    print(f"\n[{gene_name}] {len(panel) - len(missing)}/{len(panel)} species with hits",
          flush=True)
    if missing:
        print(f"  no hit yet (run step1_fallback.py): {missing}", flush=True)
    print(f"  saved {out_file}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--gene", help="single gene (e.g. AtActin)")
    group.add_argument("--all", action="store_true", help="run every gene in the panel")
    parser.add_argument("--resume", action="store_true",
                        help="skip genes that already have a hit for every species")
    args = parser.parse_args()

    ncbi.configure()
    STEP1_DIR.mkdir(parents=True, exist_ok=True)

    genes = load_genes()
    panel = list(load_species_genomes().keys())

    for gene_name in resolve_genes(args.gene):
        run_gene(gene_name, genes[gene_name], panel, resume=args.resume)

    print("\nAll done.")


if __name__ == "__main__":
    main()
