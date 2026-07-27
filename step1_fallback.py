"""
Step 1b — Targeted per-species BLAST for the species the broad search missed.

For each species still without a hit in results/step1/{gene}.json, runs blastp
restricted to that organism against refseq_protein and then the non-redundant
`nr` database, and writes the result back into the same file.

The `db` field records *why* a species has no hit, which matters for what happens
next:

    nr_targeted     resolved by this fallback
    none_targeted   confirmed absent — a clean query returned nothing.
                    These are the species to curate by hand (step1c).
    error_targeted  a query failed (network / server). NOT an absence: always
                    retried on the next run.

Results are written after every single species, so an interruption never costs
more than the one in-flight query.

Run from pipeline/:
    python step1_fallback.py --all
    python step1_fallback.py --all --resume                 # skip species already attempted
    python step1_fallback.py --all --resume --retry-failed  # also re-test confirmed absences
"""

from __future__ import annotations

import argparse
import time

from Bio.Blast import NCBIWWW, NCBIXML

from crispea import ncbi
from crispea.blast import hit_record
from crispea.paths import (
    STEP1_DIR, gene_slug, load_genes, load_json, resolve_genes, save_json,
)

HITLIST_SIZE  = 5
EVALUE        = 1e-5
DATABASES     = ("refseq_protein", "nr")
RETRY_SLEEP_S = 5
BETWEEN_DB_S  = 2

RESOLVED  = "nr_targeted"
ABSENT    = "none_targeted"
ERRORED   = "error_targeted"


def blast_one_species(protein_seq: str, species: str) -> dict:
    """blastp one bait protein against one organism; refseq_protein, then nr."""
    print(f"  [{species}] querying {' then '.join(DATABASES)} ...", flush=True)
    had_error = False

    for database in DATABASES:
        try:
            handle = NCBIWWW.qblast(
                "blastp", database, protein_seq,
                entrez_query=f"{species}[organism]",
                hitlist_size=HITLIST_SIZE, expect=EVALUE,
            )
            records = list(NCBIXML.parse(handle))
        except Exception as exc:
            print(f"    [{species}] {database} error: {exc}", flush=True)
            had_error = True
            time.sleep(RETRY_SLEEP_S)
            continue

        if records and records[0].alignments:
            hit = hit_record(records[0].alignments[0])
            print(f"    [{species}] {database}: {hit['accession']} "
                  f"({hit['identity_pct']:.1f}%)", flush=True)
            return {"db": RESOLVED, "query": species, "hits": [hit]}
        time.sleep(BETWEEN_DB_S)

    if had_error:
        print(f"    [{species}] no hit, but a query errored — not a confirmed absence",
              flush=True)
        return {"db": ERRORED, "hits": []}

    print(f"    [{species}] no hit in {' or '.join(DATABASES)}", flush=True)
    return {"db": ABSENT, "hits": []}


def species_to_retry(results: dict, resume: bool, retry_failed: bool) -> list[str]:
    """Which hit-less species this run should query."""
    missing = [sp for sp, entry in results.items() if not entry.get("hits")]
    if not resume:
        return missing
    # --resume alone skips anything a clean previous attempt settled; adding
    # --retry-failed re-tests confirmed absences too. Errored species are always
    # retried, since a failed query is not evidence of absence.
    settled = (RESOLVED,) if retry_failed else (RESOLVED, ABSENT)
    return [sp for sp in missing if results[sp].get("db") not in settled]


def run_gene(gene_name: str, protein_seq: str, resume: bool, retry_failed: bool) -> None:
    out_file = STEP1_DIR / f"{gene_slug(gene_name)}.json"
    if not out_file.exists():
        print(f"[{gene_name}] no step1 file — run step1_blast.py first")
        return

    results = load_json(out_file)
    todo = species_to_retry(results, resume, retry_failed)
    if not todo:
        print(f"[{gene_name}] nothing to retry")
        return

    print(f"\n{'=' * 68}\nGene: {gene_name} | fallback BLAST for {len(todo)} species"
          f"\n{'=' * 68}\n  {todo}", flush=True)

    for species in todo:
        results[species] = blast_one_species(protein_seq, species)
        save_json(out_file, results)      # checkpoint after every query

    n_resolved = sum(1 for sp in todo if results[sp].get("hits"))
    still_missing = [sp for sp in todo if not results[sp].get("hits")]
    print(f"\n[{gene_name}] {n_resolved} resolved, {len(still_missing)} still without a hit",
          flush=True)
    if still_missing:
        print(f"  {still_missing}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--gene", help="single gene (e.g. AtActin)")
    group.add_argument("--all", action="store_true", help="run every gene in the panel")
    parser.add_argument("--resume", action="store_true",
                        help="skip species a previous run already settled")
    parser.add_argument("--retry-failed", action="store_true",
                        help="with --resume, also re-test species confirmed absent")
    args = parser.parse_args()

    ncbi.configure()
    genes = load_genes()

    for gene_name in resolve_genes(args.gene):
        run_gene(gene_name, genes[gene_name], args.resume, args.retry_failed)

    print("\nAll done.")


if __name__ == "__main__":
    main()
