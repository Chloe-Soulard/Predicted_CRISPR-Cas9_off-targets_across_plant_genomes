"""
Step 0b — Fetch the Arabidopsis thaliana reference CDS for each bait gene.

The reference CDS is the coordinate system for the whole study: species CDSs are
aligned to it and the conserved window is defined by its position in it.

For each gene:
    1. blastp the bait protein against refseq_protein, restricted to A. thaliana
    2. fetch the top hit's GenPept record and parse its `coded_by` qualifier
    3. fetch each coding exon and concatenate them into the spliced CDS

Output: config/at_cds.json
    {"ACT1": {"protein_accession": ..., "cds": "ATG...", "cds_length": 1134,
                 "genomic_accession": ..., "exons": [[s, e], ...],
                 "strand": "+", "n_exons": 4}, ...}

Run from pipeline/:
    python fetch_at_cds.py
    python fetch_at_cds.py --gene ACT1
    python fetch_at_cds.py --resume          # skip genes already fetched
"""

from __future__ import annotations

import argparse
import time

from Bio.Blast import NCBIWWW, NCBIXML

import ncbi
from paths import AT_CDS_FILE, load_genes, load_json, resolve_genes, save_json

EVALUE = 1e-10        # stringent: we want the Arabidopsis gene itself, not a paralog


def blast_arabidopsis(gene_name: str, protein_seq: str) -> str | None:
    """Return the top A. thaliana refseq_protein accession for a bait protein."""
    print(f"  blastp {gene_name} ({len(protein_seq)} aa) vs A. thaliana "
          f"refseq_protein ...", flush=True)
    handle = NCBIWWW.qblast(
        "blastp", "refseq_protein", protein_seq,
        entrez_query="Arabidopsis thaliana[organism]", hitlist_size=3, expect=EVALUE,
    )
    records = list(NCBIXML.parse(handle))
    if not records or not records[0].alignments:
        return None

    alignment = records[0].alignments[0]
    hsp = alignment.hsps[0]
    identity = 100.0 * hsp.identities / hsp.align_length
    print(f"    -> {alignment.accession}  identity={identity:.1f}%  "
          f"bits={hsp.bits:.0f}", flush=True)
    return alignment.accession


def process_gene(gene_name: str, protein_seq: str) -> dict:
    accession = blast_arabidopsis(gene_name, protein_seq)
    if not accession:
        raise ValueError("no BLAST hit in A. thaliana refseq_protein")

    result = ncbi.fetch_cds_from_protein(accession)
    print(f"  CDS: {result['cds_length']} bp, {result['n_exons']} exon(s), "
          f"strand={result['strand']}, accession={result['genomic_accession']}",
          flush=True)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("--gene", help="run only this gene (e.g. ACT1)")
    parser.add_argument("--resume", action="store_true",
                        help="skip genes that already have a CDS")
    args = parser.parse_args()

    ncbi.configure()
    genes = load_genes()
    results = load_json(AT_CDS_FILE)

    for gene_name in resolve_genes(args.gene):
        if args.resume and "cds" in results.get(gene_name, {}):
            print(f"[{gene_name}] already done, skipping")
            continue

        print(f"\n=== {gene_name} ===", flush=True)
        try:
            results[gene_name] = process_gene(gene_name, genes[gene_name])
        except Exception as exc:
            print(f"  ERROR: {exc}", flush=True)
            results[gene_name] = {"error": str(exc)}

        save_json(AT_CDS_FILE, results)     # checkpoint after every gene
        time.sleep(1)

    n_ok  = sum(1 for v in results.values() if "cds" in v)
    n_err = sum(1 for v in results.values() if "error" in v)
    print(f"\nDone. {AT_CDS_FILE}: {n_ok} CDS, {n_err} error(s)")


if __name__ == "__main__":
    main()
