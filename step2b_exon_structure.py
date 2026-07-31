"""
Step 2b — Add real coding-exon boundaries, in CDS coordinates, to each species.

Step 2's `coded_by` accession is usually an mRNA, so it carries no genomic exon
boundaries. Guide windows must not straddle an intron, so step 3/4 needs to know
where the exon junctions fall inside each CDS. Per species:

    1. protein accession  --elink(protein -> gene)-->  GeneID
    2. efetch(db=gene, rettype=gene_table)  ->  per-transcript exon table
    3. sum the coding length of each exon of our transcript; the cumulative sums
       (excluding the last) are the junction offsets in CDS coordinates

Fields added to each results/step2/{gene}.json species entry:
    cds_junctions  CDS offsets at which each non-first coding exon begins
    n_exons        real coding-exon count (replacing step 2's interval count)
    exon_source    ncbi_gene_table | step2_exon_lengths | assumed_single_exon | error
    exon_warning   set when the coding lengths did not sum to cds_length

Two fallbacks keep species in the analysis when NCBI Gene has no record — usually
proteins annotated on GenBank/WGS assemblies:

    step2_exon_lengths    derive junctions from the `coded_by` intervals step 2
                          already parsed, accepted only if their lengths sum to
                          the CDS length
    assumed_single_exon   curated species, which have no NCBI gene model at all

Both are marked so step 3/4 can flag their windows as exon-unverified.

Run from pipeline/:
    python step2b_exon_structure.py --gene ACT1
    python step2b_exon_structure.py --all
    python step2b_exon_structure.py --all --resume
"""

from __future__ import annotations

import argparse

import ncbi
from paths import STEP2_DIR, gene_slug, load_json, resolve_genes, save_json


def junctions_from_coded_by(entry: dict) -> dict | None:
    """Derive coding-exon junctions from step 2's `coded_by` intervals.

    Those intervals are the coding exon spans, so their lengths sum to the CDS
    length. A minus-strand transcript reads the exons in reverse genomic order,
    so the lengths are reversed before accumulating.
    """
    exons = entry.get("exons")
    if not exons:
        return None
    lengths = [end - start + 1 for start, end in exons]
    if entry.get("strand") == "-":
        lengths.reverse()
    return {"cds_junctions": ncbi.cumulative_junctions(lengths),
            "n_exons": len(lengths),
            "total_coding": sum(lengths)}


def exon_structure_from_ncbi(protein_accession: str, transcript: str) -> dict:
    """Coding-exon structure for one species, or {'exon_error': reason}."""
    try:
        gene_id = ncbi.protein_to_geneid(protein_accession)
    except Exception as exc:
        return {"exon_error": f"elink failed: {exc}"}
    if not gene_id:
        return {"exon_error": "no GeneID from elink"}

    try:
        table = ncbi.fetch_gene_table(gene_id)
    except Exception as exc:
        return {"exon_error": f"gene_table fetch failed: {exc}"}

    parsed = ncbi.parse_gene_table(table, transcript)
    if parsed is None:
        return {"exon_error": f"transcript {transcript} not found or unparseable "
                              f"in gene_table", "gene_id": gene_id}
    return {**parsed, "gene_id": gene_id}


def apply_structure(entry: dict, info: dict) -> str:
    """Write an exon structure onto a species entry; return a one-line status."""
    entry["cds_junctions"] = info["cds_junctions"]
    entry["n_exons"] = info["n_exons"]
    entry["exon_source"] = "ncbi_gene_table"
    entry["gene_id"] = info["gene_id"]
    entry.pop("exon_error", None)

    # The coding lengths must reconcile with the CDS we actually reconstructed.
    # When they do not (alternative transcript, annotation update), the junctions
    # cannot be trusted, so fall back to treating the CDS as one exon.
    if info["total_coding"] != entry.get("cds_length"):
        entry["exon_warning"] = (f"coding sum {info['total_coding']} != cds_length "
                                 f"{entry.get('cds_length')}; treating as single exon")
        entry["cds_junctions"] = []
        entry["n_exons"] = 1
        return f"WARN: {entry['exon_warning']}"

    entry.pop("exon_warning", None)
    return f"OK: {info['n_exons']} coding exon(s), junctions={info['cds_junctions']}"


def run_gene(gene_name: str, resume: bool = False) -> None:
    out_file = STEP2_DIR / f"{gene_slug(gene_name)}.json"
    if not out_file.exists():
        print(f"[{gene_name}] no step2 file — run step2_fetch_genomic.py first")
        return

    data = load_json(out_file)
    print(f"\n{'=' * 68}\nGene: {gene_name} | exon structure for {len(data)} species"
          f"\n{'=' * 68}", flush=True)

    n_ok = n_curated = n_error = 0
    for species, entry in data.items():
        if "cds" not in entry:                              # no CDS to annotate
            continue
        if resume and "cds_junctions" in entry and "exon_error" not in entry:
            n_ok += 1
            continue

        # Curated species are not in NCBI Gene at all; treat as a single exon and
        # mark them so step 3/4 flags their windows as unverified.
        if entry.get("source") == "manual":
            entry["cds_junctions"] = []
            entry["n_exons"] = 1
            entry["exon_source"] = "assumed_single_exon"
            n_curated += 1
            print(f"  [{species}] curated -> assumed single exon", flush=True)
            save_json(out_file, data)
            continue

        protein    = entry.get("protein_accession", "")
        transcript = entry.get("genomic_accession", "")
        print(f"  [{species}] {protein} / {transcript} ...", flush=True)
        info = exon_structure_from_ncbi(protein, transcript)

        if "exon_error" not in info:
            print(f"      {apply_structure(entry, info)}", flush=True)
            n_ok += 1
        else:
            fallback = junctions_from_coded_by(entry)
            if fallback and fallback["total_coding"] == entry.get("cds_length"):
                entry["cds_junctions"] = fallback["cds_junctions"]
                entry["n_exons"] = fallback["n_exons"]
                entry["exon_source"] = "step2_exon_lengths"
                entry.pop("exon_error", None)
                entry.pop("exon_warning", None)
                n_ok += 1
                print(f"      FALLBACK ({info['exon_error']}): {fallback['n_exons']} "
                      f"exon(s) from coded_by, junctions={fallback['cds_junctions']}",
                      flush=True)
            else:
                entry["exon_error"] = info["exon_error"]
                entry.setdefault("cds_junctions", [])
                entry.setdefault("n_exons", 1)
                entry["exon_source"] = "error"
                n_error += 1
                print(f"      ERROR: {info['exon_error']}", flush=True)

        save_json(out_file, data)         # checkpoint after every species

    print(f"\n[{gene_name}] {n_ok} resolved, {n_curated} curated, {n_error} error(s)",
          flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--gene", help="single gene (e.g. ACT1)")
    group.add_argument("--all", action="store_true", help="run every gene in the panel")
    parser.add_argument("--resume", action="store_true",
                        help="skip species that already have exon junctions")
    args = parser.parse_args()

    ncbi.configure()
    for gene_name in resolve_genes(args.gene):
        run_gene(gene_name, resume=args.resume)

    print("\nAll done.")


if __name__ == "__main__":
    main()
