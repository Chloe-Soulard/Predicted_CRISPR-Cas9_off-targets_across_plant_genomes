"""
Step 5 — Count predicted off-target sites with CRISPOR.

Each conserved window from step 3/4 is submitted to CRISPOR against that species'
genome assembly, using SpCas9 (20 nt spacer, NGG PAM), and the result recorded per
guide: MIT and CFD specificity, total off-targets up to 4 mismatches, the mismatch
distribution, and genic / intergenic / unannotated counts from CRISPOR's
`locusDesc` annotation.

CRISPOR's batch id is a deterministic hash of sequence + genome + PAM + job name,
so a repeat submission returns the cached result instead of queueing a new job.
This step therefore submits once, records the batch id, and takes the result if it
is ready; otherwise the window is marked "pending" and picked up by the next
--resume pass. Nothing is ever resubmitted, and the shared server is never blocked
on for long. `step5_loop.py` drives the passes unattended.

Output: results/step5/{gene_slug}.json
    {"Zea mays": {"status": "ok", "found_state": "found", "found_in_genome": true,
                  "n_guides": 8, "ot_total_all_guides": 1483,
                  "batch_id": ..., "crispor_link": ..., "guides": [...]},
     "Other species": {"status": "pending", ...},
     ...}

`found_state` distinguishes a window that fully maps to the assembly ("found")
from one that is absent ("not_found") or only partly present ("partial"); the
latter two are handled by step5b_recover_genomic.py.

Run from pipeline/:
    python step5_crispor.py --gene AtActin
    python step5_crispor.py --all --resume
    python step5_crispor.py --gene AtActin --max 3     # cap windows, for a smoke test
"""

from __future__ import annotations

import argparse
import time

import crispor
import ncbi
from paths import (
    STEP34_DIR, STEP5_DIR, gene_slug, load_json, load_species_genomes,
    resolve_genes, save_json,
)

SUBMIT_DELAY_S = 3.0     # polite pause between windows
POLL_TRIES     = 2       # brief check for a fresh result; defer if still queued
POLL_DELAY_S   = 6.0


def fetch_result(session, bid: str) -> tuple[str | None, str]:
    """Poll briefly for a finished job. Returns (guides TSV or None, last HTML)."""
    html = ""
    for _ in range(POLL_TRIES):
        html = crispor.fetch_html(session, bid)
        if crispor.is_running(html):
            time.sleep(POLL_DELAY_S)
            continue
        if crispor.looks_not_in_genome(html):
            break
        guides_tsv = crispor.download_tsv(session, bid, "guides")
        if guides_tsv is not None:
            return guides_tsv, html
        time.sleep(POLL_DELAY_S)
    return None, html


def process_gene(gene: str, genomes: dict[str, str], session,
                 resume: bool, limit: int | None) -> None:
    slug = gene_slug(gene)
    step34_file = STEP34_DIR / f"{slug}.json"
    if not step34_file.exists():
        print(f"[{gene}] no step34 file — run step34_conserved.py first", flush=True)
        return

    windows = load_json(step34_file).get("exonic_windows", {})
    out_file = STEP5_DIR / f"{slug}.json"
    results = load_json(out_file)

    print(f"\n{'=' * 68}\n[{gene}] {len(windows)} windows -> CRISPOR\n{'=' * 68}",
          flush=True)

    n_done = n_pending = n_skipped = 0
    for species, window in windows.items():
        if limit is not None and n_done + n_pending >= limit:
            break

        previous = results.get(species)
        if resume and previous and previous.get("status") == "ok":
            n_skipped += 1
            continue

        sequence = (window.get("sequence") or "").upper()
        genome = genomes.get(species)
        if not sequence or not genome:
            results[species] = {"status": "error", "genome": genome,
                                "reason": "missing window sequence or CRISPOR genome"}
            save_json(out_file, results)
            continue

        provenance = {"window_source": window.get("source")}
        if window.get("exon_unverified"):
            provenance["exon_unverified"] = True

        name = crispor.job_name(gene, species)
        bid = (previous or {}).get("batch_id") or crispor.batch_id(sequence, genome, name)
        link = crispor.result_url(bid)

        # Submit only if this window has never been submitted (no batch id on file).
        if not (previous and previous.get("batch_id")):
            try:
                bid, _ = crispor.submit(session, sequence, genome, name)
                link = crispor.result_url(bid)
            except Exception as exc:
                results[species] = {"status": "error", "reason": f"submit failed: {exc}",
                                    "genome": genome, "batch_id": bid, "crispor_link": link}
                save_json(out_file, results)
                time.sleep(SUBMIT_DELAY_S)
                continue

        guides_tsv, html = fetch_result(session, bid)

        if guides_tsv is None:
            if crispor.looks_not_in_genome(html):
                results[species] = {"status": "ok", "found_in_genome": False,
                                    "found_state": "not_found", "needs_genomic_dna": True,
                                    "n_guides": 0, "guides": [], "genome": genome,
                                    "batch_id": bid, "crispor_link": link, **provenance}
                print(f"  [{species}] NOT in genome", flush=True)
                n_done += 1
            else:
                results[species] = {"status": "pending", "genome": genome,
                                    "batch_id": bid, "crispor_link": link, **provenance}
                print(f"  [{species}] pending (queued)", flush=True)
                n_pending += 1
            save_json(out_file, results)
            time.sleep(SUBMIT_DELAY_S)
            continue

        offtargets_tsv = crispor.download_tsv(session, bid, "offtargets") or ""
        record = crispor.build_record(guides_tsv, offtargets_tsv, genome=genome,
                                      bid=bid, link=link, window_source=provenance)
        results[species] = record
        save_json(out_file, results)
        print(f"  [{species}] OK: {record['n_guides']} guides, "
              f"{record['ot_total_all_guides']} off-targets "
              f"(found={record['found_in_genome']})", flush=True)
        n_done += 1
        time.sleep(SUBMIT_DELAY_S)

    print(f"\n[{gene}] this pass: {n_done} fetched, {n_pending} pending, "
          f"{n_skipped} already done", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--gene", help="single gene (e.g. AtActin)")
    group.add_argument("--all", action="store_true",
                       help="every gene that has a step34 window file")
    parser.add_argument("--resume", action="store_true",
                        help="skip windows already 'ok'; fetch pending ones, retry errors")
    parser.add_argument("--max", type=int, default=None,
                        help="cap the windows processed this run (for testing)")
    args = parser.parse_args()

    STEP5_DIR.mkdir(parents=True, exist_ok=True)
    genomes = load_species_genomes()

    genes = [g for g in resolve_genes(args.gene)
             if (STEP34_DIR / f"{gene_slug(g)}.json").exists()]
    if not genes:
        raise SystemExit("No step34 window files found — run step34_conserved.py first.")

    session = crispor.session(ncbi.contact_email())

    # One writer at a time: concurrent passes would interleave writes to the same
    # per-gene JSON and lose records.
    with crispor.submission_lock(STEP5_DIR):
        for gene in genes:
            process_gene(gene, genomes, session, resume=args.resume, limit=args.max)

    print("\nAll done.")


if __name__ == "__main__":
    main()
