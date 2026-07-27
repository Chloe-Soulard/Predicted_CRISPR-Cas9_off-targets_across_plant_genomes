"""
Step 5b — Recover windows CRISPOR could not find in its genome, and resubmit them.

A conserved window is derived from a transcript, so it can fail to match the
assembly CRISPOR hosts — typically because CRISPOR carries an older assembly
version than the one the ortholog is annotated on. For every window with
`found_in_genome == false` this step:

    1. BLASTs the window against that species' sequences at NCBI (megablast)
    2. picks the best genomic hit belonging to CRISPOR's assembly, verified
       against the UCSC GenArk scaffold list for that accession when available
    3. fetches the actual assembly sequence for the locus, plus flanks
    4. resubmits that genome-matched sequence to CRISPOR and parses the result

The recovered record replaces the original in results/step5/{gene}.json, carrying
its provenance so recovered windows stay distinguishable in the analysis:

    recovered_from_genome       true
    recovery_region             accession:start-end that was fetched
    recovery_in_assembly        did the accession appear in CRISPOR's assembly?
                                (null when GenArk has no scaffold list for it)
    recovery_window_ident_pct   identity of the window to the recovered locus
    recovery_sequence           the sequence actually submitted
    original_status             "not_found_in_genome"

Windows that still cannot be matched keep `recovered_from_genome: false` and a
`recovery_status` explaining why; those were excluded from the analysis.

Run from pipeline/:
    python step5b_recover_genomic.py --all --resume
    python step5b_recover_genomic.py --species "Nicotiana tabacum" --gene AtActin
"""

from __future__ import annotations

import argparse
import re
import time

import requests
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

from crispea import crispor, ncbi
from crispea.paths import (
    STEP5_DIR, STEP34_DIR, gene_slug, load_json, load_species_genomes, save_json,
)

FLANK_BP     = 25          # context either side, so guides at the window edge survive
HITLIST_SIZE = 15
EVALUE       = 1e-5
POLL_TRIES   = 10          # this step waits longer than step 5: far fewer windows
POLL_DELAY_S = 12

TRANSCRIPT_PREFIX = re.compile(r"^(NM_|XM_|XR_|NR_)")

_scaffold_cache: dict[str, set[str]] = {}


def is_genomic(accession: str) -> bool:
    """False for transcript records: only assembly sequence is useful here."""
    return not TRANSCRIPT_PREFIX.match(accession)


def ucsc_scaffolds(session: requests.Session, genome: str) -> set[str]:
    """Sequence names in a UCSC GenArk assembly; empty if it is not hosted there."""
    if genome in _scaffold_cache:
        return _scaffold_cache[genome]
    try:
        response = session.get("https://api.genome.ucsc.edu/list/chromosomes",
                               params={"genome": genome}, timeout=90)
        names = set(response.json().get("chromosomes", {}).keys())
    except Exception:
        names = set()
    _scaffold_cache[genome] = names
    return names


def unresolved_windows(genomes: dict[str, str]) -> list[dict]:
    """Every (gene, species) window CRISPOR did not find in its genome."""
    out = []
    for path in sorted(STEP5_DIR.glob("*.json")):
        slug = path.stem
        step34 = load_json(STEP34_DIR / f"{slug}.json")
        windows = step34.get("exonic_windows", {})
        for species, record in load_json(path).items():
            if record.get("found_in_genome") is False:
                out.append({
                    "gene_slug": slug,
                    "gene": step34.get("gene", slug),
                    "species": species,
                    "genome": record.get("genome") or genomes.get(species),
                    "window": (windows.get(species, {}).get("sequence") or "").upper(),
                })
    return out


def blast_window(sequence: str, species: str) -> list[dict]:
    """megablast a window against NCBI nt, restricted to one organism."""
    handle = NCBIWWW.qblast("blastn", "nt", sequence,
                            entrez_query=f"{species}[organism]", megablast=True,
                            hitlist_size=HITLIST_SIZE, expect=EVALUE)
    record = NCBIXML.read(handle)
    handle.close()
    return [{"accession": alignment.accession,
             "subject_start": alignment.hsps[0].sbjct_start,
             "subject_end": alignment.hsps[0].sbjct_end,
             "identities": alignment.hsps[0].identities,
             "align_length": alignment.hsps[0].align_length,
             "query_start": alignment.hsps[0].query_start,
             "query_end": alignment.hsps[0].query_end}
            for alignment in record.alignments]


def fetch_region(accession: str, low: int, high: int) -> str:
    ncbi.pause()
    with Entrez.efetch(db="nuccore", id=accession, rettype="fasta", retmode="text",
                       seq_start=low, seq_stop=high) as handle:
        return str(SeqIO.read(handle, "fasta").seq).upper()


def recover_sequence(session: requests.Session, window: dict) -> tuple[str | None, dict]:
    """Locate the window in its assembly and fetch the real genomic sequence."""
    scaffolds = ucsc_scaffolds(session, window["genome"])   # may be empty
    hits = blast_window(window["window"], window["species"])

    chosen = None
    for hit in hits:
        if not is_genomic(hit["accession"]):
            continue
        if scaffolds and hit["accession"] in scaffolds:
            chosen = {**hit, "in_assembly": True}          # best case: same assembly
            break
        if chosen is None:
            chosen = {**hit, "in_assembly": False if scaffolds else None}

    if chosen is None:
        return None, {"recovery_status": "no_genomic_hit",
                      "blast_hits": [h["accession"] for h in hits[:5]]}

    # Extend the hit back out to the full window, then add flanks.
    low, high = sorted((chosen["subject_start"], chosen["subject_end"]))
    low = max(1, low - (chosen["query_start"] - 1) - FLANK_BP)
    high = high + (len(window["window"]) - chosen["query_end"]) + FLANK_BP

    sequence = fetch_region(chosen["accession"], low, high)
    return sequence, {
        "recovery_status": "ok",
        "recovery_region": f"{chosen['accession']}:{low}-{high}",
        "recovery_in_assembly": chosen["in_assembly"],
        "recovery_window_ident_pct": round(
            100.0 * chosen["identities"] / chosen["align_length"], 1),
    }


def run_crispor(session, sequence: str, genome: str, name: str):
    """Submit a recovered sequence and wait for the result; (None, None) if queued."""
    bid, _ = crispor.submit(session, sequence, genome, name)
    link = crispor.result_url(bid)
    for _ in range(POLL_TRIES):
        if crispor.is_running(crispor.fetch_html(session, bid)):
            time.sleep(POLL_DELAY_S)
            continue
        guides = crispor.download_tsv(session, bid, "guides")
        if guides is not None:
            return bid, link, guides, crispor.download_tsv(session, bid, "offtargets") or ""
        time.sleep(POLL_DELAY_S)
    return bid, link, None, None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("--all", action="store_true",
                       help="every unresolved window (required unless --gene/--species narrows it)")
    parser.add_argument("--gene", help="restrict to one gene")
    parser.add_argument("--species", help="restrict to one species")
    parser.add_argument("--resume", action="store_true",
                        help="skip windows already recovered successfully")
    parser.add_argument("--limit", type=int, default=None,
                        help="cap the windows processed this run")
    args = parser.parse_args()

    # This step submits to CRISPOR, so make the full run an explicit choice.
    if not (args.all or args.gene or args.species):
        parser.error("pass --all, or narrow the run with --gene / --species")

    ncbi.configure()
    genomes = load_species_genomes()
    session = crispor.session(ncbi.contact_email())

    windows = unresolved_windows(genomes)
    if args.species:
        windows = [w for w in windows if w["species"] == args.species]
    if args.gene:
        windows = [w for w in windows if w["gene"] == args.gene]
    if args.limit:
        windows = windows[:args.limit]

    print(f"{len(windows)} window(s) to recover and resubmit", flush=True)

    with crispor.submission_lock(STEP5_DIR):
        for window in windows:
            species = window["species"]
            out_file = STEP5_DIR / f"{window['gene_slug']}.json"
            results = load_json(out_file)
            previous = results.get(species, {})

            if (args.resume and previous.get("recovered_from_genome")
                    and previous.get("status") == "ok"):
                print(f"  [{window['gene']} / {species}] already recovered — skip",
                      flush=True)
                continue

            print(f"  [{window['gene']} / {species} @ {window['genome']}] "
                  f"BLAST against the assembly ...", flush=True)
            try:
                sequence, meta = recover_sequence(session, window)
            except Exception as exc:
                print(f"     recovery error: {exc}", flush=True)
                continue

            if sequence is None:
                print(f"     {meta['recovery_status']} — cannot recover", flush=True)
                results[species] = {**previous, **meta, "recovered_from_genome": False}
                save_json(out_file, results)
                continue

            print(f"     recovered {meta['recovery_region']} "
                  f"(window identity {meta['recovery_window_ident_pct']}%, "
                  f"in_assembly={meta['recovery_in_assembly']}); submitting ...",
                  flush=True)

            bid, link, guides, offtargets = run_crispor(
                session, sequence, window["genome"],
                f"{window['gene']}_{species}_recovered")

            if guides is None:
                results[species] = {**previous, **meta, "recovered_from_genome": True,
                                    "status": "pending", "batch_id": bid,
                                    "crispor_link": link, "recovery_sequence": sequence}
                print(f"     pending (queued) -> {link}", flush=True)
            else:
                record = crispor.build_record(
                    guides, offtargets, genome=window["genome"], bid=bid, link=link,
                    window_source={"window_source": "recovered_genomic"})
                record.update(meta)
                record["recovered_from_genome"] = True
                record["original_status"] = "not_found_in_genome"
                record["recovery_sequence"] = sequence
                results[species] = record
                print(f"     OK: found={record['found_in_genome']} "
                      f"guides={record['n_guides']} "
                      f"off-targets={record['ot_total_all_guides']}", flush=True)

            save_json(out_file, results)
            time.sleep(2)

    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
