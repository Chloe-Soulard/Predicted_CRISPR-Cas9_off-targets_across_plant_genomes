"""CRISPOR web-service client and off-target record building.

CRISPOR derives its batch identifier deterministically from the submitted
sequence, genome, PAM and job name, so a repeated submission returns the cached
result rather than queueing a new job. The pipeline relies on this: it submits
once, records the batch id, and fetches the result on a later pass.

Shared by `step5_crispor.py` (conserved windows) and
`step5b_recover_genomic.py` (assembly-matched replacement windows).
"""

from __future__ import annotations

import base64
import hashlib
import os
import re
import time
from contextlib import contextmanager
from pathlib import Path

import requests

BASE_URL = "https://crispor.gi.ucsc.edu/crispor.py"   # crispor.tefor.net redirects here
PAM      = "NGG"                                       # Streptococcus pyogenes Cas9
VERSION  = "5.2"                                       # CRISPOR version used in the study

# CRISPOR PAM codes this pipeline has been used with. The PAM is part of the batch
# id, so results for different nucleases never collide.
PAM_NAMES = {
    "NGG":  "SpCas9, 20 bp guides",
    "TTTV": "Cas12a (Cpf1), 23 bp guides",
    "NAA":  "iSpyMacCas9, 20 bp guides",
}

# MIT and CFD specificity are Cas9 scores. CRISPOR returns -1 for every guide of
# these 5'-PAM nucleases (Cas12a/Cas12b/CasX families), so a guide's presence in
# the genome cannot be read from them. Cas9-type PAMs — including variants such as
# NAA (iSpyMacCas9) — are scored normally and are the default.
UNSCORED_PAMS = frozenset({
    "TTTV", "TTTV-21", "TTTN", "TTN", "TNN", "ATTN", "TTCN", "NGTN",
    "TYCV", "TATV", "TTTA", "TCTA", "TCCA", "CCCA", "GGTT", "YTTV", "TTYN",
})

REQUEST_TIMEOUT_S = 180


def session(contact_email: str) -> requests.Session:
    """A requests session that identifies the study to the CRISPOR maintainers."""
    s = requests.Session()
    s.headers.update({
        "User-Agent": f"CRISPEA-offtarget-study/1.0 (research; {contact_email})"
    })
    return s


def job_name(gene: str, species: str) -> str:
    """Stable CRISPOR job name; part of the deterministic batch id."""
    return f"{gene}_{re.sub(r'[^A-Za-z0-9]', '_', species)}"


def batch_id(seq: str, genome: str, name: str, pam: str = PAM) -> str:
    """Reproduce CRISPOR's own id: base64.urlsafe(sha1(seq+org+pam+name))[:20]."""
    digest = hashlib.sha1((seq + genome + pam + name).encode("latin1")).digest()[:20]
    return base64.urlsafe_b64encode(digest).decode("latin1")[:20]


def result_url(bid: str) -> str:
    return f"{BASE_URL}?batchId={bid}"


# ── HTTP ──────────────────────────────────────────────────────────────────────

def submit(session: requests.Session, seq: str, genome: str, name: str,
           pam: str = PAM) -> tuple[str, str]:
    """Submit a sequence; return (batch id, response HTML).

    The batch id is read back from the response when possible and recomputed
    locally otherwise — both routes give the same value.
    """
    response = session.get(
        BASE_URL,
        params={"seq": seq, "org": genome, "pam": pam, "name": name},
        timeout=REQUEST_TIMEOUT_S, allow_redirects=True,
    )
    match = re.search(r"batchId=([A-Za-z0-9\-_]+)", response.url + " " + response.text)
    bid = match.group(1) if match else batch_id(seq, genome, name, pam)
    return bid, response.text


def fetch_html(session: requests.Session, bid: str) -> str:
    return session.get(BASE_URL, params={"batchId": bid},
                       timeout=REQUEST_TIMEOUT_S).text


def download_tsv(session: requests.Session, bid: str, kind: str) -> str | None:
    """Download the 'guides' or 'offtargets' TSV; None if not ready or errored."""
    response = session.get(
        BASE_URL, params={"batchId": bid, "download": kind, "format": "tsv"},
        timeout=REQUEST_TIMEOUT_S,
    )
    if response.status_code != 200:
        return None
    text = response.text
    if "lead to an error" in text or "<html" in text[:200].lower():
        return None
    return text


def is_running(html: str) -> bool:
    """True while the job is still queued or computing."""
    return "Job Status" in html or "has been submitted" in html


def looks_not_in_genome(html: str) -> bool:
    """True when CRISPOR reports the submitted sequence is absent from the genome."""
    low = html.lower()
    return any(s in low for s in (
        "not present in the genome", "not found in the genome",
        "could not be found in the genome", "is not on the genome",
        "no match in the genome",
    ))


# ── TSV parsing ───────────────────────────────────────────────────────────────

# CRISPOR's off-target download omits the header row for some nucleases (Cas12a
# among them), so the first line is already data. Without this fallback the first
# off-target is mistaken for the header and every row is keyed by garbage, which
# silently zeroes the mismatch and genic/intergenic tallies.
OFFTARGET_COLUMNS = [
    "guideId", "guideSeq", "offtargetSeq", "mismatchPos", "mismatchCount",
    "mitOfftargetScore", "cfdOfftargetScore", "chrom", "start", "end",
    "strand", "locusDesc",
]


def parse_tsv(text: str, default_header: list[str] | None = None) -> list[dict]:
    """Parse a CRISPOR TSV download into a list of row dicts.

    `default_header` is used when the download carries no header line, detected by
    the first field not naming a column.
    """
    lines = [ln for ln in (text or "").splitlines() if ln.strip()]
    if not lines:
        return []

    first = lines[0].lstrip("#").split("\t")
    if default_header is not None and first[0].strip() != default_header[0]:
        header, body = default_header, lines
    else:
        header, body = first, lines[1:]

    rows = []
    for line in body:
        cells = line.split("\t")
        cells += [""] * (len(header) - len(cells))
        rows.append(dict(zip(header, cells)))
    return rows


def _number(value) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def classify_locus(locus_desc: str) -> str:
    """Classify a CRISPOR off-target from its `locusDesc` annotation.

    'exon:...', 'intron:...' and bare gene names are genic; anything beginning
    'intergenic' is intergenic; an empty annotation is unannotated.
    """
    desc = (locus_desc or "").strip().lower()
    if not desc:
        return "unannotated"
    return "intergenic" if desc.startswith("intergenic") else "genic"


def build_record(guides_tsv: str, offtargets_tsv: str, *, genome: str, bid: str,
                 link: str, window_source: dict, pam: str = PAM) -> dict:
    """Turn the two CRISPOR TSVs into one window record.

    Off-targets are grouped per guide and tallied by mismatch count (0-4) and by
    locus class. CRISPOR sets a specificity of -1 for a guide whose sequence is
    absent from the selected genome; such guides have no meaningful off-targets
    and are excluded from the window totals.

    The MIT and CFD specificity scores are defined for SpCas9 only. For any other
    nuclease CRISPOR returns -1 for every guide, so that test would discard the
    whole window; there, presence in the genome is taken from the window itself,
    which `looks_not_in_genome` already establishes before this is called.
    `in_genome_basis` records which test was used.
    """
    scored = pam.strip().upper() not in UNSCORED_PAMS
    guide_rows = parse_tsv(guides_tsv)

    offtargets_by_guide: dict[str, list[dict]] = {}
    for row in parse_tsv(offtargets_tsv, OFFTARGET_COLUMNS):
        offtargets_by_guide.setdefault(row.get("guideId", ""), []).append(row)

    guides = []
    for row in guide_rows:
        guide_id = row.get("guideId", "")
        offtargets = offtargets_by_guide.get(guide_id, [])

        by_mismatch = {str(k): 0 for k in range(5)}
        genic = intergenic = unannotated = 0
        for offtarget in offtargets:
            mismatches = offtarget.get("mismatchCount", "")
            if mismatches in by_mismatch:
                by_mismatch[mismatches] += 1
            locus_class = classify_locus(offtarget.get("locusDesc", ""))
            if locus_class == "genic":
                genic += 1
            elif locus_class == "intergenic":
                intergenic += 1
            else:
                unannotated += 1

        mit = _number(row.get("mitSpecScore"))
        cfd = _number(row.get("cfdSpecScore"))
        locus = (row.get("targetGenomeGeneLocus") or "").strip()
        count = int(_number(row.get("offtargetCount")) or 0)

        # "This guide is a repeated region, it is too unspecific" — CRISPOR gives up
        # rather than enumerate, and reports zero off-targets with a zero specificity
        # and no locus. Those zeros are unknowns, not clean guides, so they are
        # flagged and kept out of the totals; counting them as 0 would credit the
        # most repetitive guides in the most repetitive genomes as the most specific.
        repeat = count == 0 and not locus and (
            (mit == 0 and cfd == 0) if scored else True)

        guides.append({
            "guide_id":            guide_id,
            "target_seq":          row.get("targetSeq", ""),
            "on_target_locus":     row.get("targetGenomeGeneLocus", ""),
            "repeat_unspecific":   repeat,
            "in_genome":           (mit is not None and mit != -1) if scored else True,
            "mit_spec":            mit,
            "cfd_spec":            _number(row.get("cfdSpecScore")),
            "offtarget_count":     int(_number(row.get("offtargetCount")) or 0),
            "offtargets_in_table": len(offtargets),
            "ot_by_mismatch":      by_mismatch,
            "ot_genic":            genic,
            "ot_intergenic":       intergenic,
            "ot_unannotated":      unannotated,
        })

    n_in  = sum(1 for g in guides if g["in_genome"])
    n_out = len(guides) - n_in
    if not guides:
        found_state = "no_guides"
    elif n_out == 0:
        found_state = "found"       # every guide maps: the window is genomic DNA
    elif n_in == 0:
        found_state = "not_found"   # the window is absent from this assembly
    else:
        found_state = "partial"     # part of the window maps (e.g. spans a junction)

    record = {
        "genome": genome,
        "pam": pam,
        "in_genome_basis": "mit_spec" if scored else "window",
        "batch_id": bid,
        "crispor_link": link,
        "status": "ok",
        "found_state": found_state,
        "found_in_genome": found_state == "found",
        "needs_genomic_dna": found_state in ("partial", "not_found"),
        "n_guides": len(guides),
        "n_guides_in_genome": n_in,
        "n_guides_not_in_genome": n_out,
        "n_guides_repeat": sum(1 for g in guides if g["repeat_unspecific"]),
        # Totals count only guides that are actually present in the genome.
        "ot_total_all_guides": sum(g["offtarget_count"] for g in guides if g["in_genome"]),
        "guides": guides,
    }
    record.update(window_source)
    return record


# ── Single-writer lock ────────────────────────────────────────────────────────
#
# Two concurrent CRISPOR passes would interleave writes to the same
# results/step5/<gene>.json and lose records. `step5_crispor.py` takes this lock
# for the length of a run and `step5_loop.py` waits on it, so an unattended loop
# and a manual run never clash.

STALE_LOCK_S = 6 * 3600


def lock_path(directory: Path) -> Path:
    return directory / ".crispor.lock"


def lock_holder(directory: Path) -> int | None:
    """PID currently holding the lock, or None if free (or the lock is stale)."""
    path = lock_path(directory)
    try:
        age = time.time() - path.stat().st_mtime
    except FileNotFoundError:
        return None
    if age > STALE_LOCK_S:
        path.unlink(missing_ok=True)
        return None
    try:
        return int(path.read_text(encoding="utf-8").split()[0])
    except (ValueError, IndexError, OSError):
        return -1


@contextmanager
def submission_lock(directory: Path):
    """Hold the CRISPOR write lock for `directory`, or exit if already held."""
    directory.mkdir(parents=True, exist_ok=True)
    path = lock_path(directory)
    try:
        fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        holder = lock_holder(directory)
        if holder is None:                       # released or expired mid-check
            fd = os.open(path, os.O_CREAT | os.O_WRONLY)
        else:
            raise SystemExit(
                f"Another CRISPOR run (pid {holder}) is writing to {directory}.\n"
                f"Wait for it to finish, or delete {path} if it crashed."
            )
    try:
        with os.fdopen(fd, "w") as handle:
            handle.write(f"{os.getpid()} {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        yield
    finally:
        path.unlink(missing_ok=True)
