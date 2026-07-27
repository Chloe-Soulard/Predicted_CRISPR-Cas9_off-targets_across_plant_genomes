"""NCBI Entrez access: credentials, CDS reconstruction and exon structure.

Credentials come from the environment, never from the source:

    NCBI_EMAIL    required — NCBI asks every automated client to identify itself
    NCBI_API_KEY  optional — raises the rate limit from 3 to 10 requests/second

`configure()` must be called once before any Entrez request; every step that
touches NCBI does so in `main()`.
"""

from __future__ import annotations

import os
import re
import time

from Bio import Entrez, SeqIO

# Delay inserted before each Entrez request. NCBI allows 3 requests/second
# without an API key and 10 with one; both values leave headroom.
_DELAY_WITHOUT_KEY = 0.4
_DELAY_WITH_KEY    = 0.15

_delay = _DELAY_WITHOUT_KEY
_configured = False


def configure() -> None:
    """Read NCBI credentials from the environment. Idempotent."""
    global _configured, _delay
    if _configured:
        return

    email = os.environ.get("NCBI_EMAIL", "").strip()
    if not email:
        raise SystemExit(
            "NCBI_EMAIL is not set.\n"
            "NCBI requires automated clients to identify themselves. Set it once:\n"
            '  PowerShell:  $env:NCBI_EMAIL = "you@example.org"\n'
            '  bash/zsh:    export NCBI_EMAIL="you@example.org"\n'
            "Optionally also set NCBI_API_KEY (https://www.ncbi.nlm.nih.gov/account/settings/)\n"
            "to raise the rate limit from 3 to 10 requests/second."
        )
    Entrez.email = email

    api_key = os.environ.get("NCBI_API_KEY", "").strip()
    if api_key:
        Entrez.api_key = api_key
        _delay = _DELAY_WITH_KEY

    _configured = True


def contact_email() -> str:
    """The configured contact address, reused in the CRISPOR User-Agent."""
    configure()
    return Entrez.email


def pause() -> None:
    """Wait out the polite inter-request delay."""
    time.sleep(_delay)


# ── coded_by parsing ──────────────────────────────────────────────────────────

def parse_coded_by(coded_by: str) -> dict | None:
    """Parse a GenPept `coded_by` qualifier into its source interval(s).

    Handles the simple, partial (`<`/`>`), `join(...)`, `complement(...)` and
    whitespace-wrapped multi-line forms:

        NM_001036427.3:59..1192
        join(NC_003074.8:2618832..2619048,NC_003074.8:2619136..2619278)
        complement(join(NC_003075.7:100..200,NC_003075.7:300..400))

    Returns {'accession', 'exons': [[start, end], ...], 'strand'} in ascending
    genomic order (as NCBI writes them), or None if the qualifier is unparseable.
    """
    coded_by = re.sub(r"\s+", "", coded_by)          # undo GenBank line wrapping
    strand = "-" if coded_by.startswith("complement(") else "+"

    if "join(" in coded_by:
        intervals = re.findall(r"([\w.]+):<?(\d+)\.\.>?(\d+)", coded_by)
        if not intervals:
            return None
        return {"accession": intervals[0][0],
                "exons": [[int(s), int(e)] for _, s, e in intervals],
                "strand": strand}

    m = re.match(r"(?:complement\()?([\w.]+):<?(\d+)\.\.>?(\d+)\)?", coded_by)
    if not m:
        return None
    return {"accession": m.group(1),
            "exons": [[int(m.group(2)), int(m.group(3))]],
            "strand": strand}


def get_coded_by(protein_accession: str) -> str | None:
    """Fetch a GenPept record and return its CDS `coded_by` qualifier."""
    configure()
    pause()
    with Entrez.efetch(db="protein", id=protein_accession,
                       rettype="gp", retmode="text") as handle:
        record = SeqIO.read(handle, "genbank")

    for feature in record.features:
        if feature.type == "CDS":
            return feature.qualifiers.get("coded_by", [None])[0]
    return None


def fetch_interval(accession: str, start: int, end: int, strand: str) -> str:
    """Fetch one nucleotide interval, reverse-complemented for the minus strand."""
    configure()
    with Entrez.efetch(db="nucleotide", id=accession, rettype="fasta",
                       retmode="text", seq_start=start, seq_stop=end,
                       strand=2 if strand == "-" else 1) as handle:
        record = SeqIO.read(handle, "fasta")
    return str(record.seq).upper()


def splice_cds(parsed: dict) -> str:
    """Fetch and concatenate the exons of a parsed `coded_by` into a spliced CDS.

    NCBI writes `complement(join(...))` exons in ascending genomic order, but a
    minus-strand transcript reads them in the opposite order. Each exon is already
    individually reverse-complemented by `fetch_interval`, so only the
    concatenation order needs reversing — without this the reading frame is
    scrambled for every minus-strand gene.
    """
    pieces = []
    for start, end in parsed["exons"]:
        pause()
        pieces.append(fetch_interval(parsed["accession"], start, end, parsed["strand"]))
    if parsed["strand"] == "-":
        pieces.reverse()
    return "".join(pieces)


def fetch_cds_from_protein(protein_accession: str) -> dict:
    """Reconstruct the spliced CDS behind a protein accession.

    Raises ValueError if the record carries no usable `coded_by` qualifier.
    """
    coded_by = get_coded_by(protein_accession)
    if not coded_by:
        raise ValueError(f"No coded_by qualifier on {protein_accession}")

    parsed = parse_coded_by(coded_by)
    if not parsed:
        raise ValueError(f"Cannot parse coded_by: {coded_by!r}")

    cds = splice_cds(parsed)
    return {
        "protein_accession": protein_accession,
        "cds": cds,
        "cds_length": len(cds),
        "genomic_accession": parsed["accession"],
        "exons": parsed["exons"],
        "strand": parsed["strand"],
        "n_exons": len(parsed["exons"]),
    }


# ── Exon structure (NCBI Gene) ────────────────────────────────────────────────

def protein_to_geneid(protein_accession: str) -> str | None:
    """Follow the protein -> gene link to an NCBI GeneID."""
    configure()
    handle = Entrez.elink(dbfrom="protein", db="gene", id=protein_accession)
    record = Entrez.read(handle)
    handle.close()
    try:
        return record[0]["LinkSetDb"][0]["Link"][0]["Id"]
    except (IndexError, KeyError):
        return None


def fetch_gene_table(gene_id: str) -> str:
    """Fetch the NCBI Gene `gene_table` report (per-transcript exon tables)."""
    configure()
    pause()
    handle = Entrez.efetch(db="gene", id=gene_id, rettype="gene_table", retmode="text")
    text = handle.read()
    handle.close()
    return text


def parse_gene_table(text: str, target_mrna: str) -> dict | None:
    """Extract the coding-exon structure of one transcript from a `gene_table`.

    Columns are tab-separated:
        Genomic Exon | Genomic Coding | Gene Exon | Gene Coding | ExonLen | CodingLen | IntronLen
    A row with >= 6 fields is a coding exon whose coding length is field 5;
    shorter rows are UTR-only exons and contribute nothing to the CDS.

    Returns {'n_exons', 'cds_junctions', 'total_coding'} where `cds_junctions`
    are the CDS-coordinate offsets at which each non-first coding exon begins,
    or None if the transcript is absent or unparseable.
    """
    base = target_mrna.split(".")[0]

    for section in re.split(r"Exon table for\s+mRNA\s+", text)[1:]:
        lines = section.splitlines()
        if not lines or lines[0].split()[0].split(".")[0] != base:
            continue

        divider = next((i for i, l in enumerate(lines) if set(l.strip()) == {"-"}), None)
        if divider is None:
            return None

        coding_lengths = []
        for line in lines[divider + 1:]:
            if not line.strip():
                break
            fields = [f for f in re.split(r"\t+", line.strip()) if f]
            if len(fields) >= 6:
                try:
                    coding_lengths.append(int(fields[5]))
                except ValueError:
                    pass

        if not coding_lengths:
            return None
        return {"n_exons": len(coding_lengths),
                "cds_junctions": cumulative_junctions(coding_lengths),
                "total_coding": sum(coding_lengths)}
    return None


def cumulative_junctions(exon_lengths: list[int]) -> list[int]:
    """Cumulative offsets of every exon boundary except the final CDS end."""
    junctions, running = [], 0
    for length in exon_lengths[:-1]:
        running += length
        junctions.append(running)
    return junctions
