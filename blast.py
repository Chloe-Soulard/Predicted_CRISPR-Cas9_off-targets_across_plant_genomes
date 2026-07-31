"""BLAST hit parsing and species-name matching.

Shared by the broad ortholog search (`step1_blast.py`) and the per-species
fallback (`step1_fallback.py`) so both record hits in the same shape.
"""

from __future__ import annotations

import re

NCBI_PROTEIN_URL = "https://www.ncbi.nlm.nih.gov/protein/{}"


def extract_organism(title: str) -> str:
    """Pull the organism out of a BLAST title: 'actin [Zea mays]' -> 'Zea mays'."""
    m = re.search(r"\[([^\]]+)\]", title)
    return m.group(1).strip() if m else ""


def normalize_species(name: str) -> str:
    """Reduce a name to lowercase genus + species, dropping subspecies/cultivar."""
    parts = name.lower().split()
    return " ".join(parts[:2]) if len(parts) >= 2 else name.lower()


def match_species(blast_organism: str, target_species: list[str]) -> str | None:
    """Map a BLAST hit organism onto one of the panel species, or None.

    Genus + species comparison, so a hit from 'Brassica rapa subsp. pekinensis'
    matches the panel entry 'Brassica rapa'.
    """
    blast_norm = normalize_species(blast_organism)
    for species in target_species:
        if blast_norm.startswith(normalize_species(species)):
            return species
    return None


def hit_record(alignment) -> dict:
    """Summarise one BLAST alignment (its best HSP) as a JSON-serialisable dict."""
    hsp = alignment.hsps[0]
    return {
        "accession":    alignment.accession,
        "title":        alignment.title[:140],
        "organism":     extract_organism(alignment.title),
        "identity_pct": round(100.0 * hsp.identities / hsp.align_length, 2),
        "evalue":       hsp.expect,
        "bits":         hsp.bits,
        "length":       alignment.length,
        "ncbi_url":     NCBI_PROTEIN_URL.format(alignment.accession),
    }
