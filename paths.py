"""Filesystem layout, config loading and JSON helpers.

All outputs live under `results/`, one directory per step, so that every stage
can be re-run independently against the previous stage's files.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

# ── Layout ────────────────────────────────────────────────────────────────────

ROOT = Path(__file__).resolve().parent          # repository root

CONFIG_DIR   = ROOT / "config"
RESULTS_DIR  = ROOT / "results"

STEP1_DIR    = RESULTS_DIR / "step1"            # ortholog BLAST hits
STEP2_DIR    = RESULTS_DIR / "step2"            # reconstructed CDS + exon structure
STEP34_DIR   = RESULTS_DIR / "step34"           # conserved windows
STEP5_DIR    = RESULTS_DIR / "step5"            # CRISPOR off-targets
ALIGN_DIR    = RESULTS_DIR / "alignments"       # per-gene FASTA + conservation tables
R_EXPORT_DIR = RESULTS_DIR / "R_export"         # tidy CSVs + R analysis scripts

# The committed config: the gene and species panels the published results were
# produced from. at_cds.json is written by fetch_at_cds.py and
# manual_cds_overrides.json by build_manual_overrides.py.
GENE_PROTEINS_FILE    = CONFIG_DIR / "gene_proteins.json"
SPECIES_GENOMES_FILE  = CONFIG_DIR / "species_genomes.json"
SPECIES_METADATA_FILE = CONFIG_DIR / "species_metadata.json"
AT_CDS_FILE           = CONFIG_DIR / "at_cds.json"
MANUAL_OVERRIDES_FILE = CONFIG_DIR / "manual_cds_overrides.json"

# Tandem-repeat / multi-copy families: a guide in these genes also matches its own
# near-identical paralogs, so their off-target counts sit above the single-copy
# trend. Flagged in every output table; no analysis filters on them.
MULTICOPY_GENES = frozenset({"AtUBQ2", "AtHistone4"})

# Length of the conserved window submitted to CRISPOR, in bp.
WINDOW_BP = 100


# ── Helpers ───────────────────────────────────────────────────────────────────

def gene_slug(gene_name: str) -> str:
    """Filename-safe form of a gene name: 'AtGTP_EFTU' -> 'atgtp_eftu'."""
    return re.sub(r"[^a-z0-9]", "_", gene_name.lower()).strip("_")


def load_json(path: Path) -> dict:
    """Read a JSON file, returning {} when it does not exist."""
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def save_json(path: Path, data) -> None:
    """Write JSON with a stable indent. Steps call this after every species so an
    interrupted run never loses more than the one in-flight query."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")


def load_genes() -> dict[str, str]:
    """{gene name: Arabidopsis bait protein sequence}, in panel order."""
    genes = load_json(GENE_PROTEINS_FILE)
    if not genes:
        raise SystemExit(f"{GENE_PROTEINS_FILE} is missing or empty")
    return genes


def gene_names() -> list[str]:
    return list(load_genes().keys())


def load_species_genomes() -> dict[str, str]:
    """{species: CRISPOR genome identifier}."""
    genomes = load_json(SPECIES_GENOMES_FILE)
    if not genomes:
        raise SystemExit(f"{SPECIES_GENOMES_FILE} is missing or empty.")
    return genomes


def load_species_metadata() -> dict[str, dict]:
    """{species: {'genome_size_mb': float, 'ploidy': str}}."""
    return load_json(SPECIES_METADATA_FILE)


def resolve_genes(gene: str | None) -> list[str]:
    """Turn a `--gene X` argument (or None for all) into a validated gene list."""
    available = gene_names()
    if gene is None:
        return available
    if gene not in available:
        raise SystemExit(f"Unknown gene {gene!r}. Available: {', '.join(available)}")
    return [gene]


def ploidy_group(ploidy: str) -> str:
    """Collapse the four ploidy levels to the diploid / polyploid contrast."""
    return "diploid" if str(ploidy).strip().lower() == "diploid" else "polyploid"
