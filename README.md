# Predicted CRISPR-Cas9 off-targets across plant genomes

This pipeline automatically retrieves data on the number of predicted CRISPR off-target
across plant genomes (47 species), with all gRNAs possibly designed on the 100 more conserved bp of 10 different genes.

Ten conserved *Arabidopsis thaliana* housekeeping genes are used as bait.
For each gene :

-the pipeline finds the ortholog in each species using BLAST on NCBI,

-data is manually assessed. When needed, data is manually added,

-the pipeline reconstructs coding sequences, selects the most conserved region in a exon for all species, and submits each genes sequence in that window to CRISPOR (https://crispor.gi.ucsc.edu) against each species' own genome assembly.

-the pipeline then retrieves the results such as predicted off-target sites for all guides and localisation of said sites (genic/intergenic).

Verification was conducted at each step. 

Code was written using Claude code, Opus model.

---

Python 3.12
Biopython 1.87
Requests 2.34.2
Openpyxl 3.1.5
CRISPOR 5.2

need do setup e-mail adress for NCBI :


export NCBI_EMAIL="you@example.org" or $env:NCBI_EMAIL on Windows

---

## The pipeline

| # | Command | Does | Network |
|---|---|---|---|
| 0 | `python fetch_at_cds.py` | Fetches the 10 *Arabidopsis* reference CDSs | NCBI |
| 1 | `python step1_blast.py --all` | Homolog search : Does one broad blastp per gene against Viridiplantae; retrives best hit per species | NCBI |
| 1b | `python step1_fallback.py --all --resume` |Homolog search : Retries species the broad search missed, one query each against the species needed | NCBI | 
| 1c | `python make_manual_template.py` | Writes the Excel form for species NCBI does not cover for manual addition of sequences (ex: *Rubus occidentalis* and *Miscanthus sinensis*). Found sequences can also be replaced if manually found to be wrong  | no | 
| 1d | `python build_manual_overrides.py` | Converts the filled form to config | no | 
| 2 | `python step2_fetch_genomic.py --all --resume` | Reconstructs each ortholog's spliced CDS from its `coded_by` intervals | NCBI | 
| 2b | `python step2b_exon_structure.py --all --resume` | Adds coding-exon junctions in CDS coordinates | NCBI | 
| 3-4 | `python step34_conserved.py --all` | Star alignment, each species CDS is aligned to the *Arabidopsis* reference (global, NUC.4.4, gap-open −10, gap-extend −0.5), then picks the most conserved 100bp single-exon window | no | 
| 5 | `python step5_loop.py` | Submits every sequence from that window to CRISPOR (SpCas9, NGG PAM) and collects results | CRISPOR | 
| 5b | `python step5b_recover_genomic.py --all --resume` | Recovers windows CRISPOR could not find in its assembly | NCBI + CRISPOR | 
| 6a | `python step6_aggregate_tables.py` | Builds the analysis CSVs | no |



---

## Configuration files


| File | Contents |
|---|---|
| `config/gene_proteins.json` | The 10 bait gene names and *Arabidopsis* protein sequences |
| `config/species_genomes.json` | 47 species -> CRISPOR genome id (34 NCBI `GCF_`/`GCA_` accessions, 13 Ensembl Plants / Phytozome / GenArk identifiers) |
| `config/species_metadata.json` | 47 species -> genome size (Mb) and ploidy |
| `config/at_cds.json` | The 10 reference CDSs with their exon coordinates |
| `config/manual_cds_overrides.json` | 26 hand-curated entries |


---

## Output data



- `results/step1`, `step2`, `step34` and `step5` are the raw results. Used for manual review during the process.
- `results/alignments/{gene}_alignment.fasta` — the star alignment in reference coordinates, with the selected window marked
- `results/alignments/{gene}_conservation.txt` — per-position identity across the selected window


