# Honours 2026 gut microbiome processing

Snakemake workflow (Pawsey Setonix, Singularity containers) and R analysis for a longitudinal
metagenomic study of the gut microbiome (3 participants, 27 stool samples, Nanopore).

## Contents
- `Snakefile` - the workflow (read QC, host removal, Kraken2/Bracken, metaFlye, VIBRANT + geNomad, CheckV, CD-HIT, CoverM, Phold, PhaBOX2)
- `config.yaml` - parameters and database paths (reconstructed; see the CHECK notes inside)
- `sample_sheet.csv` - the 27 sample IDs
- `scripts/` - helper scripts called by the Snakefile (union of viral calls, CheckV filtering, pooling, confidence tiers, abundance table)
- `slurm_profile/` - Snakemake SLURM profile
- `Analysis.R` - downstream statistics and figures (R 4.6.1)

## Running
    snakemake --profile slurm_profile --use-singularity --cores all -n   # dry run
    snakemake --profile slurm_profile --use-singularity --cores all

Tool versions are the container image tags listed in the Snakefile and in the thesis Methods (Table 2).
