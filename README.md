# Schistosoma myeloid multiomics

This repository contains reproducible analyses for the ATAC-seq and RNA-seq datasets used in the manuscript.

## Repository structure

- `config.R` — global paths used by all scripts
- `data/` — input datasets
- `results/` — intermediate analysis results
- `figures/` — final figures
- `scripts/` — analysis scripts
- `metadata/` — sample metadata
- `docs/` — additional documentation

## ATAC-seq

Input data are stored in:

- `data/atacseq/bigwig/`
- `data/atacseq/combined/`
- `data/atacseq/peaks/`
- `data/atacseq/processed/`
- `data/atacseq/tobias/`

ATAC-seq scripts are stored in:

- `scripts/atacseq/`

ATAC-seq outputs are written to:

- `results/atacseq/`
- `figures/atacseq/`

## RNA-seq

RNA-seq scripts are stored in:

- `scripts/rnaseq/`

RNA-seq outputs are written to:

- `results/rnaseq/`

All active scripts use repository-relative paths defined in `config.R`.
