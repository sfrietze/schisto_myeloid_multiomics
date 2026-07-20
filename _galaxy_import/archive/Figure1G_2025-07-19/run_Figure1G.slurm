#!/bin/bash
#SBATCH --job-name=Figure1G
#SBATCH --output=Figure1G_%j.out
#SBATCH --error=Figure1G_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

cd /slipstream_old/home/sethandy/projects/keke

Rscript reproducibility/Figure1G/07_Figure1G_DAR_heatmap_lineplot.R
