#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(DiffBind)
  library(GenomicRanges)
})

source("config.R")

dds_female_sea <- readRDS(file.path(DAR_DIR, "female_sea_dds.rds"))
dds_male_inf   <- readRDS(file.path(DAR_DIR, "male_inf_dds.rds"))
dds_female_inf <- readRDS(file.path(DAR_DIR, "female_inf_dds.rds"))

consensus_peaks <- reduce(c(
  rowRanges(dds_female_sea),
  rowRanges(dds_male_inf),
  rowRanges(dds_female_inf)
))

samples <- read.csv(
  file.path(DAR_DIR, "all_samples.csv"),
  stringsAsFactors = FALSE
)

samples$X <- NULL
samples$PeakCaller <- "bed"
samples$bamReads <- file.path(DAR_DIR, samples$bamReads)
samples$Peaks <- file.path(DAR_DIR, samples$Peaks)

dbObj <- dba(sampleSheet = samples)

dbaObj_peaks <- dba.count(
  dbObj,
  peaks = consensus_peaks,
  minOverlap = 1
)

dir.create(
  file.path(RESULTS_DIR, "atacseq", "figure1"),
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  dbaObj_peaks,
  file.path(
    RESULTS_DIR,
    "atacseq",
    "figure1",
    "Figure1_ATAC_consensus_counts.rds"
  )
)
