#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(DiffBind)
})

source("config.R")

dbaObj_peaks <- readRDS(
  file.path(
    RESULTS_DIR,
    "atacseq",
    "figure1",
    "Figure1_ATAC_consensus_counts.rds"
  )
)

counts_mat <- dba.peakset(
  dbaObj_peaks,
  bRetrieve = TRUE,
  DataType = DBA_DATA_FRAME
)

col_data <- dbaObj_peaks$samples

counts_mat_clean <- round(counts_mat[, 4:ncol(counts_mat)])

rownames(col_data) <- col_data$SampleID
stopifnot(all(colnames(counts_mat_clean) == rownames(col_data)))

col_data$condition <- factor(col_data$Condition)
col_data$sex <- factor(
  ifelse(
    grepl("^IF|^UF|^SEA|^NOSEA", rownames(col_data)),
    "F",
    "M"
  )
)

gr_from_diffbind <- dba.peakset(dbaObj_peaks, bRetrieve = TRUE)

keep <- col_data$condition %in% c("IF", "UF")

col_data <- droplevels(col_data[keep, ])
counts_mat_clean <- counts_mat_clean[, keep]

dds_combined <- DESeqDataSetFromMatrix(
  countData = counts_mat_clean,
  colData = col_data,
  design = ~ sex + condition + sex:condition,
  rowRanges = gr_from_diffbind
)

dds_combined <- DESeq(dds_combined)

saveRDS(
  dds_combined,
  file.path(
    RESULTS_DIR,
    "atacseq",
    "figure1",
    "dds_combined.rds"
  )
)
