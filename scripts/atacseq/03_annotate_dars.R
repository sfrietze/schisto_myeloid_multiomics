#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(GenomicRanges)
})

source("config.R")

dds_combined <- readRDS(
  file.path(
    RESULTS_DIR,
    "atacseq",
    "figure1",
    "dds_combined.rds"
  )
)

dds_combined$sex <- relevel(dds_combined$sex, ref = "F")
dds_combined$condition <- relevel(dds_combined$condition, ref = "UF")

res_sex_main <- results(dds_combined, name = "sex_M_vs_F")
res_interact <- results(dds_combined, name = "sexM.conditionUF")

gr <- rowRanges(dds_combined)

mcols(gr)$log2FC_IF <- res_sex_main$log2FoldChange
mcols(gr)$padj_IF <- res_sex_main$padj
mcols(gr)$group_IF <- ifelse(
  res_sex_main$padj < 0.05 & !is.na(res_sex_main$padj),
  ifelse(
    res_sex_main$log2FoldChange > 1,
    "Male > Female (IF)",
    ifelse(
      res_sex_main$log2FoldChange < -1,
      "Female > Male (IF)",
      NA
    )
  ),
  NA
)

mcols(gr)$log2FC_UF <- res_interact$log2FoldChange
mcols(gr)$padj_UF <- res_interact$padj
mcols(gr)$group_UF <- ifelse(
  res_interact$padj < 0.05 & !is.na(res_interact$padj),
  ifelse(
    res_interact$log2FoldChange > 1,
    "Male > Female (UF)",
    ifelse(
      res_interact$log2FoldChange < -1,
      "Female > Male (UF)",
      NA
    )
  ),
  NA
)

mcols(gr)$source <- apply(
  cbind(mcols(gr)$group_IF, mcols(gr)$group_UF),
  1,
  function(x) paste(na.omit(unique(x)), collapse = ";")
)

gr_annotated <- gr[mcols(gr)$source != "", ]

saveRDS(
  gr_annotated,
  file.path(
    RESULTS_DIR,
    "atacseq",
    "figure1",
    "fig1_DARs_combined_annotated.rds"
  )
)
