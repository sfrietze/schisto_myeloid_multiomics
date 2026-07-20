#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
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

res_male_vs_female_IF <- results(
  dds_combined,
  name = "sex_M_vs_F"
)

res_male_vs_female_IF <- res_male_vs_female_IF[
  order(res_male_vs_female_IF$padj),
]

sig_dars <- res_male_vs_female_IF[
  which(res_male_vs_female_IF$padj < 0.05),
]

gr_dars_sig <- rowRanges(dds_combined)[
  as.integer(rownames(sig_dars))
]

saveRDS(
  res_male_vs_female_IF,
  file.path(
    RESULTS_DIR,
    "atacseq",
    "figure1",
    "res_male_vs_female_IF.rds"
  )
)

saveRDS(
  gr_dars_sig,
  file.path(
    RESULTS_DIR,
    "atacseq",
    "figure1",
    "gr_dars_sig_interaction.rds"
  )
)
