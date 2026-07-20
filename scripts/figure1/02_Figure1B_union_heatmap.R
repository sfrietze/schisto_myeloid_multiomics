#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(DESeq2)
})

padj_cutoff <- 0.05
lfc_cutoff <- 1

dds_file <- "results/rnaseq/figure1/Figure1A_dds.rds"
figure_dir <- "figures/main/Figure1"
results_dir <- "results/rnaseq/figure1"

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

dds <- readRDS(dds_file)

## Interaction results from the full model
res_interaction <- results(
  dds,
  name = "sexM.conditionI"
)

## Male infection comparison
dds_m <- dds[, dds$sex == "M"]
dds_m$condition <- droplevels(dds_m$condition)
design(dds_m) <- ~ condition
dds_m <- DESeq(dds_m)

## Match original U versus I direction explicitly
res_inf_m <- results(
  dds_m,
  contrast = c("condition", "U", "I")
)

## Female infection comparison
dds_f <- dds[, dds$sex == "F"]
dds_f$condition <- droplevels(dds_f$condition)
design(dds_f) <- ~ condition
dds_f <- DESeq(dds_f)

## Match original U versus I direction explicitly
res_inf_f <- results(
  dds_f,
  contrast = c("condition", "U", "I")
)

## Convert results to data frames
res_m <- as.data.frame(res_inf_m) %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::rename(
    log2FC_male = log2FoldChange,
    padj_male = padj
  )

res_f <- as.data.frame(res_inf_f) %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::rename(
    log2FC_female = log2FoldChange,
    padj_female = padj
  )

res_int <- as.data.frame(res_interaction) %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::rename(
    log2FC_int = log2FoldChange,
    padj_int = padj
  )

## Merge exactly as in the original analysis
merged <- res_m %>%
  left_join(res_f, by = "gene") %>%
  left_join(res_int, by = "gene")

## Identify DEGs
degs_male <- merged %>%
  filter(
    padj_male < padj_cutoff,
    abs(log2FC_male) > lfc_cutoff
  ) %>%
  pull(gene)

degs_female <- merged %>%
  filter(
    padj_female < padj_cutoff,
    abs(log2FC_female) > lfc_cutoff
  ) %>%
  pull(gene)

degs_interaction <- merged %>%
  filter(
    padj_int < padj_cutoff,
    abs(log2FC_int) > lfc_cutoff
  ) %>%
  pull(gene)

## Original union construction
deg_union <- union(
  union(degs_male, degs_female),
  degs_interaction
)

message("Male DEGs: ", length(degs_male))
message("Female DEGs: ", length(degs_female))
message("Interaction DEGs: ", length(degs_interaction))
message("Union DEGs: ", length(deg_union))

## rlog transformation of the complete dataset
rld <- rlog(dds, blind = FALSE)

## Original heatmap matrix
rld_mat <- assay(rld)[deg_union, , drop = FALSE]

## Column annotation
annotation_col <- as.data.frame(colData(dds)) %>%
  dplyr::select(sex, condition)

annotation_col$sex <- factor(
  annotation_col$sex,
  levels = c("M", "F")
)

annotation_col$condition <- factor(
  annotation_col$condition,
  levels = c("U", "I")
)

rownames(annotation_col) <- colnames(rld_mat)

## Original annotation colors
ann_colors <- list(
  sex = c("M" = "royalblue", "F" = "orange"),
  condition = c("U" = "gray60", "I" = "black")
)

## Save analysis objects
write.csv(
  merged,
  file.path(results_dir, "Figure1B_all_DESeq2_results.csv"),
  row.names = FALSE
)

writeLines(
  deg_union,
  file.path(results_dir, "Figure1B_DEG_union.txt")
)

saveRDS(
  rld_mat,
  file.path(results_dir, "Figure1B_rlog_union_matrix.rds")
)

## Exact heatmap analysis: row scaling and default hierarchical clustering
pdf(
  file.path(figure_dir, "Figure1B_DEG_union_heatmap.pdf"),
  width = 6,
  height = 6
)

pheatmap(
  rld_mat,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  scale = "row",
  color = colorRampPalette(
    c("royalblue", "white", "firebrick")
  )(100),
  show_rownames = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  fontsize_row = 7,
  fontsize_col = 9
)

dev.off()

png(
  file.path(figure_dir, "Figure1B_DEG_union_heatmap.png"),
  width = 6,
  height = 6,
  units = "in",
  res = 600
)

pheatmap(
  rld_mat,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  scale = "row",
  color = colorRampPalette(
    c("royalblue", "white", "firebrick")
  )(100),
  show_rownames = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  fontsize_row = 7,
  fontsize_col = 9
)

dev.off()
