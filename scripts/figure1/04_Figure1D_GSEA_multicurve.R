#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(fgsea)
  library(msigdbr)
  library(patchwork)
})

set.seed(123)

figure_dir <- file.path("figures", "main", "Figure1")
results_dir <- file.path("results", "rnaseq", "figure1")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

dds_candidates <- c(
  file.path(results_dir, "Figure1A_dds.rds"),
  file.path("results", "rnaseq", "Figure1A_dds.rds"),
  file.path("results", "Figure1A_dds.rds")
)

dds_path <- dds_candidates[file.exists(dds_candidates)][1]

if (is.na(dds_path)) {
  stop(
    "Could not locate Figure1A_dds.rds. Checked:\n",
    paste(dds_candidates, collapse = "\n")
  )
}

message("Loading: ", dds_path)

dds <- readRDS(dds_path)

interaction_coef <- "sexM.conditionI"

if (!interaction_coef %in% resultsNames(dds)) {
  stop(
    "Coefficient ", interaction_coef, " not found. Available coefficients:\n",
    paste(resultsNames(dds), collapse = "\n")
  )
}

res_interaction <- results(
  dds,
  name = interaction_coef
)

interaction_df <- as.data.frame(res_interaction) %>%
  rownames_to_column("gene") %>%
  filter(
    !is.na(gene),
    gene != "",
    !is.na(log2FoldChange),
    is.finite(log2FoldChange)
  )

write.csv(
  interaction_df,
  file.path(results_dir, "Figure1D_interaction_DESeq2_results.csv"),
  row.names = FALSE
)

ranked_df <- interaction_df %>%
  arrange(desc(abs(log2FoldChange))) %>%
  distinct(gene, .keep_all = TRUE) %>%
  arrange(desc(log2FoldChange))

gene_ranks <- ranked_df$log2FoldChange
names(gene_ranks) <- ranked_df$gene
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

message("Ranked genes: ", length(gene_ranks))

hallmark_df <- msigdbr(
  species = "Mus musculus",
  category = "H"
)

hallmark_list <- split(
  hallmark_df$gene_symbol,
  hallmark_df$gs_name
)

selected_pathways <- c(
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE"
)

missing_pathways <- setdiff(
  selected_pathways,
  names(hallmark_list)
)

if (length(missing_pathways) > 0) {
  stop(
    "Missing Hallmark pathways:\n",
    paste(missing_pathways, collapse = "\n")
  )
}

selected_sets <- hallmark_list[selected_pathways]

fgsea_res <- fgseaMultilevel(
  pathways = selected_sets,
  stats = gene_ranks,
  minSize = 15,
  maxSize = 500,
  eps = 0
)

fgsea_export <- as.data.frame(fgsea_res) %>%
  mutate(
    leadingEdge = vapply(
      leadingEdge,
      function(x) paste(x, collapse = ";"),
      character(1)
    )
  )

write.csv(
  fgsea_export,
  file.path(results_dir, "Figure1D_selected_Hallmark_GSEA_results.csv"),
  row.names = FALSE
)

pathway_labels <- c(
  HALLMARK_CHOLESTEROL_HOMEOSTASIS = "Cholesterol homeostasis",
  HALLMARK_FATTY_ACID_METABOLISM = "Fatty acid metabolism",
  HALLMARK_E2F_TARGETS = "E2F targets",
  HALLMARK_INTERFERON_ALPHA_RESPONSE = "Interferon alpha response"
)

pathway_colors <- c(
  HALLMARK_CHOLESTEROL_HOMEOSTASIS = "#E69F00",
  HALLMARK_FATTY_ACID_METABOLISM = "#009E73",
  HALLMARK_E2F_TARGETS = "#0072B2",
  HALLMARK_INTERFERON_ALPHA_RESPONSE = "#CC79A7"
)

extract_curve <- function(pathway_name) {
  enrichment_plot <- plotEnrichment(
    selected_sets[[pathway_name]],
    gene_ranks
  )

  curve_df <- ggplot_build(enrichment_plot)$data[[1]]

  tibble(
    pathway = pathway_name,
    rank = curve_df$x,
    enrichment_score = curve_df$y
  )
}

curve_df <- bind_rows(
  lapply(selected_pathways, extract_curve)
)

tick_df <- bind_rows(
  lapply(
    selected_pathways,
    function(pathway_name) {
      genes <- intersect(
        selected_sets[[pathway_name]],
        names(gene_ranks)
      )

      tibble(
        pathway = pathway_name,
        rank = match(genes, names(gene_ranks))
      )
    }
  )
)

curve_df$pathway <- factor(
  curve_df$pathway,
  levels = selected_pathways
)

tick_df$pathway <- factor(
  tick_df$pathway,
  levels = rev(selected_pathways)
)

fgsea_summary <- as.data.frame(fgsea_res) %>%
  select(pathway, NES, padj) %>%
  mutate(
    pathway_label = unname(pathway_labels[pathway])
  )

curve_annotations <- curve_df %>%
  group_by(pathway) %>%
  slice_max(
    order_by = abs(enrichment_score),
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  left_join(
    fgsea_summary,
    by = "pathway"
  ) %>%
  mutate(
    label = paste0(
      pathway_label,
      "\nNES = ",
      sprintf("%.2f", NES),
      ", FDR = ",
      format.pval(padj, digits = 2, eps = 0.001)
    )
  )

curve_plot <- ggplot(
  curve_df,
  aes(
    x = rank,
    y = enrichment_score,
    color = pathway
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.4,
    color = "grey60"
  ) +
  geom_line(
    linewidth = 1
  ) +
  geom_text_repel(
    data = curve_annotations,
    aes(
      x = rank,
      y = enrichment_score,
      label = label,
      color = pathway
    ),
    size = 3.2,
    box.padding = 0.5,
    point.padding = 0.25,
    min.segment.length = 0,
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = pathway_colors,
    labels = pathway_labels
  ) +
  scale_x_continuous(
    limits = c(1, length(gene_ranks)),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Enrichment score"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(5, 10, 0, 10)
  )

tick_plot <- ggplot(
  tick_df,
  aes(
    x = rank,
    y = pathway,
    color = pathway
  )
) +
  geom_point(
    shape = 124,
    size = 2.7,
    stroke = 0.7
  ) +
  scale_color_manual(
    values = pathway_colors
  ) +
  scale_x_continuous(
    limits = c(1, length(gene_ranks)),
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    labels = function(x) unname(pathway_labels[x])
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(0, 10, 0, 10)
  )

rank_metric_df <- tibble(
  gene = names(gene_ranks),
  rank = seq_along(gene_ranks),
  metric = as.numeric(gene_ranks)
)

label_genes <- c(
  "Isg15",
  "C3",
  "Cfb",
  "Stat1",
  "Top2a",
  "Irf7",
  "C4b",
  "C1qb",
  "Mcm5",
  "Ifih1",
  "Pcna",
  "Hk3",
  "Fabp4",
  "Mgll",
  "Dio2"
)

gene_label_df <- rank_metric_df %>%
  filter(gene %in% label_genes)

rank_plot <- ggplot(
  rank_metric_df,
  aes(
    x = rank,
    y = metric
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.4,
    color = "grey60"
  ) +
  geom_area(
    aes(
      fill = metric > 0
    ),
    alpha = 0.75,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c(
      `TRUE` = "firebrick",
      `FALSE` = "royalblue"
    )
  ) +
  geom_text_repel(
    data = gene_label_df,
    aes(
      label = gene
    ),
    size = 3,
    box.padding = 0.35,
    point.padding = 0.15,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 123
  ) +
  annotate(
    "text",
    x = length(gene_ranks) * 0.03,
    y = max(rank_metric_df$metric, na.rm = TRUE) * 0.9,
    label = "Male-biased",
    hjust = 0,
    fontface = "bold",
    size = 3.5
  ) +
  annotate(
    "text",
    x = length(gene_ranks) * 0.97,
    y = min(rank_metric_df$metric, na.rm = TRUE) * 0.9,
    label = "Female-biased",
    hjust = 1,
    fontface = "bold",
    size = 3.5
  ) +
  scale_x_continuous(
    limits = c(1, length(gene_ranks)),
    expand = c(0, 0)
  ) +
  labs(
    x = "Rank in ordered gene list",
    y = "Interaction log2 fold change"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.margin = margin(0, 10, 5, 10)
  )

figure_1d <- curve_plot /
  tick_plot /
  rank_plot +
  plot_layout(
    heights = c(3.5, 1.4, 2.2)
  )

ggsave(
  filename = file.path(
    figure_dir,
    "Figure1D_GSEA_multicurve_labeled.pdf"
  ),
  plot = figure_1d,
  width = 9,
  height = 8,
  units = "in",
  device = grDevices::cairo_pdf
)

ggsave(
  filename = file.path(
    figure_dir,
    "Figure1D_GSEA_multicurve_labeled.png"
  ),
  plot = figure_1d,
  width = 9,
  height = 8,
  units = "in",
  dpi = 300
)

saveRDS(
  list(
    gene_ranks = gene_ranks,
    selected_sets = selected_sets,
    fgsea_results = fgsea_res,
    curve_data = curve_df,
    tick_data = tick_df
  ),
  file.path(
    results_dir,
    "Figure1D_GSEA_plot_data.rds"
  )
)

message("Figure 1D complete")
message("PDF: ", file.path(
  figure_dir,
  "Figure1D_GSEA_multicurve_labeled.pdf"
))
message("PNG: ", file.path(
  figure_dir,
  "Figure1D_GSEA_multicurve_labeled.png"
))
