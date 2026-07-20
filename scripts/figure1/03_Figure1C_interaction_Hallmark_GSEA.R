#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(fgsea)
  library(msigdbr)
})

set.seed(123)

dds_file <- "results/rnaseq/figure1/Figure1A_dds.rds"
figure_dir <- "figures/main/Figure1"
results_dir <- "results/rnaseq/figure1"

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

dds <- readRDS(dds_file)

# Interaction coefficient:
# positive statistic = stronger infection response in males
# negative statistic = stronger infection response in females
interaction_name <- "sexM.conditionI"

if (!interaction_name %in% resultsNames(dds)) {
  stop(
    "Coefficient not found: ", interaction_name,
    "\nAvailable coefficients: ",
    paste(resultsNames(dds), collapse = ", ")
  )
}

res_interaction <- results(
  dds,
  name = interaction_name
)

interaction_df <- as.data.frame(res_interaction) %>%
  rownames_to_column("gene") %>%
  filter(
    !is.na(stat),
    gene != ""
  ) %>%
  group_by(gene) %>%
  slice_max(
    order_by = abs(stat),
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

# Rank genes using the DESeq2 Wald statistic
gene_list <- interaction_df$stat
names(gene_list) <- interaction_df$gene
gene_list <- sort(gene_list, decreasing = TRUE)

# Hallmark gene sets
msig_hallmark <- msigdbr(
  species = "Mus musculus",
  category = "H"
)

hallmark_list <- split(
  x = msig_hallmark$gene_symbol,
  f = msig_hallmark$gs_name
)

# Run GSEA
fgsea_res <- fgsea(
  pathways = hallmark_list,
  stats = gene_list,
  minSize = 15,
  maxSize = 500

) %>%
  arrange(padj)

# Select the 10 strongest significant pathways in each direction
plot_df <- fgsea_res %>%
  filter(
    !is.na(padj),
    padj < 0.05
  ) %>%
  mutate(
    neglog10_padj = -log10(padj),
    absNES = abs(NES),
    direction = ifelse(
      NES > 0,
      "Male-biased",
      "Female-biased"
    )
  ) %>%
  group_by(direction) %>%
  slice_max(
    order_by = absNES,
    n = 10,
    with_ties = FALSE
  ) %>%
  ungroup()

if (nrow(plot_df) == 0) {
  stop("No Hallmark pathways passed padj < 0.05.")
}

# Clean pathway labels
plot_df <- plot_df %>%
  mutate(
    pathway = gsub("^HALLMARK_", "", pathway),
    pathway = gsub("_", " ", pathway)
  ) %>%
  arrange(NES)

plot_df$pathway <- factor(
  plot_df$pathway,
  levels = plot_df$pathway
)

message("Female-biased pathways: ",
        sum(plot_df$direction == "Female-biased"))
message("Male-biased pathways: ",
        sum(plot_df$direction == "Male-biased"))

# Save complete and plotted GSEA results
fgsea_export <- as.data.frame(fgsea_res) %>%
  mutate(
    leadingEdge = vapply(
      leadingEdge,
      function(x) paste(x, collapse = ";"),
      character(1)
    )
  )

plot_export <- as.data.frame(plot_df) %>%
  mutate(
    leadingEdge = vapply(
      leadingEdge,
      function(x) paste(x, collapse = ";"),
      character(1)
    )
  )

write.csv(
  fgsea_export,
  file.path(
    results_dir,
    "Figure1C_Hallmark_GSEA_all_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  plot_export,
  file.path(
    results_dir,
    "Figure1C_Hallmark_GSEA_plotted_pathways.csv"
  ),
  row.names = FALSE
)

# Figure settings matching the manuscript panel
p_gsea <- ggplot(
  plot_df,
  aes(x = NES, y = pathway)
) +
  geom_vline(
    xintercept = seq(-3, 3, by = 1),
    color = "gray85",
    linewidth = 0.5
  ) +
  geom_hline(
    yintercept = seq_along(levels(plot_df$pathway)),
    color = "gray90",
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = 0,
    color = "black",
    linewidth = 0.8
  ) +
  geom_point(
    aes(
      size = absNES,
      color = neglog10_padj
    ),
    alpha = 0.95
  ) +
  scale_color_gradient(
    low = "dodgerblue3",
    high = "red2",
    name = expression(-log[10]~adj.~p-value)
  ) +
  scale_size(
    range = c(3, 7),
    name = "|NES|"
  ) +
  scale_x_continuous(
    breaks = c(-2, 0, 2),
    limits = c(-3.3, 3.3),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL
  ) +
  annotate(
    "text",
    x = -1.65,
    y = nrow(plot_df) + 1.15,
    label = "Female-biased",
    fontface = "bold",
    size = 4.2
  ) +
  annotate(
    "segment",
    x = -0.25,
    xend = -1.05,
    y = nrow(plot_df) + 1.15,
    yend = nrow(plot_df) + 1.15,
    arrow = arrow(
      length = unit(0.12, "inches"),
      type = "closed"
    ),
    linewidth = 0.7
  ) +
  annotate(
    "text",
    x = 1.65,
    y = nrow(plot_df) + 1.15,
    label = "Male-biased",
    fontface = "bold",
    size = 4.2
  ) +
  annotate(
    "segment",
    x = 0.25,
    xend = 1.05,
    y = nrow(plot_df) + 1.15,
    yend = nrow(plot_df) + 1.15,
    arrow = arrow(
      length = unit(0.12, "inches"),
      type = "closed"
    ),
    linewidth = 0.7
  ) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    axis.line = element_line(
      color = "black",
      linewidth = 0.6
    ),
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.6
    ),
    axis.text.y = element_text(
      size = 10,
      color = "gray35"
    ),
    axis.text.x = element_text(
      size = 11,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 13,
      color = "black"
    ),
    legend.position = "right",
    plot.margin = margin(
      t = 35,
      r = 10,
      b = 10,
      l = 10
    )
  )

ggsave(
  file.path(
    figure_dir,
    "Figure1C_interaction_Hallmark_GSEA.pdf"
  ),
  p_gsea,
  width = 6,
  height = 6,
  useDingbats = FALSE
)

ggsave(
  file.path(
    figure_dir,
    "Figure1C_interaction_Hallmark_GSEA.png"
  ),
  p_gsea,
  width = 6,
  height = 6,
  units = "in",
  dpi = 600
)
