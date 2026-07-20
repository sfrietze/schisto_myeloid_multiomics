#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

dds_file <- file.path(
  "results",
  "atacseq",
  "figure1",
  "Figure1F_ATAC_interaction_dds.rds"
)

results_dir <- file.path(
  "results",
  "atacseq",
  "figure1"
)

figure_dir <- file.path(
  "figures",
  "main",
  "Figure1"
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(dds_file)) {
  stop(
    "DESeq2 object not found. Run ",
    "scripts/atacseq/00_build_ATAC_interaction_dds.R first."
  )
}

dds <- readRDS(dds_file)

coefficient <- "sex_M_vs_F"

if (!coefficient %in% resultsNames(dds)) {
  stop(
    "Coefficient not found: ",
    coefficient,
    "\nAvailable coefficients:\n",
    paste(resultsNames(dds), collapse = "\n")
  )
}

res <- results(
  dds,
  name = coefficient,
  alpha = 0.05
)

res_df <- as.data.frame(res) %>%
  rownames_to_column("peak_id") %>%
  mutate(
    padj_plot = ifelse(
      is.na(padj),
      1,
      padj
    ),
    significance = case_when(
      padj_plot < 0.05 &
        log2FoldChange < -1 ~ "Up in Female",

      padj_plot < 0.05 &
        log2FoldChange > 1 ~ "Up in Male",

      TRUE ~ "NS"
    )
  ) %>%
  filter(
    !is.na(log2FoldChange),
    is.finite(log2FoldChange)
  )

female_count <- sum(
  res_df$significance == "Up in Female"
)

male_count <- sum(
  res_df$significance == "Up in Male"
)

message("Up in Female: ", female_count)
message("Up in Male: ", male_count)

write.csv(
  res_df,
  file.path(
    results_dir,
    "Figure1F_DAR_results.csv"
  ),
  row.names = FALSE
)

plot_df <- res_df %>%
  mutate(
    significance = factor(
      significance,
      levels = c(
        "NS",
        "Up in Female",
        "Up in Male"
      )
    )
  ) %>%
  arrange(significance)

legend_labels <- c(
  NS = "NS",
  `Up in Female` = paste0(
    "Up in Female (",
    female_count,
    ")"
  ),
  `Up in Male` = paste0(
    "Up in Male (",
    male_count,
    ")"
  )
)

figure_1f <- ggplot(
  plot_df,
  aes(
    x = log2FoldChange,
    y = -log10(padj_plot),
    color = significance
  )
) +
  geom_point(
    alpha = 0.8,
    size = 1.6
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dotted",
    linewidth = 0.6
  ) +
  scale_color_manual(
    values = c(
      NS = "grey75",
      `Up in Female` = "#2F5DA8",
      `Up in Male` = "#F02E3A"
    ),
    breaks = c(
      "Up in Female",
      "Up in Male"
    ),
    labels = c(
      paste0("Up in Female (", female_count, ")"),
      paste0("Up in Male (", male_count, ")")
    ),
    name = "DARs"
  ) +
  coord_cartesian(
    xlim = c(-10, 10)
  ) +
  labs(
    x = "log2 Fold Change (Male vs Female)",
    y = expression(-log[10] * " adjusted p-value")
  ) +
  theme_classic(
    base_size = 14
  ) +
  theme(
    legend.position = c(0.79, 0.91),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_text(
      size = 11
    ),
    legend.text = element_text(
      size = 10
    ),
    axis.title = element_text(
      size = 14
    ),
    axis.text = element_text(
      size = 11,
      color = "black"
    ),
    axis.line = element_line(
      linewidth = 0.8,
      color = "black"
    ),
    axis.ticks = element_line(
      linewidth = 0.8,
      color = "black"
    ),
    plot.margin = margin(
      10,
      10,
      10,
      10
    )
  )

ggsave(
  filename = file.path(
    figure_dir,
    "Figure1F_DAR_volcano.pdf"
  ),
  plot = figure_1f,
  width = 6.5,
  height = 6.5,
  units = "in",
  device = grDevices::cairo_pdf
)

ggsave(
  filename = file.path(
    figure_dir,
    "Figure1F_DAR_volcano.png"
  ),
  plot = figure_1f,
  width = 6.5,
  height = 6.5,
  units = "in",
  dpi = 300
)

message(
  "Saved: ",
  file.path(
    figure_dir,
    "Figure1F_DAR_volcano.pdf"
  )
)
