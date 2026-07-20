#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(ggplot2)
  library(data.table)
  library(patchwork)
})

ct2_file <- "Figure1G_ct2.rds"
out_prefix <- "Figure1G_DAR_heatmap_lineplot"

ct2 <- readRDS(ct2_file)

signal <- SummarizedExperiment::assay(ct2, "max")
row_meta <- as.data.frame(SummarizedExperiment::rowData(ct2))
col_meta <- as.data.frame(SummarizedExperiment::colData(ct2))

stopifnot(nrow(signal) == nrow(row_meta))
stopifnot(ncol(signal) == nrow(col_meta))

sample_names <- rownames(col_meta)
colnames(signal) <- sample_names

group_order <- c(
  "Female IF vs UF",
  "Male vs Female (UF)",
  "Male vs Female (IF)",
  "Male IF vs UF"
)

sample_order <- c(
  "UF_fem_inf",
  "IF_fem_inf",
  "UF_male_inf",
  "IF_male_inf"
)

keep <- !is.na(row_meta$group) &
  row_meta$group %in% group_order

signal <- signal[keep, sample_order, drop = FALSE]
row_meta <- row_meta[keep, , drop = FALSE]

row_meta$group <- factor(
  row_meta$group,
  levels = group_order
)

## Preserve the ordering already produced by sortRegions().
group_split <- split(
  seq_len(nrow(signal)),
  row_meta$group,
  drop = TRUE
)

group_counts <- vapply(group_split, length, integer(1))

cat("Groups used for plotting:\n")
print(group_counts)

## Build heatmap data.
heatmap_dt <- rbindlist(
  lapply(names(group_split), function(grp) {
    idx <- group_split[[grp]]
    mat <- signal[idx, , drop = FALSE]

    dt <- as.data.table(mat)
    dt[, region_order := seq_len(.N)]
    dt[, group := grp]

    melt(
      dt,
      id.vars = c("region_order", "group"),
      variable.name = "sample",
      value.name = "signal"
    )
  })
)

heatmap_dt[, group := factor(group, levels = group_order)]
heatmap_dt[, sample := factor(sample, levels = sample_order)]

## Separate row index within each group.
heatmap_dt[, plot_row := region_order, by = group]

heatmap_plot <- ggplot(
  heatmap_dt,
  aes(
    x = sample,
    y = plot_row,
    fill = signal
  )
) +
  geom_raster() +
  facet_grid(
    group ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_fill_gradientn(
    colours = c(
      "#577AB2",
      "#FFFFE0",
      "#BC412B"
    ),
    limits = c(0, 150),
    oob = scales::squish,
    name = "Signal"
  ) +
  scale_y_reverse(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    strip.placement = "outside",
    strip.text.y.left = element_text(
      angle = 0,
      face = "bold",
      size = 10
    ),
    strip.background = element_blank(),
    panel.spacing.y = grid::unit(0.08, "lines"),
    plot.margin = margin(5, 5, 5, 5)
  )

## The saved ct2 assay contains one summarized value per region and sample,
## not the full positional signal matrix. Therefore, the accompanying plot
## is the mean signal by group and sample rather than a metaprofile across bp.
profile_dt <- heatmap_dt[
  ,
  .(
    mean_signal = mean(signal, na.rm = TRUE),
    se_signal = sd(signal, na.rm = TRUE) / sqrt(.N)
  ),
  by = .(group, sample)
]

profile_plot <- ggplot(
  profile_dt,
  aes(
    x = sample,
    y = mean_signal,
    group = sample
  )
) +
  geom_point(size = 2.4) +
  geom_errorbar(
    aes(
      ymin = mean_signal - se_signal,
      ymax = mean_signal + se_signal
    ),
    width = 0.15
  ) +
  facet_grid(
    group ~ .,
    scales = "free_y"
  ) +
  labs(
    x = NULL,
    y = "Mean signal"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    panel.spacing.y = grid::unit(0.08, "lines"),
    plot.margin = margin(5, 5, 5, 5)
  )

combined <- heatmap_plot + profile_plot +
  plot_layout(widths = c(2.2, 1))

ggsave(
  paste0(out_prefix, ".pdf"),
  combined,
  width = 10,
  height = 10
)

ggsave(
  paste0(out_prefix, ".png"),
  combined,
  width = 10,
  height = 10,
  dpi = 300
)

write.csv(
  profile_dt,
  paste0(out_prefix, "_group_mean_signal.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  paste0(out_prefix, "_sessionInfo.txt")
)

cat("Saved:\n")
cat(paste0(out_prefix, ".pdf\n"))
cat(paste0(out_prefix, ".png\n"))
