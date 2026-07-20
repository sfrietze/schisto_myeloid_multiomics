#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

root <- normalizePath("~/projects/keke")

input_dir <- file.path(
  root,
  "reproducibility/data/atacseq/figure7C"
)

output_dir <- file.path(
  root,
  "reproducibility/results/atacseq/figure7C"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

bindetect_file <- file.path(
  input_dir,
  "bindetect_results.txt"
)

top_motifs_file <- file.path(
  input_dir,
  "top_motifs.txt"
)

if (!file.exists(bindetect_file)) {
  stop("Missing input file: ", bindetect_file)
}

if (!file.exists(top_motifs_file)) {
  stop("Missing input file: ", top_motifs_file)
}

bindetect <- read_tsv(
  bindetect_file,
  show_col_types = FALSE
)

top_motifs <- read_lines(
  top_motifs_file
)

required_columns <- c(
  "output_prefix",
  "name",
  "SEA_NOSEA_change",
  "SEA_NOSEA_pvalue"
)

missing_columns <- setdiff(
  required_columns,
  colnames(bindetect)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

q10 <- quantile(
  bindetect$SEA_NOSEA_change,
  probs = 0.10,
  na.rm = TRUE
)

q90 <- quantile(
  bindetect$SEA_NOSEA_change,
  probs = 0.90,
  na.rm = TRUE
)

bindetect <- bindetect %>%
  mutate(
    neglog10p = -log10(
      pmax(
        SEA_NOSEA_pvalue,
        .Machine$double.xmin
      )
    ),
    category = case_when(
      SEA_NOSEA_change <= q10 ~ "Negative (Top 10%)",
      SEA_NOSEA_change >= q90 ~ "Positive (Top 10%)",
      TRUE ~ "Other"
    ),
    label = if_else(
      output_prefix %in% top_motifs,
      name,
      NA_character_
    )
  )

category_levels <- c(
  "Other",
  "Negative (Top 10%)",
  "Positive (Top 10%)"
)

bindetect$category <- factor(
  bindetect$category,
  levels = category_levels
)

write_csv(
  bindetect,
  file.path(
    output_dir,
    "Figure7C_TOBIAS_volcano_data.csv"
  )
)

message("Rows: ", nrow(bindetect))
message("10th percentile cutoff: ", signif(q10, 6))
message("90th percentile cutoff: ", signif(q90, 6))
message(
  "Labeled motifs: ",
  sum(!is.na(bindetect$label))
)

p <- ggplot(
  bindetect,
  aes(
    x = SEA_NOSEA_change,
    y = neglog10p
  )
) +
  geom_point(
    aes(color = category),
    size = 2.7,
    alpha = 0.9
  ) +
  geom_text_repel(
    data = bindetect %>% filter(!is.na(label)),
    aes(label = label),
    color = "black",
    seed = 123,
    size = 5,
    box.padding = 0.35,
    point.padding = 0.20,
    min.segment.length = 0,
    segment.color = "black",
    segment.size = 0.3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Other" = "grey50",
      "Negative (Top 10%)" = "blue",
      "Positive (Top 10%)" = "red"
    ),
    breaks = c(
      "Negative (Top 10%)",
      "Positive (Top 10%)"
    ),
    name = "category"
  ) +
  scale_x_continuous(
    limits = c(-0.20, 0.27),
    breaks = seq(-0.2, 0.2, 0.1)
  ) +
  scale_y_continuous(
    limits = c(0, 122),
    breaks = seq(0, 120, 30)
  ) +
  labs(
    title = "Differential ATAC-seq footprint analysis (SEA vs NOSEA)",
    x = "SEA vs NOSEA Change",
    y = expression(-log[10](P-value))
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 14),
    axis.text = element_text(
      size = 12,
      color = "black"
    )
  )

ggsave(
  filename = file.path(
    output_dir,
    "Figure7C_TOBIAS_volcano.pdf"
  ),
  plot = p,
  width = 8.5,
  height = 6.5
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure7C_TOBIAS_volcano.png"
  ),
  plot = p,
  width = 8.5,
  height = 6.5,
  dpi = 600
)

message("Saved outputs to: ", output_dir)
