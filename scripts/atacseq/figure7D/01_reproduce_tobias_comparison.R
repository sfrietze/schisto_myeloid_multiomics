#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

source("config.R")

data_dir <- TOBIAS_DIR
results_dir <- FIGURE7D_RESULTS_DIR

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

sea_file <- file.path(
  data_dir,
  "SEA_vs_NOSEA_bindetect_results.txt"
)

ifuf_file <- file.path(
  data_dir,
  "IF_vs_UF_bindetect_results.txt"
)

sea_top_file <- file.path(
  data_dir,
  "SEA_vs_NOSEA_top_motifs.txt"
)

ifuf_top_file <- file.path(
  data_dir,
  "IF_vs_UF_top_motifs.txt"
)

required_files <- c(sea_file, ifuf_file)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required input file(s):\n",
    paste(missing_files, collapse = "\n")
  )
}

sea <- read_tsv(
  sea_file,
  show_col_types = FALSE
)

ifuf <- read_tsv(
  ifuf_file,
  show_col_types = FALSE
)

required_sea_columns <- c(
  "output_prefix",
  "name",
  "motif_id",
  "SEA_NOSEA_change",
  "SEA_NOSEA_pvalue",
  "SEA_NOSEA_highlighted"
)

required_ifuf_columns <- c(
  "output_prefix",
  "name",
  "motif_id",
  "IF_UF_change",
  "IF_UF_pvalue",
  "IF_UF_highlighted"
)

missing_sea_columns <- setdiff(required_sea_columns, colnames(sea))
missing_ifuf_columns <- setdiff(required_ifuf_columns, colnames(ifuf))

if (length(missing_sea_columns) > 0) {
  stop(
    "Missing SEA/NOSEA column(s): ",
    paste(missing_sea_columns, collapse = ", ")
  )
}

if (length(missing_ifuf_columns) > 0) {
  stop(
    "Missing IF/UF column(s): ",
    paste(missing_ifuf_columns, collapse = ", ")
  )
}

as_highlighted <- function(x) {
  tolower(as.character(x)) %in% c(
    "true",
    "t",
    "1",
    "yes",
    "y"
  )
}

sea_plot <- sea %>%
  transmute(
    output_prefix,
    motif_id,
    motif_name = name,
    SEA_NOSEA_change,
    SEA_NOSEA_pvalue,
    SEA_NOSEA_highlighted = as_highlighted(
      SEA_NOSEA_highlighted
    )
  )

ifuf_plot <- ifuf %>%
  transmute(
    output_prefix,
    motif_id,
    IF_UF_change,
    IF_UF_pvalue,
    IF_UF_highlighted = as_highlighted(
      IF_UF_highlighted
    )
  )

plot_df <- inner_join(
  sea_plot,
  ifuf_plot,
  by = c("output_prefix", "motif_id")
) %>%
  filter(
    is.finite(SEA_NOSEA_change),
    is.finite(IF_UF_change)
  ) %>%
  mutate(
    highlighted = SEA_NOSEA_highlighted |
      IF_UF_highlighted,
    quadrant = case_when(
      IF_UF_change >= 0 &
        SEA_NOSEA_change >= 0 ~ "High in Both",
      IF_UF_change < 0 &
        SEA_NOSEA_change < 0 ~ "Low in Both",
      IF_UF_change >= 0 &
        SEA_NOSEA_change < 0 ~ "Up in IF, Down in SEA",
      IF_UF_change < 0 &
        SEA_NOSEA_change >= 0 ~ "Up in SEA, Down in IF"
    ),
    plot_group = if_else(
      highlighted,
      quadrant,
      "Not highlighted"
    )
  )

read_top_motifs <- function(path) {
  if (!file.exists(path)) {
    return(character())
  }

  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]

  unique(unlist(strsplit(lines, "[\t, ]+")))
}

label_terms <- unique(c(
  read_top_motifs(sea_top_file),
  read_top_motifs(ifuf_top_file)
))

label_df <- plot_df %>%
  filter(
    output_prefix %in% label_terms |
      motif_id %in% label_terms |
      motif_name %in% label_terms
  )

plot_df$plot_group <- factor(
  plot_df$plot_group,
  levels = c(
    "Not highlighted",
    "High in Both",
    "Low in Both",
    "Up in IF, Down in SEA",
    "Up in SEA, Down in IF"
  )
)

group_colors <- c(
  "Not highlighted" = "grey80",
  "High in Both" = "mediumpurple2",
  "Low in Both" = "blue",
  "Up in IF, Down in SEA" = "red",
  "Up in SEA, Down in IF" = "orange"
)

p <- ggplot(
  plot_df,
  aes(
    x = IF_UF_change,
    y = SEA_NOSEA_change
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_point(
    aes(color = plot_group),
    size = 2.2,
    alpha = 0.85
  ) +
  geom_text_repel(
    data = label_df,
    aes(label = motif_name),
    color = "black",
    size = 3.2,
    box.padding = 0.4,
    point.padding = 0.25,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 7
  ) +
  scale_color_manual(
    values = group_colors,
    drop = FALSE
  ) +
  labs(
    x = "IF vs UF Change",
    y = "SEA vs NOSEA Change",
    color = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    axis.title = element_text(face = "bold")
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        size = 3,
        alpha = 1
      )
    )
  )

write_csv(
  plot_df,
  file.path(
    results_dir,
    "Supplementary_Data_Figure7D_TOBIAS_comparison.csv"
  )
)

ggsave(
  filename = file.path(
    results_dir,
    "Figure7D_TOBIAS_comparison.pdf"
  ),
  plot = p,
  width = 7,
  height = 6,
  units = "in"
)

ggsave(
  filename = file.path(
    results_dir,
    "Figure7D_TOBIAS_comparison.png"
  ),
  plot = p,
  width = 7,
  height = 6,
  units = "in",
  dpi = 300
)

message("Merged motifs: ", nrow(plot_df))
message("Highlighted motifs: ", sum(plot_df$highlighted))
message("Labeled motifs: ", nrow(label_df))
message("Results written to: ", results_dir)
