suppressPackageStartupMessages({
  library(GenomicRanges)
  library(data.table)
  library(ssvTracks)
  library(seqsetvis)
  library(ggplot2)
})

source("config.R")

output_dir <- FIGURE1_RESULTS_DIR
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Nrf1 locus, mm10.
custom_region <- GRanges(
  seqnames = "chr6",
  ranges = IRanges(
    start = 30034000,
    end = 30058000
  )
)
genome(custom_region) <- "mm10"

cfg_bw <- data.table(
  sample = c(
    "UF_fem_inf",
    "IF_fem_inf",
    "UF_male_inf",
    "IF_male_inf"
  ),
  file = c(
    file.path(BIGWIG_DIR, "UF.pooled.bigWig"),
    file.path(BIGWIG_DIR, "IF.pooled.bigWig"),
    file.path(BIGWIG_DIR, "UNINF.pooled.bigWig"),
    file.path(BIGWIG_DIR, "INF.pooled.bigWig")
  )
)

missing_files <- cfg_bw[!file.exists(file), file]

if (length(missing_files) > 0) {
  stop(
    "Missing bigWig files:\n",
    paste(missing_files, collapse = "\n")
  )
}

cfg_bw[, group := fifelse(
  grepl("fem", sample),
  "Female",
  "Male"
)]

cfg_bw[, group := factor(
  group,
  levels = c("Female", "Male")
)]

cfg_bw[, id := factor(
  sample,
  levels = c(
    "UF_fem_inf",
    "IF_fem_inf",
    "UF_male_inf",
    "IF_male_inf"
  )
)]

cfg_bw[, nspline := 10]
cfg_bw[, nMovingAverage := 20]
cfg_bw[, ceiling_value := 350]

signal_colors <- c(
  UF_fem_inf = "gray40",
  IF_fem_inf = "darkred",
  UF_male_inf = "gray40",
  IF_male_inf = "darkred"
)

plot_overlay <- track_chip(
  signal_files = cfg_bw,
  query_gr = custom_region,
  fetch_fun = ssvFetchBigwig,
  fill_VAR = "id",
  color_VAR = "id",
  facet_VAR = "group",
  fill_mapping = signal_colors,
  color_mapping = signal_colors,
  nspline = 10,
  nMovingAverage = 20,
  ceiling_value = 350,
  y_label = "ATAC-seq (RPKM)"
) +
  labs(
    x = NULL,
    y = "ATAC-seq (RPKM)"
  ) +
  theme(
    strip.text = element_text(
      size = 13,
      face = "plain",
      hjust = 0
    ),
    legend.position = "right"
  )

gtf_file <- path.expand(
  "~/reference/gencode.vM25.annotation.gtf"
)

if (!file.exists(gtf_file)) {
  stop("Missing GTF file: ", gtf_file)
}

p_gene <- track_gene_reference(
  ref = gtf_file,
  query_gr = custom_region,
  show_tss = TRUE,
  tss_arrow_size = 0.2,
  minus_strand_color = "black",
  plus_strand_color = "black"
)

final_plot <- assemble_tracks(
  list(
    plot_overlay,
    p_gene
  ),
  query_gr = custom_region,
  rel_heights = c(1.5, 1)
)

output_file <- file.path(
  output_dir,
  "Figure1H_Nrf1_ATAC_overlay.pdf"
)

ggsave(
  filename = output_file,
  plot = final_plot,
  width = 10,
  height = 4,
  units = "in",
  device = cairo_pdf
)

message("Saved: ", output_file)
