#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(data.table)
  library(chiptsne2)
  library(rtracklayer)
  library(zoo)
})

source("config.R")

gr_dars_annotated <- readRDS(
  file.path(
    RESULTS_DIR,
    "atacseq",
    "figure1",
    "fig1_DARs_combined_annotated.rds"
  )
)

bw_files <- list(
  UF_fem_inf  = file.path(FEMALE_INF_DIR, "UF.pooled.bigWig"),
  IF_fem_inf  = file.path(FEMALE_INF_DIR, "IF.pooled.bigWig"),
  UF_male_inf = file.path(MALE_INF_DIR, "UNINF.pooled.bigWig"),
  IF_male_inf = file.path(MALE_INF_DIR, "INF.pooled.bigWig")
)

stopifnot(all(file.exists(unlist(bw_files))))

bw_cfg <- data.table(
  path = unlist(bw_files),
  name = names(bw_files)
)

fetch_cfg <- FetchConfig(
  bw_cfg,
  read_mode = "bigwig",
  view_size = 2500,
  window_size = 50
)

centerOnSignalMax <- function(gr, bw_path, window = 2500, smooth_window = 5) {
  gr_exp <- resize(gr, width = window, fix = "center")
  centers <- numeric(length(gr_exp))

  for (i in seq_along(gr_exp)) {
    region <- gr_exp[i]

    sig <- tryCatch(
      import(bw_path, which = region, as = "NumericList"),
      error = function(e) NULL
    )

    if (is.null(sig) || length(sig) == 0 || length(sig[[1]]) == 0) {
      centers[i] <- start(region) + window / 2
      next
    }

    smoothed <- zoo::rollmean(
      as.numeric(sig[[1]]),
      k = smooth_window,
      fill = NA,
      align = "center"
    )

    max_idx <- which.max(smoothed)

    centers[i] <- ifelse(
      is.na(max_idx),
      start(region) + window / 2,
      start(region) + max_idx - 1
    )
  }

  trim(
    GRanges(
      seqnames = seqnames(gr_exp),
      ranges = IRanges(
        start = centers - window / 2,
        width = window
      ),
      strand = strand(gr_exp)
    )
  )
}

gr_up_IF <- gr_dars_annotated[
  which(mcols(gr_dars_annotated)$group_IF == "Male > Female (IF)")
]

gr_up_UF <- gr_dars_annotated[
  which(mcols(gr_dars_annotated)$group_UF == "Male > Female (UF)")
]

gr_down_IF <- gr_dars_annotated[
  which(mcols(gr_dars_annotated)$group_IF == "Female > Male (IF)")
]

gr_down_UF <- gr_dars_annotated[
  which(mcols(gr_dars_annotated)$group_UF == "Female > Male (UF)")
]

gr_up <- c(gr_up_IF, gr_up_UF)
gr_down <- c(gr_down_IF, gr_down_UF)

names(gr_up) <- paste0("up_region_", seq_along(gr_up))
names(gr_down) <- paste0("down_region_", seq_along(gr_down))

gr_up_centered <- centerOnSignalMax(
  gr_up,
  bw_files$IF_fem_inf
)

gr_down_centered <- centerOnSignalMax(
  gr_down,
  bw_files$UF_fem_inf
)

names(gr_up_centered) <- names(gr_up)
names(gr_down_centered) <- names(gr_down)

gr_centered <- c(gr_up_centered, gr_down_centered)

olap <- findOverlaps(gr_centered, gr_dars_annotated)

mcols(gr_centered)$log2FC_IF <- NA_real_
mcols(gr_centered)$log2FC_UF <- NA_real_
mcols(gr_centered)$source <- NA_character_

mcols(gr_centered)$log2FC_IF[queryHits(olap)] <-
  mcols(gr_dars_annotated)$log2FC_IF[subjectHits(olap)]

mcols(gr_centered)$log2FC_UF[queryHits(olap)] <-
  mcols(gr_dars_annotated)$log2FC_UF[subjectHits(olap)]

mcols(gr_centered)$source[queryHits(olap)] <-
  mcols(gr_dars_annotated)$source[subjectHits(olap)]

priority_order <- c(
  "Male > Female (IF)",
  "Female > Male (IF)",
  "Male > Female (UF)",
  "Female > Male (UF)"
)

assign_primary_group <- function(source_str) {
  sources <- unlist(strsplit(source_str, ";"))

  for (group in priority_order) {
    if (group %in% sources) {
      return(group)
    }
  }

  NA_character_
}

gr_centered$primary_group <- vapply(
  gr_centered$source,
  assign_primary_group,
  character(1)
)

group_counts <- table(gr_centered$primary_group)
valid_groups <- names(group_counts[group_counts >= 20])

gr_filtered <- gr_centered[
  gr_centered$primary_group %in% valid_groups
]

ct2 <- ChIPtsne2.from_FetchConfig(
  fetch_cfg,
  query_gr = gr_filtered
)

rowData(ct2)$group <- gr_filtered$primary_group

ct2 <- sortRegions(
  ct2,
  sort_strategy = "sort",
  group_VAR = "group"
)

saveRDS(
  ct2,
  file.path(
    RESULTS_DIR,
    "atacseq",
    "figure1",
    "Figure1G_heatmap_object.rds"
  )
)
