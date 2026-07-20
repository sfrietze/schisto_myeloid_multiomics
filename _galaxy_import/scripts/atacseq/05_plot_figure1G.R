#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(chiptsne2)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(cowplot)
  library(ggplotify)
})

source("config.R")

ct2 <- readRDS(
  file.path(
    RESULTS_DIR,
    "atacseq",
    "figure1",
    "Figure1G_heatmap_object.rds"
  )
)

flatten_plot_list <- function(x) {
  if (
    inherits(x, "ggplot") ||
    inherits(x, "grob") ||
    inherits(x, "gtable")
  ) {
    return(list(x))
  }

  if (is.list(x)) {
    return(unlist(
      lapply(x, flatten_plot_list),
      recursive = FALSE
    ))
  }

  list()
}

patched_get_rel_widths <- function(my_plots, sync_width = TRUE) {
  my_plots <- flatten_plot_list(my_plots)

  stopifnot(length(my_plots) > 0)

  plot_grobs <- lapply(
    my_plots,
    function(x) {
      if (inherits(x, "ggplot")) {
        ggplot2::ggplotGrob(x)
      } else {
        x
      }
    }
  )

  if (!sync_width || length(plot_grobs) == 1) {
    return(rep(1, length(plot_grobs)))
  }

  widths <- vapply(
    plot_grobs,
    function(x) {
      if (!is.null(x$widths)) {
        sum(as.numeric(x$widths))
      } else {
        1
      }
    },
    numeric(1)
  )

  widths / max(widths)
}

assignInNamespace(
  "get_rel_widths",
  patched_get_rel_widths,
  ns = "chiptsne2"
)

group_levels <- c(
  "Male > Female (IF)",
  "Female > Male (IF)",
  "Male > Female (UF)",
  "Female > Male (UF)"
)

rowData(ct2)$group <- factor(
  rowData(ct2)$group,
  levels = group_levels
)


published_group_order <- c(
    "Female > Male (IF)",
    "Female > Male (UF)",
    "Male > Female (IF)",
    "Male > Female (UF)"
)

group_order_index <- match(
    as.character(rowData(ct2)$group),
    published_group_order
)

if (anyNA(group_order_index)) {
    stop(
        "Unexpected group labels: ",
        paste(
            unique(as.character(rowData(ct2)$group)[is.na(group_order_index)]),
            collapse = ", "
        )
    )
}

ct2 <- ct2[order(group_order_index), ]

rowData(ct2)$group <- factor(
    as.character(rowData(ct2)$group),
    levels = published_group_order
)

ht <- plotSignalHeatmap(
  ct2,
  group_VARS = "group",
  sort_strategy = "none",
  relative_heatmap_width = 0.4,
  heatmap_fill_limits = c(0, 150),
  heatmap_colors = c("#577AB2", "#FFFFE0", "#BC412B")
)

cd <- colData(ct2)
cd$sample <- rownames(cd)
cd$sex_group <- ifelse(
  grepl("^UF_fem|^IF_fem", cd$sample),
  "Female",
  "Male"
)
colData(ct2) <- cd

ct2$sample <- factor(
  ct2$sample,
  levels = c(
    "UF_fem_inf",
    "IF_fem_inf",
    "UF_male_inf",
    "IF_male_inf"
  )
)

ct2$sex_group <- factor(
  ct2$sex_group,
  levels = c("Female", "Male")
)

custom_colors <- c(
  UF_fem_inf  = "gray",
  IF_fem_inf  = "darkred",
  UF_male_inf = "gray",
  IF_male_inf = "darkred"
)

lp <- plotSignalLinePlot(
  ct2,
  group_VAR = "group",
  color_VAR = "sample",
  facet_VAR = "sex_group",
  moving_average_window = 5,
  n_splines = 5,
  linewidth = 1.2
) +
  scale_color_manual(values = custom_colors) +
  facet_grid(group ~ sex_group, scales = "free_y") +
  labs(
    x = "Distance from Region Center (bp)",
    y = "Signal Intensity"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )


dir.create(
    file.path("results", "atacseq", "figure1"),
    recursive = TRUE,
    showWarnings = FALSE
)

pdf(
    file.path(
        "results",
        "atacseq",
        "figure1",
        "Figure1G_heatmap.pdf"
    ),
    width = 10,
    height = 9
)
print(ht)
dev.off()

pdf(
    file.path(
        "results",
        "atacseq",
        "figure1",
        "Figure1G_lineplot.pdf"
    ),
    width = 8,
    height = 9
)
print(lp)
dev.off()

