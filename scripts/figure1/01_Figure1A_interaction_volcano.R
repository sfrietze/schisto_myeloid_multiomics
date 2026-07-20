#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(DESeq2)
    library(readr)
    library(dplyr)
    library(tibble)
    library(ggplot2)
    library(ggrepel)
})

repo <- normalizePath(".")

count_file <- file.path(repo, "data/rnaseq/raw/keke_counts.csv")
results_dir <- file.path(repo, "results/rnaseq/figure1")
figure_dir <- file.path(repo, "figures/main/Figure1")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

counts_input <- read_csv(
    count_file,
    show_col_types = FALSE,
    name_repair = "minimal"
)

gene_names <- trimws(as.character(counts_input[[1]]))
counts_df <- counts_input[, -1, drop = FALSE]

counts_df[] <- lapply(
    counts_df,
    function(x) suppressWarnings(as.numeric(x))
)

valid_rows <- !is.na(gene_names) &
    gene_names != "" &
    rowSums(is.na(counts_df)) == 0

message(
    "Removing ",
    sum(!valid_rows),
    " blank or nonnumeric rows from the count file."
)

gene_names <- gene_names[valid_rows]
counts_df <- counts_df[valid_rows, , drop = FALSE]

count_matrix <- as.matrix(counts_df)
storage.mode(count_matrix) <- "numeric"

rownames(count_matrix) <- make.unique(gene_names)
count_matrix <- round(count_matrix)

if (anyNA(count_matrix)) {
    stop("NA values remain in the count matrix after cleaning.")
}

if (any(count_matrix < 0)) {
    stop("Negative values were found in the count matrix.")
}

message(
    "Imported ",
    nrow(count_matrix),
    " genes and ",
    ncol(count_matrix),
    " samples."
)

sample_names <- colnames(count_matrix)

sample_metadata <- tibble(
    sample = sample_names,
    sex = sub("_.*$", "", sample_names),
    condition = sub("^[FM]_([UI])_.*$", "\\1", sample_names)
) |>
    mutate(
        sex = factor(sex, levels = c("F", "M")),
        condition = factor(condition, levels = c("U", "I"))
    ) |>
    column_to_rownames("sample")

stopifnot(
    identical(
        rownames(sample_metadata),
        colnames(count_matrix)
    )
)

dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sample_metadata,
    design = ~ sex + condition + sex:condition
)

dds <- DESeq(dds)

print(resultsNames(dds))

interaction_name <- "sexM.conditionI"

if (!interaction_name %in% resultsNames(dds)) {
    stop(
        "Expected coefficient not found: ",
        interaction_name,
        "\nAvailable coefficients:\n",
        paste(resultsNames(dds), collapse = "\n")
    )
}

res <- results(
    dds,
    name = interaction_name
)

volcano_df <- as.data.frame(res) |>
    rownames_to_column("gene") |>
    mutate(
        interaction_log2FC = log2FoldChange,
        neglog10_padj = -log10(padj),
        direction = case_when(
            !is.na(padj) &
                padj < 0.05 &
                interaction_log2FC > 1 ~ "Up in Male",
            !is.na(padj) &
                padj < 0.05 &
                interaction_log2FC < -1 ~ "Up in Female",
            TRUE ~ "Not significant"
        )
    )

finite_max <- max(
    volcano_df$neglog10_padj[is.finite(volcano_df$neglog10_padj)],
    na.rm = TRUE
)

volcano_df <- volcano_df |>
    mutate(
        neglog10_padj = ifelse(
            is.infinite(neglog10_padj),
            finite_max,
            neglog10_padj
        )
    )

n_female <- sum(volcano_df$direction == "Up in Female")
n_male <- sum(volcano_df$direction == "Up in Male")

message("Up in Female: ", n_female)
message("Up in Male: ", n_male)

label_genes <- bind_rows(
    volcano_df |>
        filter(direction == "Up in Female") |>
        arrange(padj) |>
        slice_head(n = 10),
    volcano_df |>
        filter(direction == "Up in Male") |>
        arrange(padj) |>
        slice_head(n = 10)
)

plot_df <- volcano_df |>
    filter(
        !is.na(interaction_log2FC),
        !is.na(neglog10_padj),
        is.finite(interaction_log2FC),
        is.finite(neglog10_padj)
    )

volcano_plot <- ggplot(
    plot_df,
    aes(
        x = interaction_log2FC,
        y = neglog10_padj
    )
) +
    geom_point(
        data = filter(plot_df, direction == "Not significant"),
        color = "gray80",
        alpha = 0.65,
        size = 1.3
    ) +
    geom_point(
        data = filter(plot_df, direction != "Not significant"),
        aes(color = direction),
        alpha = 0.85,
        size = 1.7
    ) +
    geom_vline(
        xintercept = c(-1, 1),
        linetype = "dashed",
        linewidth = 0.4
    ) +
    geom_hline(
        yintercept = -log10(0.05),
        linetype = "dotted",
        linewidth = 0.4
    ) +
    geom_text_repel(
        data = label_genes,
        aes(label = gene),
        color = "black",
        size = 3,
        min.segment.length = 0,
        max.overlaps = Inf,
        seed = 123
    ) +
    scale_color_manual(
        values = c(
            "Up in Female" = "#3B6FB6",
            "Up in Male" = "#E64B35"
        ),
        labels = c(
            "Up in Female" = paste0("Up in Female (", n_female, ")"),
            "Up in Male" = paste0("Up in Male (", n_male, ")")
        )
    ) +
    labs(
        x = "Sex-biased infection response (log2 fold change)",
        y = expression(-log[10] ~ "adjusted p-value"),
        color = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        legend.position = c(0.72, 0.83),
        legend.background = element_blank()
    )

write_csv(
    volcano_df,
    file.path(
        results_dir,
        "Figure1A_interaction_results.csv"
    )
)

saveRDS(
    dds,
    file.path(
        results_dir,
        "Figure1A_dds.rds"
    )
)

ggsave(
    file.path(
        figure_dir,
        "Figure1A_interaction_volcano.pdf"
    ),
    plot = volcano_plot,
    width = 5,
    height = 5,
    useDingbats = FALSE
)

ggsave(
    file.path(
        figure_dir,
        "Figure1A_interaction_volcano.png"
    ),
    plot = volcano_plot,
    width = 5,
    height = 5,
    dpi = 600
)
