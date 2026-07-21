#!/usr/bin/env Rscript
# See README.md for upstream ATAC-seq preprocessing provenance.

suppressPackageStartupMessages({
  library(DESeq2)
  library(GenomicRanges)
  library(SummarizedExperiment)
})

input_file <- file.path(
  "data",
  "atacseq",
  "processed",
  "ATAC_consensus_counts_metadata_bundle.rds"
)

output_dir <- file.path(
  "results",
  "atacseq",
  "figure1"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

bundle <- readRDS(input_file)

required_objects <- c(
  "counts",
  "peaks",
  "metadata"
)

missing_objects <- setdiff(
  required_objects,
  names(bundle)
)

if (length(missing_objects) > 0) {
  stop(
    "Missing bundle components: ",
    paste(missing_objects, collapse = ", ")
  )
}

count_matrix <- as.matrix(bundle$counts)
peak_table <- as.data.frame(bundle$peaks)
sample_metadata <- as.data.frame(bundle$metadata)

storage.mode(count_matrix) <- "numeric"

if (!"sample" %in% colnames(sample_metadata)) {
  sample_metadata$sample <- rownames(sample_metadata)
}

if (!all(colnames(count_matrix) %in% sample_metadata$sample)) {
  stop(
    "Count-matrix sample names do not match metadata sample names."
  )
}

sample_metadata <- sample_metadata[
  match(
    colnames(count_matrix),
    sample_metadata$sample
  ),
  ,
  drop = FALSE
]

rownames(sample_metadata) <- sample_metadata$sample

condition_column <- intersect(
  c("condition", "Condition"),
  colnames(sample_metadata)
)[1]

sex_column <- intersect(
  c("sex", "Sex"),
  colnames(sample_metadata)
)[1]

if (is.na(condition_column)) {
  stop("Could not identify the condition column in the metadata.")
}

if (is.na(sex_column)) {
  stop("Could not identify the sex column in the metadata.")
}

sample_metadata$condition <- as.character(
  sample_metadata[[condition_column]]
)

sample_metadata$sex <- as.character(
  sample_metadata[[sex_column]]
)

sample_metadata$condition <- factor(
  sample_metadata$condition,
  levels = c("IF", "UF")
)

sample_metadata$sex <- factor(
  sample_metadata$sex,
  levels = c("F", "M")
)

if (any(is.na(sample_metadata$condition))) {
  stop(
    "Unexpected condition values: ",
    paste(
      unique(
        sample_metadata[[condition_column]]
      ),
      collapse = ", "
    )
  )
}

if (any(is.na(sample_metadata$sex))) {
  stop(
    "Unexpected sex values: ",
    paste(
      unique(
        sample_metadata[[sex_column]]
      ),
      collapse = ", "
    )
  )
}

keep_samples <- sample_metadata$condition %in% c("IF", "UF")

sample_metadata <- sample_metadata[
  keep_samples,
  ,
  drop = FALSE
]

count_matrix <- count_matrix[
  ,
  rownames(sample_metadata),
  drop = FALSE
]

count_matrix <- round(count_matrix)

if (any(is.na(count_matrix))) {
  stop("NA values remain in the count matrix.")
}

if (any(count_matrix < 0)) {
  stop("Negative values were found in the count matrix.")
}

required_peak_columns <- c(
  "chromosome",
  "start",
  "end"
)

missing_peak_columns <- setdiff(
  required_peak_columns,
  colnames(peak_table)
)

if (length(missing_peak_columns) > 0) {
  stop(
    "Missing peak columns: ",
    paste(missing_peak_columns, collapse = ", ")
  )
}

if (!"strand" %in% colnames(peak_table)) {
  peak_table$strand <- "*"
}

peak_ranges <- GRanges(
  seqnames = peak_table$chromosome,
  ranges = IRanges(
    start = as.integer(peak_table$start),
    end = as.integer(peak_table$end)
  ),
  strand = peak_table$strand
)

if (nrow(count_matrix) != length(peak_ranges)) {
  stop(
    "Peak count rows and genomic ranges differ: ",
    nrow(count_matrix),
    " versus ",
    length(peak_ranges)
  )
}

peak_ids <- if ("peak_id" %in% colnames(peak_table)) {
  as.character(peak_table$peak_id)
} else {
  paste0(
    peak_table$chromosome,
    ":",
    peak_table$start,
    "-",
    peak_table$end
  )
}

rownames(count_matrix) <- peak_ids
names(peak_ranges) <- peak_ids

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_metadata,
  design = ~ sex + condition + sex:condition,
  rowRanges = peak_ranges
)

# No additional count filtering was applied in the original analysis.
dds <- DESeq(dds)

output_file <- file.path(
  output_dir,
  "Figure1F_ATAC_interaction_dds.rds"
)

saveRDS(
  dds,
  output_file
)

write.csv(
  as.data.frame(colData(dds)),
  file.path(
    output_dir,
    "Figure1F_ATAC_sample_metadata.csv"
  ),
  row.names = TRUE
)

writeLines(
  capture.output({
    cat("Input file:\n")
    cat(input_file, "\n\n")

    cat("Design:\n")
    print(design(dds))

    cat("\nSex levels:\n")
    print(levels(dds$sex))

    cat("\nCondition levels:\n")
    print(levels(dds$condition))

    cat("\nSample groups:\n")
    print(table(dds$sex, dds$condition))

    cat("\nDimensions:\n")
    print(dim(dds))

    cat("\nCoefficient names:\n")
    print(resultsNames(dds))

    cat("\nSession information:\n")
    print(sessionInfo())
  }),
  file.path(
    output_dir,
    "Figure1F_ATAC_dds_build_provenance.txt"
  )
)

message("Saved: ", output_file)
message("Coefficient names:")
message(paste(resultsNames(dds), collapse = "\n"))
