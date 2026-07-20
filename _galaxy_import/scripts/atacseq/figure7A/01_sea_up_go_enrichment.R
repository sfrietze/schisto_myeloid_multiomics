suppressPackageStartupMessages({
  library(DESeq2)
  library(GenomicRanges)
  library(ChIPseeker)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(ggplot2)
})

root <- normalizePath("~/projects/keke")
outdir <- file.path(root, "reproducibility/results/atacseq/figure7A")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------
# 1. Load the SEA-specific consensus-count object
# ------------------------------------------------------------------
dds <- readRDS(file.path(root, "DAR/female_sea_dds.rds"))

# ------------------------------------------------------------------
# 2. Use the four samples retained in the original repository analysis
# ------------------------------------------------------------------
samples_to_keep <- c(
  "NOSEA_REP2.mLb.clN",
  "NOSEA_REP3.mLb.clN",
  "SEA_REP1.mLb.clN",
  "SEA_REP2.mLb.clN"
)

missing_samples <- setdiff(samples_to_keep, colnames(dds))

if (length(missing_samples) > 0) {
  stop(
    "Missing expected samples: ",
    paste(missing_samples, collapse = ", ")
  )
}

dds_sea <- dds[, samples_to_keep]

# The saved object uses inherited IF/UF labels, so reconstruct the
# experimental condition explicitly from the sample names.
dds_sea$condition <- factor(
  ifelse(
    grepl("^SEA_", colnames(dds_sea)),
    "SEA",
    "NOSEA"
  ),
  levels = c("NOSEA", "SEA")
)

design(dds_sea) <- ~ condition

message("Samples used:")
print(
  data.frame(
    sample = colnames(dds_sea),
    condition = dds_sea$condition
  )
)

# ------------------------------------------------------------------
# 3. Run DESeq2
# ------------------------------------------------------------------
dds_sea <- DESeq(dds_sea)

res <- results(
  dds_sea,
  contrast = c("condition", "SEA", "NOSEA"),
  alpha = 0.05
)

res_df <- as.data.frame(res)
res_df$peak_id <- rownames(res_df)

write.csv(
  res_df,
  file.path(outdir, "female_SEA_vs_NOSEA_DESeq2_all.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# 4. Use all SEA-vs-NOSEA differential peaks
# ------------------------------------------------------------------
dars <- readRDS(
  file.path(root, "DAR/female_sea_2v2_DESeq2_DARs.rds")
)

message("SEA-vs-NOSEA DARs used for enrichment: ", length(dars))

saveRDS(
  dars,
  file.path(outdir, "female_SEA_vs_NOSEA_DARs.rds")
)

# ------------------------------------------------------------------
# 5. Annotate DARs to nearest genes
# ------------------------------------------------------------------
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

peak_annotation <- annotatePeak(
  dars,
  TxDb = txdb,
  annoDb = "org.Mm.eg.db",
  tssRegion = c(-3000, 3000),
  verbose = FALSE
)

annotation_df <- as.data.frame(peak_annotation)

write.csv(
  annotation_df,
  file.path(outdir, "female_SEA_vs_NOSEA_DARs_annotated.csv"),
  row.names = FALSE
)

dar_gene_ids <- unique(na.omit(as.character(annotation_df$geneId)))

message("Unique genes assigned to DARs: ", length(dar_gene_ids))

# ------------------------------------------------------------------
# 7. Export the annotated differential-accessibility peaks
# ------------------------------------------------------------------
supp_dars <- annotation_df

preferred_columns <- c(
  "seqnames",
  "start",
  "end",
  "width",
  "strand",
  "baseMean",
  "log2FC",
  "padj",
  "annotation",
  "geneChr",
  "geneStart",
  "geneEnd",
  "geneLength",
  "geneStrand",
  "geneId",
  "transcriptId",
  "distanceToTSS",
  "ENSEMBL",
  "SYMBOL",
  "GENENAME"
)

preferred_columns <- intersect(
  preferred_columns,
  colnames(supp_dars)
)

remaining_columns <- setdiff(
  colnames(supp_dars),
  preferred_columns
)

supp_dars <- supp_dars[
  ,
  c(preferred_columns, remaining_columns),
  drop = FALSE
]

write.csv(
  supp_dars,
  file.path(
    outdir,
    "Supplementary_Data_Figure7A_DARs_annotated.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------
# 8. GO Biological Process enrichment
#
# The complete result table is retained. Redundant GO terms are
# collapsed only for visualization.
# ------------------------------------------------------------------
ego <- enrichGO(
  gene = dar_gene_ids,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  minGSSize = 5,
  readable = TRUE
)

ego_full_df <- as.data.frame(ego)

if (nrow(ego_full_df) == 0) {
  stop("No GO Biological Process terms were returned")
}

ego_full_df <- ego_full_df[
  order(
    ego_full_df$p.adjust,
    ego_full_df$pvalue,
    -ego_full_df$Count
  ),
  ,
  drop = FALSE
]

write.csv(
  ego_full_df,
  file.path(
    outdir,
    "Supplementary_Data_Figure7A_GO_BP_full.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------
# 9. Collapse semantically redundant GO Biological Process terms
#
# cutoff = 0.70 groups terms with high semantic similarity.
# The term with the smallest adjusted p-value is retained.
# ------------------------------------------------------------------
ego_simplified <- simplify(
  ego,
  cutoff = 0.70,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang"
)

ego_simplified_df <- as.data.frame(ego_simplified)

if (nrow(ego_simplified_df) == 0) {
  stop("No GO terms remained after redundancy reduction")
}

ego_simplified_df <- ego_simplified_df[
  order(
    ego_simplified_df$p.adjust,
    ego_simplified_df$pvalue,
    -ego_simplified_df$Count
  ),
  ,
  drop = FALSE
]

write.csv(
  ego_simplified_df,
  file.path(
    outdir,
    "Supplementary_Data_Figure7A_GO_BP_nonredundant.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------
# 10. Select and export the exact terms shown in Figure 7A
# ------------------------------------------------------------------
n_terms <- min(10, nrow(ego_simplified_df))

plot_df <- ego_simplified_df[
  seq_len(n_terms),
  ,
  drop = FALSE
]

ratio_parts <- strsplit(
  as.character(plot_df$GeneRatio),
  "/",
  fixed = TRUE
)

plot_df$GeneRatioNumeric <- vapply(
  ratio_parts,
  function(x) {
    if (length(x) != 2) {
      return(NA_real_)
    }

    as.numeric(x[1]) / as.numeric(x[2])
  },
  numeric(1)
)

plot_df$minus_log10_adjusted_pvalue <- -log10(
  pmax(plot_df$p.adjust, .Machine$double.xmin)
)

plot_df$Description <- factor(
  plot_df$Description,
  levels = rev(plot_df$Description)
)

write.csv(
  plot_df,
  file.path(
    outdir,
    "Supplementary_Data_Figure7A_GO_BP_plotted_terms.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------
# 11. Revised Figure 7A
# ------------------------------------------------------------------
p <- ggplot(
  plot_df,
  aes(
    x = minus_log10_adjusted_pvalue,
    y = Description,
    size = Count,
    fill = GeneRatioNumeric
  )
) +
  geom_point(
    shape = 21,
    color = "black",
    stroke = 0.8,
    alpha = 0.95
  ) +
  scale_fill_gradientn(
    colors = c(
      "blue",
      "purple",
      "deeppink",
      "red"
    ),
    name = "GeneRatio"
  ) +
  scale_size_continuous(
    name = "Count",
    range = c(3, 10)
  ) +
  labs(
    x = expression(-log[10]("adjusted p-value")),
    y = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_blank(),
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.6
    ),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      color = "grey25",
      size = 11
    ),
    axis.text.y = element_text(
      color = "grey30",
      size = 11
    ),
    axis.title.x = element_text(
      color = "black",
      size = 12
    ),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.margin = margin(
      t = 12,
      r = 12,
      b = 12,
      l = 18
    )
  )

ggsave(
  file.path(
    outdir,
    "Figure7A_GO_BP_nonredundant.pdf"
  ),
  plot = p,
  width = 8.5,
  height = 6.5,
  device = cairo_pdf
)

ggsave(
  file.path(
    outdir,
    "Figure7A_GO_BP_nonredundant.png"
  ),
  plot = p,
  width = 8.5,
  height = 6.5,
  dpi = 600
)

message("Figure 7A outputs written to: ", outdir)

message(
  "Annotated DARs: ",
  nrow(supp_dars)
)

message(
  "Complete GO BP terms: ",
  nrow(ego_full_df)
)

message(
  "Nonredundant GO BP terms: ",
  nrow(ego_simplified_df)
)

message(
  "Terms plotted: ",
  nrow(plot_df)
)
