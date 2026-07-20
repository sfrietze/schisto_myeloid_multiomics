suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ggplot2)
})

source("config.R")

outdir <- FIGURE7A_RESULTS_DIR

annotation_file <- file.path(
  outdir,
  "Supplementary_Data_Figure7A_DARs_annotated.csv"
)

if (!file.exists(annotation_file)) {
  annotation_file <- file.path(
    outdir,
    "female_SEA_vs_NOSEA_DARs_annotated.csv"
  )
}

if (!file.exists(annotation_file)) {
  stop("Annotated DAR CSV not found")
}

annotation_df <- read.csv(
  annotation_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

gene_ids <- unique(
  na.omit(as.character(annotation_df$geneId))
)

gene_ids <- gene_ids[gene_ids != ""]

message("Genes used for enrichment: ", length(gene_ids))

ego <- enrichGO(
  gene = gene_ids,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  minGSSize = 5,
  readable = TRUE
)

full_df <- as.data.frame(ego)
full_df <- full_df[
  order(full_df$p.adjust, full_df$pvalue),
  ,
  drop = FALSE
]

write.csv(
  full_df,
  file.path(
    outdir,
    "Supplementary_Data_Figure7A_GO_BP_full.csv"
  ),
  row.names = FALSE
)

sig_df <- full_df[
  !is.na(full_df$p.adjust) & full_df$p.adjust < 0.05,
  ,
  drop = FALSE
]

message("Significant GO terms before simplification: ", nrow(sig_df))

if (nrow(sig_df) == 0) {
  stop("No GO terms have adjusted p-value < 0.05")
}

ego_sig <- ego
ego_sig@result <- ego_sig@result[
  ego_sig@result$ID %in% sig_df$ID,
  ,
  drop = FALSE
]

ego_simple <- simplify(
  ego_sig,
  cutoff = 0.70,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang"
)

simple_df <- as.data.frame(ego_simple)
simple_df <- simple_df[
  order(simple_df$p.adjust, simple_df$pvalue),
  ,
  drop = FALSE
]

write.csv(
  simple_df,
  file.path(
    outdir,
    "Supplementary_Data_Figure7A_GO_BP_nonredundant.csv"
  ),
  row.names = FALSE
)

plot_df <- head(simple_df, 10)

ratio_split <- strsplit(
  as.character(plot_df$GeneRatio),
  "/",
  fixed = TRUE
)

plot_df$GeneRatioNumeric <- vapply(
  ratio_split,
  function(x) as.numeric(x[1]) / as.numeric(x[2]),
  numeric(1)
)

plot_df$minus_log10_padj <- -log10(plot_df$p.adjust)

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

p <- ggplot(
  plot_df,
  aes(
    x = minus_log10_padj,
    y = Description,
    size = Count,
    fill = GeneRatioNumeric
  )
) +
  geom_point(
    shape = 21,
    color = "black",
    stroke = 0.8
  ) +
  scale_fill_gradientn(
    colors = c("blue", "purple", "deeppink", "red"),
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
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.6
    ),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "grey25"),
    legend.position = "right",
    plot.margin = margin(12, 12, 12, 18)
  )

ggsave(
  file.path(outdir, "Figure7A_GO_BP_nonredundant.pdf"),
  p,
  width = 8.5,
  height = 6.5
)

ggsave(
  file.path(outdir, "Figure7A_GO_BP_nonredundant.png"),
  p,
  width = 8.5,
  height = 6.5,
  dpi = 600
)

message("Complete GO terms: ", nrow(full_df))
message("Nonredundant significant terms: ", nrow(simple_df))
message("Finished")
