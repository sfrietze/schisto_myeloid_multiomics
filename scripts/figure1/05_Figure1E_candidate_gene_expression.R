#!/usr/bin/env Rscript

suppressPackageStartupMessages({
library(DESeq2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
})

figure_dir <- file.path("figures", "main", "Figure1")
results_dir <- file.path("results", "rnaseq", "figure1")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

dds_candidates <- c(
file.path(results_dir, "Figure1A_dds.rds"),
file.path("results", "rnaseq", "Figure1A_dds.rds"),
file.path("results", "Figure1A_dds.rds")
)

dds_path <- dds_candidates[file.exists(dds_candidates)][1]

if (is.na(dds_path)) {
stop(
"Could not locate Figure1A_dds.rds. Checked:\n",
paste(dds_candidates, collapse = "\n")
)
}

dds <- readRDS(dds_path)

genes_to_plot <- c(
"Dio2",
"Fabp4",
"Hk3",
"Mgll",
"Scin",
"Slc1a3"
)

missing_genes <- setdiff(genes_to_plot, rownames(dds))

if (length(missing_genes) > 0) {
stop(
"Genes not found in the DESeq2 object:\n",
paste(missing_genes, collapse = "\n")
)
}

normalized_counts <- counts(
dds,
normalized = TRUE
)

metadata <- as.data.frame(colData(dds)) %>%
rownames_to_column("sample") %>%
mutate(
sex = as.character(sex),
condition = as.character(condition),
group = case_when(
sex == "F" & condition == "U" ~ "Fem_UF",
sex == "F" & condition == "I" ~ "Fem_IF",
sex == "M" & condition == "U" ~ "Male_UF",
sex == "M" & condition == "I" ~ "Male_IF",
TRUE ~ paste(sex, condition, sep = "_")
),
group = factor(
group,
levels = c(
"Fem_UF",
"Fem_IF",
"Male_UF",
"Male_IF"
)
)
)

plot_df <- normalized_counts[genes_to_plot, , drop = FALSE] %>%
as.data.frame() %>%
rownames_to_column("gene") %>%
pivot_longer(
cols = -gene,
names_to = "sample",
values_to = "normalized_count"
) %>%
left_join(
metadata %>%
select(sample, sex, condition, group),
by = "sample"
) %>%
mutate(
gene = factor(
gene,
levels = genes_to_plot
)
)

if (any(is.na(plot_df$group))) {
stop(
"One or more samples could not be assigned to the expected groups.\n",
paste(
unique(
plot_df$sample[is.na(plot_df$group)]
),
collapse = "\n"
)
)
}

write.csv(
plot_df,
file.path(
results_dir,
"Figure1E_candidate_gene_normalized_counts.csv"
),
row.names = FALSE
)

summary_df <- plot_df %>%
group_by(gene, group) %>%
summarise(
mean = mean(normalized_count),
sd = sd(normalized_count),
n = n(),
.groups = "drop"
)

write.csv(
summary_df,
file.path(
results_dir,
"Figure1E_candidate_gene_summary.csv"
),
row.names = FALSE
)

figure_1e <- ggplot(
summary_df,
aes(
x = group,
y = mean,
fill = group
)
) +
geom_col(
width = 0.68,
color = "black",
linewidth = 0.8
) +
geom_errorbar(
aes(
ymin = pmax(mean - sd, 0),
ymax = mean + sd
),
width = 0.18,
linewidth = 0.8
) +
geom_point(
data = plot_df,
aes(
x = group,
y = normalized_count
),
inherit.aes = FALSE,
position = position_jitter(
width = 0.08,
height = 0
),
size = 2.1,
alpha = 0.8
) +
facet_wrap(
~ gene,
ncol = 2,
scales = "free_y",
labeller = label_parsed
) +
scale_fill_manual(
values = c(
Fem_UF = "grey65",
Fem_IF = "firebrick4",
Male_UF = "grey65",
Male_IF = "firebrick4"
)
) +
scale_y_continuous(
expand = expansion(
mult = c(0, 0.08)
)
) +
labs(
x = NULL,
y = "Normalized Counts"
) +
theme_classic(
base_size = 12
) +
theme(
legend.position = "none",
strip.background = element_blank(),
strip.text = element_text(
face = "italic",
size = 13
),
axis.text.x = element_text(
angle = 50,
hjust = 1,
vjust = 1,
size = 10
),
axis.text.y = element_text(
size = 10
),
axis.title.y = element_text(
size = 13,
margin = margin(r = 12)
),
axis.line = element_line(
linewidth = 0.8
),
axis.ticks = element_line(
linewidth = 0.8
),
panel.grid.major.y = element_line(
color = "grey80",
linewidth = 0.5
),
panel.grid.minor = element_blank(),
panel.spacing = unit(
0.7,
"lines"
),
plot.margin = margin(
10,
10,
10,
10
)
)

ggsave(
filename = file.path(
figure_dir,
"Figure1E_candidate_gene_expression.pdf"
),
plot = figure_1e,
width = 6,
height = 8,
units = "in",
device = grDevices::cairo_pdf
)

ggsave(
filename = file.path(
figure_dir,
"Figure1E_candidate_gene_expression.png"
),
plot = figure_1e,
width = 6,
height = 8,
units = "in",
dpi = 300
)

message(
"Figure 1E complete: ",
file.path(
figure_dir,
"Figure1E_candidate_gene_expression.pdf"
)
)
