suppressPackageStartupMessages({
    library(chiptsne2)
})

ct2 <- readRDS("Figure1G_ct2.rds")

cat("Loaded ct2\n")

ht <- plotSignalHeatmap(
    ct2,
    group_VARS = "group",
    sort_strategy = "none",
    relative_heatmap_width = 0.4,
    heatmap_fill_limits = c(0,150),
    heatmap_colors = c(
        "#577AB2",
        "#FFFFE0",
        "#BC412B"
    )
)

cat("Heatmap OK\n")

pdf("test_heatmap.pdf", width=7, height=8)
print(ht)
dev.off()

cat("PDF written\n")
