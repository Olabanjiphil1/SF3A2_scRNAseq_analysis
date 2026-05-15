# ============================================================
# Asexual pseudotime heatmap
# NF54 vs SF3A2 along ER->LS trajectory.
# ============================================================

source("short_reads/scripts/00_setup.R")
load_required_packages(c("tradeSeq", "pheatmap", "gridExtra", "rlang"))

sce_asex <- readRDS(file.path(config$paths$short_reads$objects, "sce_asexual_fitGAM.rds"))

gene_lookup <- unlist(config$genes$asexual_heatmap)
genes <- intersect(names(gene_lookup), rownames(sce_asex))

df <- tradeSeq::predictSmooth(sce_asex, gene = genes, nPoints = 50, tidy = TRUE)

conditions <- unique(df$condition)
NF54_name <- get_condition_name(conditions, "NF54")
SF3A2_name <- get_condition_name(conditions, "SF3A2")

NF54_mat <- zscale_rows(make_smoother_matrix(df, NF54_name, lineage_id = 1))
SF3A2_mat <- zscale_rows(make_smoother_matrix(df, SF3A2_name, lineage_id = 1))

peak_bin <- function(m) apply(m, 1, which.max)
ord <- order(peak_bin(NF54_mat), decreasing = FALSE)

pal <- colorRampPalette(c("skyblue", "lightyellow", "red"))(100)
brks <- seq(-2, 2, length.out = 101)

png(file.path(config$paths$short_reads$figures, "pseudotime_heatmap_asexual_NF54_vs_SF3A2.png"),
    width = 3000, height = 2400, res = 300)

p1 <- pheatmap::pheatmap(NF54_mat[ord, ], cluster_rows = FALSE, cluster_cols = FALSE,
                         labels_row = make_gene_labels(NF54_mat, gene_lookup, ord),
                         show_colnames = FALSE, fontsize_row = 8,
                         color = pal, breaks = brks,
                         main = "NF54 asexual pseudotime", silent = TRUE)

p2 <- pheatmap::pheatmap(SF3A2_mat[ord, ], cluster_rows = FALSE, cluster_cols = FALSE,
                         labels_row = make_gene_labels(NF54_mat, gene_lookup, ord),
                         show_colnames = FALSE, fontsize_row = 8,
                         color = pal, breaks = brks,
                         main = "SF3A2 asexual pseudotime", silent = TRUE)

gridExtra::grid.arrange(p1$gtable, p2$gtable, ncol = 2)
dev.off()
