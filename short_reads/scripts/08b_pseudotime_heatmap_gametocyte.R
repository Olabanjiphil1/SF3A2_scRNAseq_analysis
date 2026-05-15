# ============================================================
# Gametocyte pseudotime heatmap
# 2x2 layout:
#   NF54 SP->EM
#   SF3A2 SP->EM
#   NF54 SP->EF
#   SF3A2 SP->EF
# ============================================================

source("short_reads/scripts/00_setup.R")
load_required_packages(c("tradeSeq", "pheatmap", "gridExtra", "rlang"))

sce_gam <- readRDS(file.path(config$paths$short_reads$objects, "sce_gametocyte_fitGAM.rds"))

gene_lookup <- unlist(config$genes$gametocyte_heatmap)

sig_file <- file.path(config$paths$short_reads$tables, "tradeseq_conditionTest_gametocyte_significant_genes.csv")
if (file.exists(sig_file)) {
  genes <- read.csv(sig_file)$gene
  genes <- intersect(genes, rownames(sce_gam))
  if (length(genes) == 0) genes <- intersect(names(gene_lookup), rownames(sce_gam))
} else {
  genes <- intersect(names(gene_lookup), rownames(sce_gam))
}

df <- tradeSeq::predictSmooth(sce_gam, gene = genes, nPoints = 50, tidy = TRUE)

conditions <- unique(df$condition)
NF54_name <- get_condition_name(conditions, "NF54")
SF3A2_name <- get_condition_name(conditions, "SF3A2")

NF54_L1 <- make_smoother_matrix(df, NF54_name, lineage_id = 1)
SF3A2_L1 <- make_smoother_matrix(df, SF3A2_name, lineage_id = 1)
NF54_L2 <- make_smoother_matrix(df, NF54_name, lineage_id = 2)
SF3A2_L2 <- make_smoother_matrix(df, SF3A2_name, lineage_id = 2)

mat_L1 <- cbind(zscale_rows(NF54_L1), zscale_rows(SF3A2_L1))
mat_L2 <- cbind(zscale_rows(NF54_L2), zscale_rows(SF3A2_L2))

peak_bin <- function(m) apply(m, 1, which.max)
ord_L1 <- order(peak_bin(mat_L1[, 1:50, drop = FALSE]), decreasing = TRUE)
ord_L2 <- order(peak_bin(mat_L2[, 1:50, drop = FALSE]), decreasing = TRUE)

pal <- colorRampPalette(c("skyblue", "lightyellow", "red"))(100)
brks <- seq(-2, 4, length.out = 101)

plot_ht <- function(mat, ord, cols, title) {
  pheatmap::pheatmap(
    mat[ord, cols, drop = FALSE],
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    labels_row = make_gene_labels(mat, gene_lookup, ord),
    show_colnames = FALSE,
    fontsize_row = 5,
    color = pal,
    breaks = brks,
    main = title,
    silent = TRUE
  )$gtable
}

png(file.path(config$paths$short_reads$figures, "pseudotime_heatmap_gametocyte_lineages.png"),
    width = 3000, height = 3000, res = 300)

gridExtra::grid.arrange(
  plot_ht(mat_L1, ord_L1, 1:50, "NF54 SP→EM"),
  plot_ht(mat_L1, ord_L1, 51:100, "SF3A2 SP→EM"),
  plot_ht(mat_L2, ord_L2, 1:50, "NF54 SP→EF"),
  plot_ht(mat_L2, ord_L2, 51:100, "SF3A2 SP→EF"),
  ncol = 2
)
dev.off()
