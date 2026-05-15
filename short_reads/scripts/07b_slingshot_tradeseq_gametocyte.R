# ============================================================
# Sexual/gametocyte Slingshot + tradeSeq workflow
# Stages: SP, EF, EM
# This is intentionally separate from asexual pseudotime.
# Expected input object can be an induced/sexual integrated object with:
#   - celltypes = SP / EF / EM
#   - orig.ident = NF54 / SF3A2
#   - reduction = umap.harmony or umap_cca
# ============================================================

source("short_reads/scripts/00_setup.R")
load_required_packages(c("SingleCellExperiment", "slingshot", "tradeSeq", "SummarizedExperiment", "scales"))

# Prefer a dedicated induced-sexual object if present.
sexual_file <- file.path(config$paths$short_reads$objects, "integrated_induced_NF54_SF3A2.rds")
fallback_file <- file.path(config$paths$short_reads$objects, "NF54_SF3A2_gametocyte_only.rds")

if (file.exists(sexual_file)) {
  comm <- readRDS(sexual_file)
} else {
  warning("Dedicated induced-sexual object not found; using gametocyte-only fallback object.")
  comm <- readRDS(fallback_file)
}

sexual_stages <- unlist(config$clusters$gametocyte_stages_induced)
comm_gams <- subset(comm, subset = celltypes %in% sexual_stages)

DefaultAssay(comm_gams) <- "RNA"
comm_gams <- safe_join_layers(comm_gams)

sce_gam <- Seurat::as.SingleCellExperiment(comm_gams, assay = "RNA")

# Prefer Harmony UMAP if available because the original sexual workflow used umap.harmony.
if ("umap.harmony" %in% names(comm_gams@reductions)) {
  SingleCellExperiment::reducedDims(sce_gam)$UMAP <- Seurat::Embeddings(comm_gams, "umap.harmony")
} else if ("umap_cca" %in% names(comm_gams@reductions)) {
  SingleCellExperiment::reducedDims(sce_gam)$UMAP <- Seurat::Embeddings(comm_gams, "umap_cca")
} else if ("umap" %in% names(comm_gams@reductions)) {
  SingleCellExperiment::reducedDims(sce_gam)$UMAP <- Seurat::Embeddings(comm_gams, "umap")
} else {
  stop("No UMAP reduction found for gametocyte Slingshot.")
}

sce_gam <- slingshot::slingshot(
  sce_gam,
  reducedDim = "UMAP",
  clusterLabels = "celltypes",
  start.clus = "SP",
  approx_points = 150
)

saveRDS(sce_gam, file.path(config$paths$short_reads$objects, "sce_gametocyte_slingshot.rds"))

# Optional trajectory plot
umap_coords <- SingleCellExperiment::reducedDims(sce_gam)$UMAP
pt_mat <- slingshot::slingPseudotime(sce_gam)
combined_pt <- do.call(pmax, c(as.data.frame(pt_mat), list(na.rm = TRUE)))

png(file.path(config$paths$short_reads$figures, "gametocyte_umap_with_slingshot_trajectory.png"),
    width = 5, height = 4, units = "in", res = 300)
palette <- grDevices::colorRampPalette(c("navy", "yellow"))(100)
breaks <- seq(min(combined_pt, na.rm = TRUE), max(combined_pt, na.rm = TRUE), length.out = 101)
bins <- findInterval(combined_pt, vec = breaks, all.inside = TRUE)
cols <- scales::alpha(palette[bins], 0.8)
cols[is.na(combined_pt)] <- "grey80"
plot(umap_coords, col = cols, pch = 16, asp = 1,
     xlab = "UMAP1", ylab = "UMAP2", main = "Gametocyte Slingshot pseudotime")
lines(slingshot::SlingshotDataSet(sce_gam), lwd = 2)
dev.off()

condition <- factor(SingleCellExperiment::colData(sce_gam)$orig.ident)

sce_gam <- tradeSeq::fitGAM(
  sce_gam,
  conditions = condition,
  nknots = 6,
  verbose = TRUE
)

saveRDS(sce_gam, file.path(config$paths$short_reads$objects, "sce_gametocyte_fitGAM.rds"))

cond_res <- tradeSeq::conditionTest(sce_gam, l2fc = log2(2)) %>%
  as.data.frame()

cond_res$padj <- p.adjust(cond_res$pvalue, method = "fdr")
cond_res$gene <- rownames(cond_res)
cond_res <- cond_res[, c("gene", setdiff(colnames(cond_res), "gene"))]

write.csv(cond_res,
          file.path(config$paths$short_reads$tables, "tradeseq_conditionTest_gametocyte_all_genes.csv"),
          row.names = FALSE)

sig_res <- cond_res %>%
  filter(!is.na(padj), padj <= 0.05)

write.csv(sig_res,
          file.path(config$paths$short_reads$tables, "tradeseq_conditionTest_gametocyte_significant_genes.csv"),
          row.names = FALSE)

message("Gametocyte significant pseudotime genes: ", nrow(sig_res))
