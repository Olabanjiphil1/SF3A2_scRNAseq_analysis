# ============================================================
# Asexual-stage Slingshot + tradeSeq workflow
# Stages: ER, LR, ET, LT, ES, LS
# ============================================================

source("short_reads/scripts/00_setup.R")
load_required_packages(c("SingleCellExperiment", "slingshot", "tradeSeq", "SummarizedExperiment"))

obj <- readRDS(file.path(config$paths$short_reads$objects, "NF54_SF3A2_annotated_final.rds"))

asexual_stages <- unlist(config$clusters$asexual_stages)
obj_asex <- subset(obj, subset = celltypes %in% asexual_stages)

DefaultAssay(obj_asex) <- "RNA"
obj_asex <- safe_join_layers(obj_asex)

sce_asex <- Seurat::as.SingleCellExperiment(obj_asex, assay = "RNA")

# Use integrated UMAP from Seurat object as Slingshot reducedDim.
if ("umap_cca" %in% names(obj_asex@reductions)) {
  SingleCellExperiment::reducedDims(sce_asex)$UMAP <- obj_asex@reductions$umap_cca@cell.embeddings
}

sce_asex <- slingshot::slingshot(
  sce_asex,
  reducedDim = "UMAP",
  clusterLabels = "celltypes",
  start.clus = "ER",
  end.clus = "LS",
  approx_points = 150
)

saveRDS(sce_asex, file.path(config$paths$short_reads$objects, "sce_asexual_slingshot.rds"))

condition <- factor(SingleCellExperiment::colData(sce_asex)$orig.ident)

sce_asex <- tradeSeq::fitGAM(
  sce_asex,
  conditions = condition,
  nknots = 6,
  verbose = TRUE
)

saveRDS(sce_asex, file.path(config$paths$short_reads$objects, "sce_asexual_fitGAM.rds"))

cond_res <- tradeSeq::conditionTest(sce_asex, l2fc = log2(2)) %>%
  as.data.frame()

cond_res$padj <- p.adjust(cond_res$pvalue, method = "fdr")
cond_res$gene <- rownames(cond_res)
cond_res <- cond_res[, c("gene", setdiff(colnames(cond_res), "gene"))]

write.csv(cond_res,
          file.path(config$paths$short_reads$tables, "tradeseq_conditionTest_asexual_all_genes.csv"),
          row.names = FALSE)

sig_res <- cond_res %>%
  filter(!is.na(padj), padj <= 0.05)

write.csv(sig_res,
          file.path(config$paths$short_reads$tables, "tradeseq_conditionTest_asexual_significant_genes.csv"),
          row.names = FALSE)

message("Asexual significant pseudotime genes: ", nrow(sig_res))
