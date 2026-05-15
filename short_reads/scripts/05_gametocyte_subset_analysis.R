source("short_reads/scripts/00_setup.R")

obj <- readRDS(file.path(config$paths$short_reads$objects, "NF54_SF3A2_annotated_final.rds"))

Idents(obj) <- obj$celltypes
obj_gam <- subset(obj, idents = unlist(config$clusters$gametocyte_stages_short))
DefaultAssay(obj_gam) <- "RNA"
obj_gam <- safe_join_layers(obj_gam)

obj_gam <- NormalizeData(obj_gam, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
obj_gam <- FindVariableFeatures(obj_gam, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
obj_gam <- ScaleData(obj_gam, features = rownames(obj_gam), verbose = FALSE)
obj_gam <- RunPCA(obj_gam, features = VariableFeatures(obj_gam), verbose = FALSE)

gam_dims <- 1:min(15, config$seurat$dims)
obj_gam <- IntegrateLayers(obj_gam, method = CCAIntegration, orig.reduction = "pca",
                           new.reduction = "integrated.cca", verbose = FALSE)
obj_gam <- FindNeighbors(obj_gam, reduction = "integrated.cca", dims = gam_dims, verbose = FALSE)
obj_gam <- FindClusters(obj_gam, resolution = 0.2, cluster.name = "gam_cca_clusters", verbose = FALSE)
obj_gam <- RunUMAP(obj_gam, reduction = "integrated.cca", dims = gam_dims,
                   reduction.name = "umap_gam_cca", verbose = FALSE)

# Adjust this map manually if cluster IDs change after reintegration.
gam_cluster_map <- c("0" = "EG", "1" = "FG", "2" = "MG")
Idents(obj_gam) <- "gam_cca_clusters"
if (all(levels(obj_gam) %in% names(gam_cluster_map))) {
  obj_gam <- RenameIdents(obj_gam, gam_cluster_map[levels(obj_gam)])
  obj_gam$celltypes <- as.character(Idents(obj_gam))
}

saveRDS(obj_gam, file.path(config$paths$short_reads$objects, "NF54_SF3A2_gametocyte_only.rds"))

p <- DimPlot(obj_gam, reduction = "umap_gam_cca", group.by = "celltypes",
             split.by = "orig.ident", label = TRUE, pt.size = 0.4,
             cols = unlist(config$colors$clusters)) + NoLegend()
save_plot(p, file.path(config$paths$short_reads$figures, "umap_gametocyte_only_split.png"),
          width = 8, height = 5)
