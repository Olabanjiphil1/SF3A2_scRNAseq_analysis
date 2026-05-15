source("short_reads/scripts/00_setup.R")

NF54 <- readRDS(file.path(config$paths$short_reads$objects, "NF54_singlets.rds"))
SF3A2 <- readRDS(file.path(config$paths$short_reads$objects, "SF3A2_singlets.rds"))

obj <- merge(NF54, SF3A2)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst",
                            nfeatures = config$seurat$variable_features_integrated,
                            verbose = FALSE)
obj <- ScaleData(obj, features = rownames(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), verbose = FALSE)

obj <- IntegrateLayers(obj, method = CCAIntegration, orig.reduction = "pca",
                       new.reduction = "integrated.cca", verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = get_dims(config), verbose = FALSE)
obj <- FindClusters(obj, resolution = config$seurat$resolution,
                    cluster.name = "cca_clusters", verbose = FALSE)
obj <- RunUMAP(obj, reduction = "integrated.cca", dims = get_dims(config),
               reduction.name = "umap_cca", verbose = FALSE)

saveRDS(obj, file.path(config$paths$short_reads$objects, "NF54_SF3A2_integrated_cca.rds"))

p <- DimPlot(obj, reduction = "umap_cca", group.by = "orig.ident") +
  ggtitle("NF54 and SF3A2 CCA-integrated UMAP")
save_plot(p, file.path(config$paths$short_reads$figures, "umap_integrated_by_sample.png"))
