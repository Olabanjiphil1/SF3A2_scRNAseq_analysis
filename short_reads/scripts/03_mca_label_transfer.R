source("short_reads/scripts/00_setup.R")

obj <- readRDS(file.path(config$paths$short_reads$objects, "NF54_SF3A2_integrated_cca.rds"))

MCA_metadata <- read.csv(config$paths$short_reads$mca_metadata, row.names = 1)
MCA_counts <- read.csv(config$paths$short_reads$mca_counts, row.names = 1)

MCA <- CreateSeuratObject(counts = MCA_counts, meta.data = MCA_metadata)
DefaultAssay(MCA) <- "RNA"
MCA <- NormalizeData(MCA, verbose = FALSE)
Idents(MCA) <- "STAGE_HR2"

MCA_subset <- subset(MCA, idents = unlist(config$mca$exclude_stages), invert = TRUE)

obj <- safe_join_layers(obj)
DefaultAssay(obj) <- "RNA"

shared_features <- intersect(rownames(MCA_subset), rownames(obj))
MCA_subset <- subset(MCA_subset, features = shared_features)
obj <- subset(obj, features = shared_features)

MCA_subset <- FindVariableFeatures(MCA_subset, nfeatures = 1000, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = 1000, verbose = FALSE)
MCA_subset <- ScaleData(MCA_subset, verbose = FALSE)
MCA_subset <- RunPCA(MCA_subset, npcs = 30, features = VariableFeatures(MCA_subset), verbose = FALSE)

anchors <- FindTransferAnchors(reference = MCA_subset, query = obj,
                               dims = get_dims(config), reference.reduction = "pca")

predicted <- TransferData(anchorset = anchors, refdata = Idents(MCA_subset), dims = get_dims(config))
obj$predicted.id <- predicted$predicted.id
obj$prediction.score.max <- predicted$prediction.score.max

saveRDS(obj, file.path(config$paths$short_reads$objects, "NF54_SF3A2_MCA_annotated.rds"))

write.csv(as.data.frame.matrix(table(obj$orig.ident, obj$predicted.id)),
          file.path(config$paths$short_reads$tables, "MCA_prediction_summary.csv"))

p <- DimPlot(obj, reduction = "umap_cca", group.by = "predicted.id",
             split.by = "orig.ident", label = FALSE) +
  ggtitle("MCA-transferred stage labels")
save_plot(p, file.path(config$paths$short_reads$figures, "umap_MCA_transferred_labels.png"),
          width = 8, height = 5)
