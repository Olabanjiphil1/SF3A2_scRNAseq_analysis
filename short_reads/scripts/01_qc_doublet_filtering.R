source("short_reads/scripts/00_setup.R")
load_required_packages(c("DoubletFinder"))

load_10x_sample <- function(sample_name, data_dir, min_cells, min_features) {
  if (!dir.exists(data_dir)) stop("Input directory not found: ", data_dir)
  counts <- Seurat::Read10X(data.dir = data_dir)
  obj <- Seurat::CreateSeuratObject(counts = counts, project = sample_name,
                                    min.cells = min_cells, min.features = min_features)
  obj$sample <- sample_name
  obj
}

add_percent_mt <- function(obj, pattern = "MIT") {
  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = pattern)
  obj
}

filter_cells <- function(obj, min_features, max_features, max_percent_mt) {
  subset(obj, subset = nFeature_RNA > min_features &
           nFeature_RNA < max_features &
           percent.mt < max_percent_mt)
}

preprocess_seurat <- function(obj, nfeatures, dims, resolution) {
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  obj <- ScaleData(obj, features = rownames(obj), verbose = FALSE)
  obj <- RunPCA(obj, features = VariableFeatures(obj), verbose = FALSE)
  obj <- FindNeighbors(obj, dims = dims, verbose = FALSE)
  obj <- FindClusters(obj, resolution = resolution, verbose = FALSE)
  obj <- RunUMAP(obj, dims = dims, verbose = FALSE)
  obj
}

run_doubletfinder <- function(obj, dims, pN, pK, expected_doublet_rate) {
  homotypic_prop <- DoubletFinder::modelHomotypic(obj$seurat_clusters)
  n_exp <- round(expected_doublet_rate * ncol(obj))
  n_exp_adj <- round(n_exp * (1 - homotypic_prop))

  df_fun <- if ("doubletFinder" %in% getNamespaceExports("DoubletFinder")) {
    DoubletFinder::doubletFinder
  } else {
    DoubletFinder::doubletFinder_v3
  }

  obj <- df_fun(obj, PCs = dims, pN = pN, pK = pK, nExp = n_exp_adj,
                reuse.pANN = FALSE, sct = FALSE)

  df_col <- tail(grep("^DF.classifications", colnames(obj@meta.data), value = TRUE), 1)
  if (length(df_col) == 0) stop("DoubletFinder classification column was not found.")

  obj$doublet_status <- obj@meta.data[[df_col]]
  message("Doublet summary:")
  print(table(obj$doublet_status))

  subset(obj, subset = doublet_status == "Singlet")
}

sample_paths <- list(
  NF54 = config$paths$short_reads$nf54_10x,
  SF3A2 = config$paths$short_reads$sf3a2_10x
)

for (sample_name in names(sample_paths)) {
  message("Processing ", sample_name)

  obj <- load_10x_sample(
    sample_name,
    sample_paths[[sample_name]],
    config$qc$min_cells,
    config$qc$min_features
  )

  obj <- add_percent_mt(obj, config$qc$mt_pattern)

  obj <- filter_cells(
    obj,
    min_features = config$qc$min_nFeature_RNA,
    max_features = config$qc$max_nFeature_RNA,
    max_percent_mt = config$qc$max_percent_mt[[sample_name]]
  )

  obj <- preprocess_seurat(
    obj,
    nfeatures = config$seurat$variable_features_sample,
    dims = get_dims(config),
    resolution = config$seurat$resolution
  )

  obj <- run_doubletfinder(
    obj,
    dims = get_dims(config),
    pN = config$doubletfinder$pN,
    pK = config$doubletfinder$pK[[sample_name]],
    expected_doublet_rate = config$doubletfinder$expected_rate
  )

  saveRDS(obj, file.path(config$paths$short_reads$objects, paste0(sample_name, "_singlets.rds")))
}
