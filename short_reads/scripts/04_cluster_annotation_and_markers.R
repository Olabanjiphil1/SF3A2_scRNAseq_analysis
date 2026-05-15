source("short_reads/scripts/00_setup.R")
load_required_packages(c("scales"))

obj <- readRDS(file.path(config$paths$short_reads$objects, "NF54_SF3A2_MCA_annotated.rds"))

Idents(obj) <- "cca_clusters"
cluster_map <- unlist(config$clusters$full_cluster_map)
obj <- RenameIdents(obj, cluster_map[levels(obj)])
obj$celltypes <- as.character(Idents(obj))

cluster_colors <- unlist(config$colors$clusters)

p <- DimPlot(obj, reduction = "umap_cca", group.by = "celltypes", split.by = "orig.ident",
             cols = cluster_colors, label = TRUE, pt.size = 0.4) + NoLegend()
save_plot(p, file.path(config$paths$short_reads$figures, "umap_annotated_clusters_split.png"),
          width = 8, height = 5)

stage_order <- unlist(config$clusters$full_stage_order)
prop_df <- obj@meta.data %>%
  mutate(cluster = factor(celltypes, levels = stage_order)) %>%
  group_by(orig.ident, cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(freq = n / sum(n))

write.csv(prop_df, file.path(config$paths$short_reads$tables, "cluster_proportions_by_sample.csv"),
          row.names = FALSE)

p_prop <- ggplot(prop_df, aes(x = orig.ident, y = freq, fill = cluster)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = cluster_colors) +
  labs(x = NULL, y = "Proportion of cells", fill = "Cluster") +
  theme_minimal(base_size = 14)
save_plot(p_prop, file.path(config$paths$short_reads$figures, "cluster_proportions_by_sample.png"),
          width = 5, height = 4)

saveRDS(obj, file.path(config$paths$short_reads$objects, "NF54_SF3A2_annotated_final.rds"))
