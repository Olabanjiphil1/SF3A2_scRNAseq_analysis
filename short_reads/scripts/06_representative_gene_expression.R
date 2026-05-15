source("short_reads/scripts/00_setup.R")

obj_gam <- readRDS(file.path(config$paths$short_reads$objects, "NF54_SF3A2_gametocyte_only.rds"))

genes <- unlist(config$genes$representative)
gene_ids <- intersect(names(genes), rownames(obj_gam))
genes <- genes[gene_ids]

expr_data <- FetchData(obj_gam, vars = c(gene_ids, "celltypes", "orig.ident"))

expr_long <- expr_data %>%
  pivot_longer(cols = all_of(gene_ids), names_to = "GeneID", values_to = "Expression") %>%
  mutate(GeneName = genes[GeneID], Genotype = orig.ident) %>%
  group_by(GeneID) %>%
  mutate(threshold = 0.25 * max(Expression, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(Expression >= threshold, celltypes %in% unlist(config$clusters$gametocyte_stages_short)) %>%
  mutate(celltypes = factor(celltypes, levels = unlist(config$clusters$gametocyte_stages_short)))

sig_df <- expr_long %>%
  group_by(GeneName, celltypes) %>%
  summarise(
    p_val = tryCatch(wilcox.test(Expression ~ Genotype)$p.value, error = function(e) NA_real_),
    y_pos = max(Expression, na.rm = TRUE) * 1.1,
    .groups = "drop"
  ) %>%
  mutate(signif_label = pvalue_to_stars(p_val))

write.csv(expr_long, file.path(config$paths$short_reads$tables, "representative_gene_expression_filtered.csv"), row.names = FALSE)
write.csv(sig_df, file.path(config$paths$short_reads$tables, "representative_gene_expression_wilcoxon.csv"), row.names = FALSE)

p <- ggplot(expr_long, aes(x = Genotype, y = Expression, fill = Genotype)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, color = "black") +
  facet_grid(GeneName ~ celltypes, scales = "free_y") +
  geom_segment(data = sig_df, aes(x = 1, xend = 2, y = y_pos, yend = y_pos),
               inherit.aes = FALSE, linewidth = 0.6) +
  geom_text(data = sig_df, aes(x = 1.5, y = y_pos * 1.03, label = signif_label),
            inherit.aes = FALSE, size = 5, fontface = "bold") +
  scale_fill_manual(values = unlist(config$colors$genotype)) +
  labs(title = "Representative gene expression >=25% max",
       x = "Genotype", y = "Expression (log-normalized)") +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1))

save_plot(p, file.path(config$paths$short_reads$figures, "representative_gene_expression_boxplots.png"),
          width = 9, height = 7)
