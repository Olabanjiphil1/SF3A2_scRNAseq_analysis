#!/usr/bin/env Rscript
# ============================================================
# Target-gene cis-context analysis for selected junctions.
# Designed for genes such as PF3D7_1438800 and PF3D7_1343000.
# Requires gene-annotated junction rows and PBID FL support columns.
# ============================================================

suppressPackageStartupMessages({
  library(yaml)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

config_file <- Sys.getenv("SF3A2_CONFIG", unset = "config/config.yml")
if (!file.exists(config_file)) config_file <- "config/config_template.yml"
config <- yaml::read_yaml(config_file)

source("R/longread_helpers.R")

outdir <- file.path(config$paths$long_reads$tables, "target_gene_cis_context")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

nf54 <- read_tsv(config$paths$long_reads$nf54_gene_annotated_junctions, show_col_types = FALSE)
sf3a2 <- read_tsv(config$paths$long_reads$sf3a2_gene_annotated_junctions, show_col_types = FALSE)

split_gene_rows <- function(df) {
  df %>%
    filter(!is.na(gene_id_from_gtf)) %>%
    mutate(gene_id_split = strsplit(gene_id_from_gtf, ";", fixed = TRUE)) %>%
    tidyr::unnest_longer(gene_id_split) %>%
    mutate(gene_id_split = trimws(gene_id_split))
}

nf54_rows <- nf54 %>% mutate(strain = "NF54") %>% split_gene_rows()
sf3a2_rows <- sf3a2 %>% mutate(strain = "SF3A2") %>% split_gene_rows()

target_genes <- unlist(config$long_reads$target_genes)

gene_rows <- bind_rows(nf54_rows, sf3a2_rows) %>%
  filter(gene_id_split %in% target_genes) %>%
  mutate(
    junction_category = tolower(as.character(junction_category)),
    start_site_category = tolower(as.character(start_site_category)),
    end_site_category = tolower(as.character(end_site_category)),
    novel_subtype = make_novel_subtype(junction_category, start_site_category, end_site_category),
    stage_value = as.character(stage),
    donor_pos = if_else(strand == "+", genomic_start_coord, genomic_end_coord),
    acceptor_pos = if_else(strand == "+", genomic_end_coord, genomic_start_coord)
  )

min_support <- config$long_reads$min_pbid_fl_support
min_cpm <- config$long_reads$min_pbid_fl_cpm

junction_summary <- gene_rows %>%
  group_by(gene_id_split, strain, junction_id) %>%
  summarise(
    junction_category = first(junction_category),
    novel_subtype = first(novel_subtype),
    chrom = first(chrom),
    strand = first(strand),
    genomic_start_coord = first(genomic_start_coord),
    genomic_end_coord = first(genomic_end_coord),
    donor_pos = first(donor_pos),
    acceptor_pos = first(acceptor_pos),
    max_pbid_fl_support = suppressWarnings(max(pbid_fl_support, na.rm = TRUE)),
    max_pbid_fl_cpm = suppressWarnings(max(pbid_fl_cpm, na.rm = TRUE)),
    n_pbids = n_distinct(isoform),
    stages = paste(sort(unique(stage_value)), collapse = ";"),
    .groups = "drop"
  ) %>%
  mutate(
    max_pbid_fl_support = ifelse(is.infinite(max_pbid_fl_support), NA, max_pbid_fl_support),
    max_pbid_fl_cpm = ifelse(is.infinite(max_pbid_fl_cpm), NA, max_pbid_fl_cpm)
  ) %>%
  filter(!is.na(max_pbid_fl_support), !is.na(max_pbid_fl_cpm)) %>%
  filter(max_pbid_fl_support >= min_support, max_pbid_fl_cpm >= min_cpm)

write_tsv(junction_summary, file.path(outdir, "target_gene_supported_junction_summary.tsv"))

# Prepare target-gene acceptor and branchpoint BED files.
windows <- list(
  acceptor_bpaware = unlist(config$long_reads$windows$acceptor_bpaware),
  branchpoint_focus = unlist(config$long_reads$windows$branchpoint_focus)
)

write_bed <- function(df, outfile) {
  bed <- df %>%
    transmute(
      chrom,
      bed_start,
      bed_end,
      name = paste(gene_id_split, strain, junction_id, window_type, stages, sep = "|"),
      score = 0,
      strand
    )
  write_tsv(bed, outfile, col_names = FALSE)
}

for (w in names(windows)) {
  win <- make_site_windows(junction_summary, "acceptor_pos", windows[[w]][1], windows[[w]][2], w)
  write_bed(win, file.path(outdir, paste0("target_gene_", w, ".bed")))
}
