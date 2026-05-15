#!/usr/bin/env Rscript
# ============================================================
# Prepare BED inputs for cis-element and motif analysis.
# Positive set: SF3A2 supported novel junctions.
# Background set: NF54 supported known/canonical junctions.
# Windows: donor, acceptor_bpaware, branchpoint_focus.
# ============================================================

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
})

config_file <- Sys.getenv("SF3A2_CONFIG", unset = "config/config.yml")
if (!file.exists(config_file)) config_file <- "config/config_template.yml"
config <- yaml::read_yaml(config_file)

source("R/longread_helpers.R")

outdir <- config$paths$long_reads$motif
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

read_mapped <- function(file, strain_label) {
  dat <- readr::read_tsv(file, show_col_types = FALSE)

  required <- c(
    "isoform", "junction_number", "chrom", "strand",
    "genomic_start_coord", "genomic_end_coord",
    "junction_category", "start_site_category", "end_site_category",
    "splice_site", "pbid_fl_support", "pbid_fl_cpm", "stage"
  )

  missing <- setdiff(required, names(dat))
  if (length(missing) > 0) stop(strain_label, " missing columns: ", paste(missing, collapse = ", "))

  dat %>%
    mutate(
      strain = strain_label,
      junction_category = tolower(as.character(junction_category)),
      start_site_category = tolower(as.character(start_site_category)),
      end_site_category = tolower(as.character(end_site_category)),
      splice_site_clean = clean_splice_site(splice_site),
      canonical = splice_site_clean %in% c("GTAG", "GCAG", "ATAC", "CTAC", "CTGC", "GTAT"),
      pbid_fl_support = suppressWarnings(as.numeric(pbid_fl_support)),
      pbid_fl_cpm = suppressWarnings(as.numeric(pbid_fl_cpm)),
      novel_subtype = make_novel_subtype(junction_category, start_site_category, end_site_category),
      junction_id = paste(chrom, strand, genomic_start_coord, genomic_end_coord, sep = "_"),
      donor_pos = if_else(strand == "+", genomic_start_coord, genomic_end_coord),
      acceptor_pos = if_else(strand == "+", genomic_end_coord, genomic_start_coord)
    )
}

nf54 <- read_mapped(config$paths$long_reads$nf54_junctions_pbid, "NF54")
sf3a2 <- read_mapped(config$paths$long_reads$sf3a2_junctions_pbid, "SF3A2")

all_junc <- bind_rows(nf54, sf3a2) %>%
  filter(!is.na(pbid_fl_support), pbid_fl_support >= config$long_reads$min_pbid_fl_support) %>%
  filter(!is.na(pbid_fl_cpm), pbid_fl_cpm >= config$long_reads$min_pbid_fl_cpm)

positive <- all_junc %>%
  filter(strain == "SF3A2", junction_category == "novel")

if (config$long_reads$positive_mode != "all_novel") {
  positive <- positive %>% filter(novel_subtype == config$long_reads$positive_mode)
}

background <- all_junc %>%
  filter(strain == "NF54", junction_category == "known", canonical)

collapse_junctions <- function(df) {
  df %>%
    group_by(strain, junction_id) %>%
    summarise(
      chrom = first(chrom),
      strand = first(strand),
      donor_pos = first(donor_pos),
      acceptor_pos = first(acceptor_pos),
      junction_category = first(junction_category),
      novel_subtype = first(novel_subtype),
      max_pbid_fl_support = max(pbid_fl_support, na.rm = TRUE),
      max_pbid_fl_cpm = max(pbid_fl_cpm, na.rm = TRUE),
      stages = paste(sort(unique(stage)), collapse = ";"),
      .groups = "drop"
    )
}

positive <- collapse_junctions(positive)
background <- collapse_junctions(background)

write_tsv(positive, file.path(outdir, "SF3A2_positive_junctions.tsv"))
write_tsv(background, file.path(outdir, "NF54_background_junctions.tsv"))

write_bed <- function(df, outfile) {
  bed <- df %>%
    transmute(
      chrom,
      bed_start,
      bed_end,
      name = paste(strain, junction_id, window_type, stages, sep = "|"),
      score = 0,
      strand
    )
  readr::write_tsv(bed, outfile, col_names = FALSE)
}

windows <- config$long_reads$windows

window_specs <- list(
  donor = list(site_col = "donor_pos", rel = unlist(windows$donor)),
  acceptor_bpaware = list(site_col = "acceptor_pos", rel = unlist(windows$acceptor_bpaware)),
  branchpoint_focus = list(site_col = "acceptor_pos", rel = unlist(windows$branchpoint_focus))
)

for (w in names(window_specs)) {
  spec <- window_specs[[w]]

  pos_win <- make_site_windows(positive, spec$site_col, spec$rel[1], spec$rel[2], w)
  bg_win <- make_site_windows(background, spec$site_col, spec$rel[1], spec$rel[2], w)

  write_bed(pos_win, file.path(outdir, paste0("SF3A2_positive_", w, ".bed")))
  write_bed(bg_win, file.path(outdir, paste0("NF54_background_", w, ".bed")))
}

message("Wrote cis-element BED inputs to: ", outdir)
