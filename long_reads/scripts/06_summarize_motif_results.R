#!/usr/bin/env Rscript
# ============================================================
# Summarize FASTA composition and FIMO motif enrichment.
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
summary_dir <- file.path(config$paths$long_reads$tables, "motif_summary")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

count_fasta_sequences <- function(fasta_file) {
  if (!file.exists(fasta_file)) return(0L)
  sum(startsWith(readLines(fasta_file, warn = FALSE), ">"))
}

read_fasta_tbl <- function(fasta_file, dataset_name, window_name) {
  if (!file.exists(fasta_file)) return(tibble())
  lines <- readLines(fasta_file, warn = FALSE)
  hdr_idx <- which(startsWith(lines, ">"))
  if (length(hdr_idx) == 0) return(tibble())

  seq_end_idx <- c(hdr_idx[-1] - 1L, length(lines))

  tibble(header = sub("^>", "", lines[hdr_idx]), start_idx = hdr_idx, end_idx = seq_end_idx) %>%
    rowwise() %>%
    mutate(sequence = paste(lines[(start_idx + 1L):end_idx], collapse = ""),
           dataset = dataset_name,
           window = window_name) %>%
    ungroup() %>%
    select(header, sequence, dataset, window)
}

calc_seq_metrics <- function(seq_tbl) {
  seq_tbl %>%
    mutate(
      seq_upper = toupper(sequence),
      seq_len = nchar(seq_upper),
      a_count = stringr::str_count(seq_upper, "A"),
      c_count = stringr::str_count(seq_upper, "C"),
      g_count = stringr::str_count(seq_upper, "G"),
      t_count = stringr::str_count(seq_upper, "T"),
      a_fraction = a_count / seq_len,
      pyrimidine_fraction = (c_count + t_count) / seq_len,
      gc_fraction = (g_count + c_count) / seq_len
    )
}

read_fimo_safe <- function(file, window_name, dataset_name) {
  if (!file.exists(file)) return(tibble())
  dat <- tryCatch(readr::read_tsv(file, comment = "#", show_col_types = FALSE), error = function(e) tibble())
  if (nrow(dat) == 0) return(tibble())

  names(dat) <- names(dat) %>% gsub("-", "_", .) %>% gsub(" ", "_", .) %>% tolower()
  if (!"motif_alt_id" %in% names(dat)) dat$motif_alt_id <- NA_character_
  dat %>% mutate(window = window_name, dataset = dataset_name)
}

windows <- c("donor", "acceptor_bpaware", "branchpoint_focus")

seq_tbl <- purrr::map_dfr(windows, function(w) {
  bind_rows(
    read_fasta_tbl(file.path(outdir, paste0("SF3A2_positive_", w, ".fa")), "SF3A2_positive", w),
    read_fasta_tbl(file.path(outdir, paste0("NF54_background_", w, ".fa")), "NF54_background", w)
  )
})

seq_metrics <- calc_seq_metrics(seq_tbl)
write_tsv(seq_metrics, file.path(summary_dir, "sequence_composition_per_window.tsv"))

seq_summary <- seq_metrics %>%
  group_by(window, dataset) %>%
  summarise(
    n = n(),
    median_a_fraction = median(a_fraction, na.rm = TRUE),
    median_pyrimidine_fraction = median(pyrimidine_fraction, na.rm = TRUE),
    median_gc_fraction = median(gc_fraction, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(seq_summary, file.path(summary_dir, "sequence_composition_summary.tsv"))

fimo_tbl <- purrr::map_dfr(windows, function(w) {
  bind_rows(
    read_fimo_safe(file.path(outdir, paste0("FIMO_", w, "_positive"), "fimo.tsv"), w, "SF3A2_positive"),
    read_fimo_safe(file.path(outdir, paste0("FIMO_", w, "_background"), "fimo.tsv"), w, "NF54_background")
  )
})

write_tsv(fimo_tbl, file.path(summary_dir, "all_fimo_hits.tsv"))

message("Motif summary written to: ", summary_dir)
