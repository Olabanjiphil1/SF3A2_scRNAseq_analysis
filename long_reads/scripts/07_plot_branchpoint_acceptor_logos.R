#!/usr/bin/env Rscript
# ============================================================
# Plot NF54 vs SF3A2 sequence logos for branchpoint/acceptor windows.
# Finds the most divergent local block before plotting.
# ============================================================

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(ggplot2)
  library(ggseqlogo)
  library(patchwork)
})

config_file <- Sys.getenv("SF3A2_CONFIG", unset = "config/config.yml")
if (!file.exists(config_file)) config_file <- "config/config_template.yml"
config <- yaml::read_yaml(config_file)

source("R/longread_helpers.R")

motif_dir <- config$paths$long_reads$motif
outdir <- config$paths$long_reads$figures
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

plot_window_pair <- function(nf54_fa, sf3a2_fa, window_title, k = 9) {
  nf54_seqs <- read_fasta_sequences(nf54_fa) |> keep_modal_length()
  sf3a2_seqs <- read_fasta_sequences(sf3a2_fa) |> keep_modal_length()

  common_len <- min(as.integer(names(table(nchar(nf54_seqs)))),
                    as.integer(names(table(nchar(sf3a2_seqs)))))
  nf54_seqs <- nf54_seqs[nchar(nf54_seqs) == common_len]
  sf3a2_seqs <- sf3a2_seqs[nchar(sf3a2_seqs) == common_len]

  block <- find_most_divergent_block(nf54_seqs, sf3a2_seqs, k = k)

  nf54_block <- substr(nf54_seqs, block$start, block$end)
  sf3a2_block <- substr(sf3a2_seqs, block$start, block$end)

  p_nf54 <- ggseqlogo::ggseqlogo(nf54_block, method = "prob") +
    labs(title = "NF54 matched background",
         subtitle = paste(window_title, "| positions", block$start, "-", block$end),
         x = NULL, y = "Base probability") +
    theme_classic(base_size = 13)

  p_sf3a2 <- ggseqlogo::ggseqlogo(sf3a2_block, method = "prob") +
    labs(title = "SF3A2 positive junctions",
         subtitle = paste(window_title, "| positions", block$start, "-", block$end),
         x = NULL, y = "Base probability") +
    theme_classic(base_size = 13)

  p_nf54 / p_sf3a2
}

windows <- c("acceptor_bpaware", "branchpoint_focus")

for (w in windows) {
  p <- plot_window_pair(
    nf54_fa = file.path(motif_dir, paste0("NF54_background_", w, ".fa")),
    sf3a2_fa = file.path(motif_dir, paste0("SF3A2_positive_", w, ".fa")),
    window_title = w,
    k = 9
  )

  ggsave(file.path(outdir, paste0("logo_", w, "_NF54_vs_SF3A2.png")),
         p, width = 8, height = 5, dpi = 600, bg = "white")
  ggsave(file.path(outdir, paste0("logo_", w, "_NF54_vs_SF3A2.pdf")),
         p, width = 8, height = 5, bg = "white")
}
