#!/usr/bin/env Rscript
# ============================================================
# SQANTI3 junctions.txt comparison: NF54 vs SF3A2
# Publication-focused junction-level summary.
# ============================================================

suppressPackageStartupMessages({
  library(yaml)
  library(tidyverse)
  library(janitor)
  library(scales)
  library(patchwork)
})

config_file <- Sys.getenv("SF3A2_CONFIG", unset = "config/config.yml")
if (!file.exists(config_file)) config_file <- "config/config_template.yml"
config <- yaml::read_yaml(config_file)

source("R/longread_helpers.R")

outdir <- file.path(config$paths$long_reads$tables, "junction_comparison")
figdir <- file.path(config$paths$long_reads$figures, "junction_comparison")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

read_sqanti_junctions <- function(file, strain_label) {
  dat <- readr::read_tsv(file, show_col_types = FALSE, na = c("", "NA", ".", "null")) %>%
    janitor::clean_names()

  required <- c(
    "isoform", "junction_number", "chrom", "strand",
    "genomic_start_coord", "genomic_end_coord",
    "junction_category", "start_site_category", "end_site_category",
    "diff_to_ref_start_site", "diff_to_ref_end_site",
    "bite_junction", "splice_site", "rts_junction", "indel_near_junct"
  )

  missing <- setdiff(required, names(dat))
  if (length(missing) > 0) {
    stop(strain_label, " missing columns: ", paste(missing, collapse = ", "))
  }

  dat %>%
    mutate(
      strain = strain_label,
      junction_category = tolower(as.character(junction_category)),
      start_site_category = tolower(as.character(start_site_category)),
      end_site_category = tolower(as.character(end_site_category)),
      bite_junction = to_logical(bite_junction),
      rts_junction = to_logical(rts_junction),
      indel_near_junct = to_logical(indel_near_junct),
      diff_to_ref_start_site = suppressWarnings(as.numeric(diff_to_ref_start_site)),
      diff_to_ref_end_site = suppressWarnings(as.numeric(diff_to_ref_end_site)),
      abs_diff_to_ref_start_site = abs(diff_to_ref_start_site),
      abs_diff_to_ref_end_site = abs(diff_to_ref_end_site),
      splice_site_clean = clean_splice_site(splice_site),
      splice_class = case_when(
        splice_site_clean %in% c("GTAG", "GCAG", "ATAC") ~ "Canonical",
        is.na(splice_site_clean) | splice_site_clean == "" ~ NA_character_,
        TRUE ~ "Non-canonical"
      ),
      junction_id = paste(chrom, strand, genomic_start_coord, genomic_end_coord, sep = "_"),
      site_combo = case_when(
        start_site_category == "known" & end_site_category == "known" ~ "Both sites annotated",
        start_site_category == "known" & end_site_category != "known" ~ "Start annotated only",
        start_site_category != "known" & end_site_category == "known" ~ "End annotated only",
        TRUE ~ "Neither site annotated"
      )
    )
}

nf54 <- read_sqanti_junctions(config$paths$long_reads$nf54_junctions_txt, "NF54")
sf3a2 <- read_sqanti_junctions(config$paths$long_reads$sf3a2_junctions_txt, "SF3A2")
junc <- bind_rows(nf54, sf3a2)

write_tsv(junc, file.path(outdir, "all_junction_rows_cleaned.tsv"))

junc_unique <- junc %>%
  group_by(strain, junction_id) %>%
  summarise(
    chrom = first_non_na(chrom),
    strand = first_non_na(strand),
    genomic_start_coord = first_non_na(genomic_start_coord),
    genomic_end_coord = first_non_na(genomic_end_coord),
    junction_category = first_non_na(junction_category),
    start_site_category = first_non_na(start_site_category),
    end_site_category = first_non_na(end_site_category),
    site_combo = first_non_na(site_combo),
    splice_site_clean = first_non_na(splice_site_clean),
    splice_class = first_non_na(splice_class),
    bite_junction = any(bite_junction %in% TRUE, na.rm = TRUE),
    rts_junction = any(rts_junction %in% TRUE, na.rm = TRUE),
    indel_near_junct = any(indel_near_junct %in% TRUE, na.rm = TRUE),
    abs_diff_to_ref_start_site = first_non_na(abs_diff_to_ref_start_site),
    abs_diff_to_ref_end_site = first_non_na(abs_diff_to_ref_end_site),
    n_isoforms_supporting = n_distinct(isoform),
    .groups = "drop"
  ) %>%
  mutate(min_shift_to_ref = safe_min_shift(abs_diff_to_ref_start_site, abs_diff_to_ref_end_site))

junction_occurrence <- junc_unique %>%
  distinct(strain, junction_id) %>%
  count(junction_id, name = "n_strains")

junc_unique <- junc_unique %>%
  left_join(junction_occurrence, by = "junction_id") %>%
  mutate(
    sharing_status = if_else(n_strains == 2, "Shared", "Strain-specific"),
    owner = if_else(n_strains == 2, "Shared", strain)
  )

write_tsv(junc_unique, file.path(outdir, "unique_junctions_cleaned.tsv"))

summary_unique <- junc_unique %>%
  count(strain, junction_category) %>%
  group_by(strain) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

write_tsv(summary_unique, file.path(outdir, "summary_unique_known_vs_novel.tsv"))

p <- ggplot(summary_unique, aes(x = strain, y = prop, fill = junction_category)) +
  geom_col(width = 0.6) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Proportion of unique junctions", fill = "Junction category") +
  theme_pub()

save_both(p, figdir, "unique_known_vs_novel_junctions", width = 5, height = 4)
