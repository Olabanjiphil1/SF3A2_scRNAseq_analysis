#!/usr/bin/env Rscript
# ============================================================
# Strict SQANTI3 junction comparison using PBID FL support.
# Novel = at least one unannotated splice site.
# Excludes likely artifacts and selected structural categories.
# ============================================================

suppressPackageStartupMessages({
  library(yaml)
  library(tidyverse)
  library(readxl)
  library(janitor)
  library(writexl)
  library(scales)
})

config_file <- Sys.getenv("SF3A2_CONFIG", unset = "config/config.yml")
if (!file.exists(config_file)) config_file <- "config/config_template.yml"
config <- yaml::read_yaml(config_file)

source("R/longread_helpers.R")

outdir <- file.path(config$paths$long_reads$tables, "junction_comparison_strict_novelty")
figdir <- file.path(config$paths$long_reads$figures, "junction_comparison_strict_novelty")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

read_mapped_junction_xlsx <- function(xlsx_file, strain_label) {
  dat <- readxl::read_excel(xlsx_file, sheet = "junction_rows") %>%
    janitor::clean_names()

  required <- c(
    "isoform", "junction_number", "chrom", "strand",
    "genomic_start_coord", "genomic_end_coord",
    "junction_category", "start_site_category", "end_site_category",
    "bite_junction", "splice_site", "rts_junction", "indel_near_junct",
    "pbid_fl_support", "pbid_fl_cpm"
  )

  missing <- setdiff(required, names(dat))
  if (length(missing) > 0) {
    stop(strain_label, " missing required columns: ", paste(missing, collapse = ", "))
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
      pbid_fl_support = suppressWarnings(as.numeric(pbid_fl_support)),
      pbid_fl_cpm = suppressWarnings(as.numeric(pbid_fl_cpm)),
      stage = if ("stage" %in% names(.)) trimws(as.character(stage)) else NA_character_,
      structural_category = if ("structural_category" %in% names(.)) clean_structural_category(structural_category) else NA_character_,
      splice_site_clean = clean_splice_site(splice_site),
      splice_class = case_when(
        splice_site_clean %in% c("GTAG", "GCAG", "ATAC", "CTAC", "CTGC", "GTAT") ~ "Canonical_or_reverse-complement",
        is.na(splice_site_clean) | splice_site_clean == "" ~ NA_character_,
        TRUE ~ "Non-canonical"
      ),
      junction_id = paste(chrom, strand, genomic_start_coord, genomic_end_coord, sep = "_"),
      novelty_strict = strict_novelty_call(junction_category, start_site_category, end_site_category)
    )
}

nf54 <- read_mapped_junction_xlsx(config$paths$long_reads$nf54_junctions_pbid_xlsx, "NF54")
sf3a2 <- read_mapped_junction_xlsx(config$paths$long_reads$sf3a2_junctions_pbid_xlsx, "SF3A2")

excluded_structural <- unlist(config$long_reads$excluded_structural_categories)

junc <- bind_rows(nf54, sf3a2) %>%
  filter(!is.na(stage), stage != "") %>%
  filter(!is.na(pbid_fl_support), pbid_fl_support >= config$long_reads$strict_min_pbid_fl_support) %>%
  filter(is.na(rts_junction) | rts_junction == FALSE) %>%
  filter(is.na(indel_near_junct) | indel_near_junct == FALSE) %>%
  filter(is.na(bite_junction) | bite_junction == FALSE) %>%
  filter(!is.na(novelty_strict)) %>%
  filter(is.na(structural_category) | !structural_category %in% excluded_structural)

write_tsv(junc, file.path(outdir, "row_level_strictQC_stageFiltered.tsv"))

junc_unique <- junc %>%
  group_by(strain, junction_id) %>%
  summarise(
    chrom = first_non_na(chrom),
    strand = first_non_na(strand),
    genomic_start_coord = first_non_na(genomic_start_coord),
    genomic_end_coord = first_non_na(genomic_end_coord),
    novelty_strict = first_non_na(novelty_strict),
    splice_class = first_non_na(splice_class),
    max_pbid_fl_support = max(pbid_fl_support, na.rm = TRUE),
    max_pbid_fl_cpm = max(pbid_fl_cpm, na.rm = TRUE),
    stages = paste(sort(unique(stage)), collapse = ";"),
    n_isoforms = n_distinct(isoform),
    .groups = "drop"
  )

write_tsv(junc_unique, file.path(outdir, "unique_junctions_strictQC.tsv"))
writexl::write_xlsx(list(row_level = junc, unique_junctions = junc_unique),
                    file.path(outdir, "junction_comparison_strictQC.xlsx"))

summary_tbl <- junc_unique %>%
  count(strain, novelty_strict) %>%
  group_by(strain) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

write_tsv(summary_tbl, file.path(outdir, "summary_strict_known_vs_novel.tsv"))

p <- ggplot(summary_tbl, aes(x = strain, y = prop, fill = novelty_strict)) +
  geom_col(width = 0.6) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Proportion of unique supported junctions", fill = "Strict novelty") +
  theme_pub()

save_both(p, figdir, "strict_supported_known_vs_novel", width = 5, height = 4)
