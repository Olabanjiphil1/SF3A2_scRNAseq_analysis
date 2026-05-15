
#!/usr/bin/env Rscript
## ============================================================
## 09_sf3a2_strict_novel_sp_ef_em_summary.R
##
## Purpose:
##   Final long-read target-junction summary module for the SF3A2 project.
##
##   This script summarizes SF3A2-specific strict-novel splice junctions
##   across SP, EF, and EM gametocyte stages. It extracts donor,
##   acceptor, branchpoint-proximal, and matched reference sequence windows,
##   assigns target-gene context, compares REF vs SF3A2 junction sequences,
##   and predicts local coding-frame effects.
##
## Notes for GitHub use:
##   - The original analysis logic is preserved below.
##   - Before running, edit the path variables in STEP 1 or adapt them to
##     config/config.yml.
##   - Required inputs include:
##       1. PlasmoDB GFF/GTF annotation
##       2. PlasmoDB genome FASTA
##       3. SF3A2 PBID-supported junction table
##       4. NF54 PBID-supported junction table
##
## Suggested run:
##   Rscript long_reads/scripts/09_sf3a2_strict_novel_sp_ef_em_summary.R
## ============================================================


## ============================================================
## SF3A2-specific strict-novel splice-junction summary
## SP / EF / EM stages
##
## No Rsamtools version
##
## Summarizes:
## - SF3A2-specific strict-novel junctions across SP, EF, EM
## - junction identity
## - long-read PBID support
## - splice-site novelty
## - donor/acceptor/branchpoint-proximal sequence windows
## - predicted coding-frame effects
## - nearest reference exon-exon context
## - direct REF vs SF3A2 junction sequence comparison
## ============================================================


## ============================================================
## STEP 0. Load packages
## ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(Biostrings)
  library(GenomicRanges)
  library(rtracklayer)
})


## ============================================================
## STEP 1. Paths and analysis settings
## ============================================================

ref_gtf <- "/Users/olatunbosunaringbangba/Documents/Research/single_cell_RNASeq/short_reads/analysis/final/longreads/SF3A2/stage_bams_SF3A2/PlasmoDB-68_Pfalciparum3D7 3.gff"

genome_FASTA <- "/Users/olatunbosunaringbangba/Documents/Research/single_cell_RNASeq/short_reads/analysis/final/longreads/SF3A2/stage_bams_SF3A2/PlasmoDB-68_Pfalciparum3D7_Genome.fasta"

sf3a2_file <- "/Users/olatunbosunaringbangba/Documents/Research/single_cell_RNASeq/short_reads/analysis/final/longreads/junction_FL_per_PBID-B/SF3A2_junction_rows_with_PBID_FL_support_noArtifacts_noNAstage.tsv"

nf54_file <- "/Users/olatunbosunaringbangba/Documents/Research/single_cell_RNASeq/short_reads/analysis/final/longreads/junction_FL_per_PBID-B/NF54_junction_rows_with_PBID_FL_support_noArtifacts_noNAstage.tsv"

out_dir <- "/Users/olatunbosunaringbangba/Documents/Research/single_cell_RNASeq/short_reads/analysis/final/longreads/SF3A2_strict_novel_target_junction_sequences"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

target_stages <- c("SP", "EF", "EM")

## Strict-novel definition:
## strict novel = junction_category == "novel"
## AND at least one splice site is novel.
strict_requires_junction_category_novel <- TRUE
strict_requires_at_least_one_novel_site <- TRUE

## Optional artifact filters if columns exist.
exclude_rts_junction <- TRUE
exclude_indel_near_junction <- TRUE
exclude_bite_junction <- TRUE

## Junction-file coordinate convention.
## In your SQANTI junction file, donor_site/genomic_start_coord often behaves like
## the first intronic base on the + strand.
## Therefore, convert observed donor to the true exon boundary.
observed_donor_coord_type <- "first_intronic_base"  # options: "first_intronic_base", "exon_boundary"

## Acceptor is treated as the first exonic base.
observed_acceptor_coord_type <- "first_exonic_base" # options: "first_exonic_base", "last_intronic_base"

## Sequence window settings
donor_exonic_nt <- 5
donor_intronic_nt <- 20

acceptor_intronic_nt <- 60
acceptor_exonic_nt <- 10

acceptor_core_upstream_nt <- 20
acceptor_core_exonic_nt <- 5

bp_upstream_from <- 60
bp_upstream_to <- 15

## Direct REF-vs-SF3A2 junction comparison.
## REF = last 15 nt upstream exon + first 15 nt downstream exon.
## SF3A2 = same reference neighborhood but using mutant donor/acceptor.
ref_junction_nt_each_side <- 15

target_genes <- c(
  "PF3D7_1438800",
  "PF3D7_1312800",
  "PF3D7_1466500",
  "PF3D7_1441300",
  "PF3D7_0114000",
  "PF3D7_1327100",
  "PF3D7_1122200",
  "PF3D7_1401900",
  "PF3D7_1120600",
  "PF3D7_1221300",
  "PF3D7_0205300",
  "PF3D7_1430200",
  "PF3D7_0709400",
  "PF3D7_0109000",
  "PF3D7_0311700",
  "PF3D7_1025500",
  "PF3D7_1025000",
  "PF3D7_1322900",
  "PF3D7_0306200",
  "PF3D7_0205100",
  "PF3D7_0216700",
  "PF3D7_1016300",
  "PF3D7_1202100",
  "PF3D7_1368400",
  "PF3D7_0523700",
  "PF3D7_0831800",
  "PF3D7_1221100",
  "PF3D7_1446100",
  "PF3D7_0922600",
  "PF3D7_0916200",
  "PF3D7_1020200",
  "PF3D7_1343000",
  "PF3D7_1361300",
  "PF3D7_1240600",
  "PF3D7_0922300",
  "PF3D7_1415100",
  "PF3D7_1347600",
  "PF3D7_0605600",
  "PF3D7_1369300"
)


## ============================================================
## STEP 2. General helper functions
## ============================================================

extract_pf_gene_id <- function(x) {
  x <- as.character(x)

  stringr::str_extract(x, "PF3D7[_-][0-9]+") |>
    stringr::str_replace("-", "_")
}

as_bool <- function(x) {
  if (is.logical(x)) {
    return(x)
  }

  x <- tolower(as.character(x))

  x %in% c("true", "t", "1", "yes", "y")
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

safe_int <- function(x) {
  suppressWarnings(as.integer(x))
}

clean_dna <- function(x) {
  x <- toupper(as.character(x))

  if (length(x) == 0 || is.na(x[1]) || x[1] == "") {
    return(NA_character_)
  }

  gsub("[^ACGT]", "", x)
}

base_fraction <- function(seq, bases = c("A")) {
  seq <- clean_dna(seq)

  if (is.na(seq) || nchar(seq) == 0) {
    return(NA_real_)
  }

  chars <- strsplit(seq, "")[[1]]

  mean(chars %in% bases)
}

translate_codon <- function(codon) {
  codon <- clean_dna(codon)

  if (is.na(codon) || nchar(codon) != 3) {
    return(" ")
  }

  as.character(
    Biostrings::translate(
      Biostrings::DNAString(codon),
      if.fuzzy.codon = "X"
    )
  )
}

split_into_codons <- function(seq) {
  seq <- clean_dna(seq)

  if (is.na(seq) || nchar(seq) < 3) {
    return(character())
  }

  trim_len <- nchar(seq) - (nchar(seq) %% 3)
  seq <- substr(seq, 1, trim_len)

  starts <- seq(1, nchar(seq), by = 3)

  substr(seq, starts, starts + 2)
}

format_dna_codons_with_boundary <- function(seq, boundary_after_nt) {
  seq <- clean_dna(seq)

  if (is.na(seq) || nchar(seq) < 3) {
    return(NA_character_)
  }

  codons <- split_into_codons(seq)

  codon_starts <- seq(1, by = 3, length.out = length(codons))
  codon_ends <- codon_starts + 2

  codon_labels <- purrr::pmap_chr(
    list(codons, codon_starts, codon_ends),
    function(codon, start, end) {
      if (!is.na(boundary_after_nt) && boundary_after_nt >= start && boundary_after_nt < end) {
        cut <- boundary_after_nt - start + 1

        paste0(substr(codon, 1, cut), "/", substr(codon, cut + 1, 3))

      } else if (!is.na(boundary_after_nt) && boundary_after_nt == end) {
        paste0(codon, "/")

      } else if (!is.na(boundary_after_nt) && boundary_after_nt == start - 1) {
        paste0("/", codon)

      } else {
        codon
      }
    }
  )

  paste(codon_labels, collapse = " ")
}

format_aa_codons <- function(seq) {
  codons <- split_into_codons(seq)

  if (length(codons) == 0) {
    return(NA_character_)
  }

  aa <- purrr::map_chr(codons, translate_codon)

  paste(aa, collapse = "   ")
}

aa_has_stop <- function(aa_line) {
  !is.na(aa_line) && stringr::str_detect(aa_line, "\\*")
}

first_stop_codon_index <- function(aa_line) {
  if (is.na(aa_line) || !stringr::str_detect(aa_line, "\\*")) {
    return(NA_integer_)
  }

  aa_vec <- unlist(strsplit(gsub("\\s+", "", aa_line), ""))

  which(aa_vec == "*")[1]
}

adjust_observed_donor_to_exon_boundary <- function(donor_site, strand, coord_type) {
  donor_site <- safe_int(donor_site)
  strand <- as.character(strand)

  if (is.na(donor_site) || !(strand %in% c("+", "-"))) {
    return(NA_integer_)
  }

  if (coord_type == "exon_boundary") {
    return(donor_site)
  }

  if (coord_type == "first_intronic_base") {
    if (strand == "+") return(donor_site - 1L)
    if (strand == "-") return(donor_site + 1L)
  }

  stop("Unknown observed_donor_coord_type: ", coord_type, call. = FALSE)
}

adjust_observed_acceptor_to_exon_boundary <- function(acceptor_site, strand, coord_type) {
  acceptor_site <- safe_int(acceptor_site)
  strand <- as.character(strand)

  if (is.na(acceptor_site) || !(strand %in% c("+", "-"))) {
    return(NA_integer_)
  }

  if (coord_type == "first_exonic_base") {
    return(acceptor_site)
  }

  if (coord_type == "last_intronic_base") {
    if (strand == "+") return(acceptor_site + 1L)
    if (strand == "-") return(acceptor_site - 1L)
  }

  stop("Unknown observed_acceptor_coord_type: ", coord_type, call. = FALSE)
}

compute_phase_shift <- function(ref_donor, ref_acceptor, mut_donor, mut_acceptor, strand) {
  ref_donor <- safe_int(ref_donor)
  ref_acceptor <- safe_int(ref_acceptor)
  mut_donor <- safe_int(mut_donor)
  mut_acceptor <- safe_int(mut_acceptor)
  strand <- as.character(strand)

  if (
    any(is.na(c(ref_donor, ref_acceptor, mut_donor, mut_acceptor))) ||
      !(strand %in% c("+", "-"))
  ) {
    return(tibble(
      upstream_exon_change_nt = NA_integer_,
      downstream_exon_change_nt = NA_integer_,
      net_exonic_change_nt = NA_integer_,
      computed_phase_shift = NA_integer_,
      computed_frameshift = NA
    ))
  }

  if (strand == "+") {
    upstream_change <- mut_donor - ref_donor
    downstream_change <- ref_acceptor - mut_acceptor
  } else {
    upstream_change <- ref_donor - mut_donor
    downstream_change <- mut_acceptor - ref_acceptor
  }

  net_change <- upstream_change + downstream_change
  phase_shift <- ((net_change %% 3) + 3) %% 3

  tibble(
    upstream_exon_change_nt = as.integer(upstream_change),
    downstream_exon_change_nt = as.integer(downstream_change),
    net_exonic_change_nt = as.integer(net_change),
    computed_phase_shift = as.integer(phase_shift),
    computed_frameshift = phase_shift != 0
  )
}


## ============================================================
## STEP 3. FASTA setup and sequence extraction helpers
## No Rsamtools version
## ============================================================

genome <- Biostrings::readDNAStringSet(genome_FASTA)

## Keep first FASTA header token only.
names(genome) <- sub("\\s.*$", "", names(genome))

fa_seqnames <- names(genome)
fa_lengths <- Biostrings::width(genome)
names(fa_lengths) <- fa_seqnames

message("Loaded genome FASTA sequences: ", length(genome))

match_fasta_seqname <- function(chrom) {
  chrom <- as.character(chrom)

  if (chrom %in% fa_seqnames) {
    return(chrom)
  }

  alt1 <- stringr::str_replace_all(chrom, "_", "-")
  if (alt1 %in% fa_seqnames) {
    return(alt1)
  }

  alt2 <- stringr::str_replace_all(chrom, "-", "_")
  if (alt2 %in% fa_seqnames) {
    return(alt2)
  }

  partial_hit <- fa_seqnames[stringr::str_detect(fa_seqnames, fixed(chrom))]

  if (length(partial_hit) == 1) {
    return(partial_hit)
  }

  stop("Chromosome not found in FASTA: ", chrom, call. = FALSE)
}

fetch_oriented_seq <- function(chrom, start, end, strand) {
  chrom2 <- match_fasta_seqname(chrom)

  start <- safe_int(start)
  end <- safe_int(end)
  strand <- as.character(strand)

  if (is.na(start) || is.na(end)) {
    return(NA_character_)
  }

  if (start > end) {
    tmp <- start
    start <- end
    end <- tmp
  }

  start <- max(1L, start)
  end <- min(as.integer(fa_lengths[[chrom2]]), end)

  if (start > end) {
    return(NA_character_)
  }

  seq <- Biostrings::subseq(
    genome[[chrom2]],
    start = start,
    end = end
  )

  if (strand == "-") {
    seq <- Biostrings::reverseComplement(seq)
  }

  as.character(seq)
}

site_window_seq <- function(chrom, site, strand, upstream_before_site, downstream_after_site) {
  site <- safe_int(site)
  strand <- as.character(strand)

  if (is.na(site) || !(strand %in% c("+", "-"))) {
    return(NA_character_)
  }

  if (strand == "+") {
    fetch_oriented_seq(
      chrom = chrom,
      start = site - upstream_before_site,
      end = site + downstream_after_site,
      strand = strand
    )
  } else {
    fetch_oriented_seq(
      chrom = chrom,
      start = site - downstream_after_site,
      end = site + upstream_before_site,
      strand = strand
    )
  }
}

upstream_region_seq <- function(chrom, site, strand, from_upstream, to_upstream) {
  site <- safe_int(site)
  strand <- as.character(strand)

  if (is.na(site) || !(strand %in% c("+", "-"))) {
    return(NA_character_)
  }

  if (strand == "+") {
    fetch_oriented_seq(
      chrom = chrom,
      start = site - from_upstream,
      end = site - to_upstream,
      strand = strand
    )
  } else {
    fetch_oriented_seq(
      chrom = chrom,
      start = site + to_upstream,
      end = site + from_upstream,
      strand = strand
    )
  }
}

extract_splice_motif_from_genome <- function(chrom, donor_boundary, acceptor_boundary, strand) {
  donor_boundary <- safe_int(donor_boundary)
  acceptor_boundary <- safe_int(acceptor_boundary)
  strand <- as.character(strand)

  if (
    any(is.na(c(donor_boundary, acceptor_boundary))) ||
      !(strand %in% c("+", "-"))
  ) {
    return(NA_character_)
  }

  if (strand == "+") {
    donor_dinuc <- fetch_oriented_seq(
      chrom,
      donor_boundary + 1,
      donor_boundary + 2,
      strand
    )

    acceptor_dinuc <- fetch_oriented_seq(
      chrom,
      acceptor_boundary - 2,
      acceptor_boundary - 1,
      strand
    )

  } else {
    donor_dinuc <- fetch_oriented_seq(
      chrom,
      donor_boundary - 2,
      donor_boundary - 1,
      strand
    )

    acceptor_dinuc <- fetch_oriented_seq(
      chrom,
      acceptor_boundary + 1,
      acceptor_boundary + 2,
      strand
    )
  }

  paste0(donor_dinuc, "-", acceptor_dinuc)
}

make_ref_and_mut_junction_seq <- function(chrom,
                                          strand,
                                          ref_donor,
                                          ref_acceptor,
                                          mut_donor,
                                          mut_acceptor,
                                          nt_each_side = 15) {
  ref_donor <- safe_int(ref_donor)
  ref_acceptor <- safe_int(ref_acceptor)
  mut_donor <- safe_int(mut_donor)
  mut_acceptor <- safe_int(mut_acceptor)
  strand <- as.character(strand)

  if (
    any(is.na(c(ref_donor, ref_acceptor, mut_donor, mut_acceptor))) ||
      !(strand %in% c("+", "-"))
  ) {
    return(tibble(
      ref_left_seq = NA_character_,
      ref_right_seq = NA_character_,
      sf3a2_left_seq = NA_character_,
      sf3a2_right_seq = NA_character_,
      ref_junction_seq = NA_character_,
      sf3a2_junction_seq = NA_character_,
      ref_junction_boundary_after_nt = NA_integer_,
      sf3a2_junction_boundary_after_nt = NA_integer_,
      ref_junction_len = NA_integer_,
      sf3a2_junction_len = NA_integer_,
      junction_len_change = NA_integer_
    ))
  }

  if (strand == "+") {
    ref_left_seq <- fetch_oriented_seq(
      chrom,
      ref_donor - nt_each_side + 1,
      ref_donor,
      strand
    )

    ref_right_seq <- fetch_oriented_seq(
      chrom,
      ref_acceptor,
      ref_acceptor + nt_each_side - 1,
      strand
    )

    ## Same reference neighborhood, but cut with SF3A2 donor/acceptor.
    sf3a2_left_seq <- fetch_oriented_seq(
      chrom,
      ref_donor - nt_each_side + 1,
      mut_donor,
      strand
    )

    sf3a2_right_seq <- fetch_oriented_seq(
      chrom,
      mut_acceptor,
      ref_acceptor + nt_each_side - 1,
      strand
    )

  } else {
    ref_left_seq <- fetch_oriented_seq(
      chrom,
      ref_donor,
      ref_donor + nt_each_side - 1,
      strand
    )

    ref_right_seq <- fetch_oriented_seq(
      chrom,
      ref_acceptor - nt_each_side + 1,
      ref_acceptor,
      strand
    )

    sf3a2_left_seq <- fetch_oriented_seq(
      chrom,
      mut_donor,
      ref_donor + nt_each_side - 1,
      strand
    )

    sf3a2_right_seq <- fetch_oriented_seq(
      chrom,
      ref_acceptor - nt_each_side + 1,
      mut_acceptor,
      strand
    )
  }

  ref_junction_seq <- paste0(ref_left_seq, ref_right_seq)
  sf3a2_junction_seq <- paste0(sf3a2_left_seq, sf3a2_right_seq)

  ref_boundary <- nchar(clean_dna(ref_left_seq))
  sf3a2_boundary <- nchar(clean_dna(sf3a2_left_seq))

  ref_len <- nchar(clean_dna(ref_junction_seq))
  sf3a2_len <- nchar(clean_dna(sf3a2_junction_seq))

  tibble(
    ref_left_seq = ref_left_seq,
    ref_right_seq = ref_right_seq,
    sf3a2_left_seq = sf3a2_left_seq,
    sf3a2_right_seq = sf3a2_right_seq,
    ref_junction_seq = ref_junction_seq,
    sf3a2_junction_seq = sf3a2_junction_seq,
    ref_junction_boundary_after_nt = as.integer(ref_boundary),
    sf3a2_junction_boundary_after_nt = as.integer(sf3a2_boundary),
    ref_junction_len = as.integer(ref_len),
    sf3a2_junction_len = as.integer(sf3a2_len),
    junction_len_change = as.integer(sf3a2_len - ref_len)
  )
}


## ============================================================
## STEP 4. Read and standardize junction files
## ============================================================

standardize_junction_file <- function(file, sample_name) {
  raw <- readr::read_tsv(file, show_col_types = FALSE)

  required <- c(
    "isoform",
    "chrom",
    "strand",
    "genomic_start_coord",
    "genomic_end_coord",
    "junction_category",
    "start_site_category",
    "end_site_category",
    "stage"
  )

  missing <- setdiff(required, colnames(raw))

  if (length(missing) > 0) {
    stop(
      "Missing required columns in ",
      file,
      ": ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  out <- raw %>%
    mutate(
      sample = sample_name,
      pbid = as.character(isoform),
      chrom = as.character(chrom),
      strand = as.character(strand),
      stage = as.character(stage),

      genomic_start_coord = safe_int(genomic_start_coord),
      genomic_end_coord = safe_int(genomic_end_coord),

      ## Transcript-oriented donor and acceptor.
      donor_site = ifelse(strand == "+", genomic_start_coord, genomic_end_coord),
      acceptor_site = ifelse(strand == "+", genomic_end_coord, genomic_start_coord),

      donor_site = safe_int(donor_site),
      acceptor_site = safe_int(acceptor_site),

      mut_donor_boundary = purrr::pmap_int(
        list(donor_site, strand),
        function(donor_site, strand) {
          adjust_observed_donor_to_exon_boundary(
            donor_site,
            strand,
            observed_donor_coord_type
          )
        }
      ),

      mut_acceptor_boundary = purrr::pmap_int(
        list(acceptor_site, strand),
        function(acceptor_site, strand) {
          adjust_observed_acceptor_to_exon_boundary(
            acceptor_site,
            strand,
            observed_acceptor_coord_type
          )
        }
      ),

      junction_key = paste(
        chrom,
        strand,
        mut_donor_boundary,
        mut_acceptor_boundary,
        sep = ":"
      ),

      gene_id_from_associated_transcript = if ("associated_transcript" %in% colnames(.)) {
        extract_pf_gene_id(associated_transcript)
      } else {
        NA_character_
      },

      strict_novel = case_when(
        strict_requires_junction_category_novel &&
          strict_requires_at_least_one_novel_site ~
          junction_category == "novel" &
          (start_site_category == "novel" | end_site_category == "novel"),

        strict_requires_junction_category_novel &&
          !strict_requires_at_least_one_novel_site ~
          junction_category == "novel",

        TRUE ~ junction_category == "novel"
      )
    )

  if (exclude_rts_junction && "rts_junction" %in% colnames(out)) {
    out <- out %>%
      filter(!(tolower(as.character(rts_junction)) %in% c("true", "t", "1", "yes", "y")))
  }

  if (exclude_indel_near_junction && "indel_near_junct" %in% colnames(out)) {
    out <- out %>%
      filter(!(tolower(as.character(indel_near_junct)) %in% c("true", "t", "1", "yes", "y")))
  }

  if (exclude_bite_junction && "bite_junction" %in% colnames(out)) {
    out <- out %>%
      filter(!(tolower(as.character(bite_junction)) %in% c("true", "t", "1", "yes", "y")))
  }

  out
}

sf3a2_junc <- standardize_junction_file(sf3a2_file, "SF3A2")
nf54_junc <- standardize_junction_file(nf54_file, "NF54")

message("SF3A2 rows loaded: ", nrow(sf3a2_junc))
message("NF54 rows loaded: ", nrow(nf54_junc))


## ============================================================
## STEP 5. Build target-gene genomic ranges from reference annotation
## ============================================================

make_ref_target_gene_ranges <- function(ref_gtf, target_genes) {
  ref_gr <- rtracklayer::import(ref_gtf)

  ref_df <- as.data.frame(ref_gr)

  meta_cols <- setdiff(
    colnames(ref_df),
    c("seqnames", "start", "end", "width", "strand")
  )

  ref_df <- ref_df %>%
    mutate(
      annotation_text = apply(
        select(., all_of(meta_cols)),
        1,
        function(x) paste(as.character(x), collapse = ";")
      ),
      gene_id_clean = extract_pf_gene_id(annotation_text),
      type_clean = if ("type" %in% colnames(.)) as.character(type) else NA_character_
    ) %>%
    filter(
      !is.na(gene_id_clean),
      gene_id_clean %in% target_genes
    )

  if (nrow(ref_df) == 0) {
    stop("No target genes found in reference annotation.", call. = FALSE)
  }

  chosen <- ref_df %>%
    filter(type_clean == "gene")

  if (nrow(chosen) == 0) {
    chosen <- ref_df %>%
      filter(type_clean %in% c("mRNA", "transcript"))
  }

  if (nrow(chosen) == 0) {
    chosen <- ref_df %>%
      filter(type_clean %in% c("exon", "CDS"))
  }

  if (nrow(chosen) == 0) {
    chosen <- ref_df
  }

  gr <- GenomicRanges::GRanges(
    seqnames = chosen$seqnames,
    ranges = IRanges::IRanges(start = chosen$start, end = chosen$end),
    strand = chosen$strand
  )

  gr$gene_id <- chosen$gene_id_clean

  grl <- split(gr, gr$gene_id)

  red <- GenomicRanges::reduce(grl, ignore.strand = FALSE)

  gene_ranges <- unlist(red, use.names = FALSE)
  gene_ranges$gene_id <- rep(names(red), elementNROWS(red))

  gene_ranges
}

target_gene_ranges <- make_ref_target_gene_ranges(ref_gtf, target_genes)

message("Target-gene ranges created: ", length(target_gene_ranges))
message("Target genes represented: ", length(unique(target_gene_ranges$gene_id)))


## ============================================================
## STEP 6. Rescue/assign gene IDs by coordinate overlap
## ============================================================

assign_gene_by_coordinate <- function(junc_tbl, target_gene_ranges) {
  junc_gr <- GenomicRanges::GRanges(
    seqnames = junc_tbl$chrom,
    ranges = IRanges::IRanges(
      start = pmin(junc_tbl$mut_donor_boundary, junc_tbl$mut_acceptor_boundary, na.rm = TRUE),
      end = pmax(junc_tbl$mut_donor_boundary, junc_tbl$mut_acceptor_boundary, na.rm = TRUE)
    ),
    strand = junc_tbl$strand
  )

  junc_gr$row_id_for_gene_overlap <- seq_len(nrow(junc_tbl))

  hits <- GenomicRanges::findOverlaps(
    junc_gr,
    target_gene_ranges,
    ignore.strand = FALSE
  )

  hit_tbl <- tibble(
    row_id_for_gene_overlap = junc_gr$row_id_for_gene_overlap[queryHits(hits)],
    gene_id_by_coord = target_gene_ranges$gene_id[subjectHits(hits)]
  ) %>%
    distinct(row_id_for_gene_overlap, gene_id_by_coord)

  junc_tbl %>%
    mutate(row_id_for_gene_overlap = row_number()) %>%
    left_join(hit_tbl, by = "row_id_for_gene_overlap") %>%
    mutate(
      gene_id = case_when(
        !is.na(gene_id_from_associated_transcript) ~ gene_id_from_associated_transcript,
        !is.na(gene_id_by_coord) ~ gene_id_by_coord,
        TRUE ~ NA_character_
      ),

      gene_assignment_method = case_when(
        !is.na(gene_id_from_associated_transcript) ~ "associated_transcript",
        is.na(gene_id_from_associated_transcript) & !is.na(gene_id_by_coord) ~ "coordinate_overlap",
        TRUE ~ "unassigned"
      )
    )
}

sf3a2_junc <- assign_gene_by_coordinate(sf3a2_junc, target_gene_ranges)
nf54_junc <- assign_gene_by_coordinate(nf54_junc, target_gene_ranges)

sf3a2_junc %>%
  count(gene_assignment_method, sort = TRUE) %>%
  print(n = 20)


## ============================================================
## STEP 7. Filter to SF3A2-specific strict-novel junctions
## ============================================================

nf54_all_junction_keys <- unique(nf54_junc$junction_key)

sf3a2_strict_only <- sf3a2_junc %>%
  filter(
    gene_id %in% target_genes,
    stage %in% target_stages,
    strict_novel,
    !junction_key %in% nf54_all_junction_keys
  ) %>%
  mutate(
    sf3a2_only_strict_novel = TRUE
  )

message("SF3A2-only strict-novel rows: ", nrow(sf3a2_strict_only))
message("Unique SF3A2-only strict-novel junctions: ", n_distinct(sf3a2_strict_only$junction_key))
message("Genes represented: ", n_distinct(sf3a2_strict_only$gene_id))

sf3a2_strict_only %>%
  count(gene_id, stage, sort = TRUE) %>%
  print(n = 100)

if (nrow(sf3a2_strict_only) == 0) {
  stop("No SF3A2-specific strict-novel junctions found under the current filters.", call. = FALSE)
}


## ============================================================
## STEP 8. Build reference exon table and nearest exon-exon context
## ============================================================

make_reference_exon_table <- function(ref_gtf, target_genes) {
  ref_gr <- rtracklayer::import(ref_gtf)

  ref_df <- as.data.frame(ref_gr)

  meta_cols <- setdiff(
    colnames(ref_df),
    c("seqnames", "start", "end", "width", "strand")
  )

  ref_df <- ref_df %>%
    mutate(
      annotation_text = apply(
        select(., all_of(meta_cols)),
        1,
        function(x) paste(as.character(x), collapse = ";")
      ),
      gene_id_clean = extract_pf_gene_id(annotation_text),
      type_clean = if ("type" %in% colnames(.)) as.character(type) else NA_character_
    ) %>%
    filter(
      type_clean == "exon",
      !is.na(gene_id_clean),
      gene_id_clean %in% target_genes
    )

  if (nrow(ref_df) == 0) {
    stop("No target-gene exon rows found in reference GFF/GTF.", call. = FALSE)
  }

  tx_candidates <- c("transcript_id", "Parent", "ID", "Name")

  tx_col <- tx_candidates[tx_candidates %in% colnames(ref_df)][1]

  if (is.na(tx_col)) {
    ref_df <- ref_df %>%
      mutate(transcript_id_clean = gene_id_clean)
  } else {
    ref_df <- ref_df %>%
      mutate(transcript_id_clean = as.character(.data[[tx_col]]))
  }

  tx_choice <- ref_df %>%
    group_by(gene_id_clean, transcript_id_clean) %>%
    summarise(
      total_exon_bp = sum(width, na.rm = TRUE),
      exon_count = n(),
      .groups = "drop"
    ) %>%
    arrange(gene_id_clean, desc(total_exon_bp), desc(exon_count), transcript_id_clean) %>%
    group_by(gene_id_clean) %>%
    slice_head(n = 1) %>%
    ungroup()

  ref_exons <- ref_df %>%
    inner_join(tx_choice, by = c("gene_id_clean", "transcript_id_clean")) %>%
    transmute(
      gene_id = gene_id_clean,
      transcript_id = transcript_id_clean,
      chrom = as.character(seqnames),
      strand = as.character(strand),
      exon_start = safe_int(start),
      exon_end = safe_int(end)
    )

  ref_exons <- ref_exons %>%
    group_by(gene_id, transcript_id) %>%
    group_modify(function(.x, .y) {
      strd <- .x$strand[1]

      if (strd == "+") {
        .x <- .x %>% arrange(exon_start)
      } else {
        .x <- .x %>% arrange(desc(exon_start))
      }

      .x %>%
        mutate(
          exon_number = row_number(),
          exon_label = paste0("Exon ", exon_number),
          donor_ref_boundary_coord = ifelse(strand == "+", exon_end, exon_start),
          acceptor_ref_boundary_coord = ifelse(strand == "+", exon_start, exon_end)
        )
    }) %>%
    ungroup()

  ref_exons
}

ref_exons <- make_reference_exon_table(ref_gtf, target_genes)

annotate_reference_context_one <- function(gene_id, donor_boundary, acceptor_boundary, ref_exons) {
  ex <- ref_exons %>%
    filter(gene_id == !!gene_id)

  if (nrow(ex) == 0) {
    return(tibble(
      ref_transcript_id = NA_character_,

      donor_ref_exon_number = NA_integer_,
      donor_ref_exon_start = NA_integer_,
      donor_ref_exon_end = NA_integer_,
      donor_ref_boundary_coord = NA_integer_,
      donor_ref_boundary_offset_nt = NA_integer_,

      acceptor_ref_exon_number = NA_integer_,
      acceptor_ref_exon_start = NA_integer_,
      acceptor_ref_exon_end = NA_integer_,
      acceptor_ref_boundary_coord = NA_integer_,
      acceptor_ref_boundary_offset_nt = NA_integer_,

      reference_exon_exon_context = NA_character_,
      reference_exon_exon_context_with_coords = NA_character_
    ))
  }

  donor_hit <- ex %>%
    mutate(offset = abs(safe_int(donor_boundary) - donor_ref_boundary_coord)) %>%
    arrange(offset) %>%
    slice_head(n = 1)

  acceptor_hit <- ex %>%
    mutate(offset = abs(safe_int(acceptor_boundary) - acceptor_ref_boundary_coord)) %>%
    arrange(offset) %>%
    slice_head(n = 1)

  tibble(
    ref_transcript_id = donor_hit$transcript_id[1],

    donor_ref_exon_number = donor_hit$exon_number[1],
    donor_ref_exon_start = donor_hit$exon_start[1],
    donor_ref_exon_end = donor_hit$exon_end[1],
    donor_ref_boundary_coord = donor_hit$donor_ref_boundary_coord[1],
    donor_ref_boundary_offset_nt = safe_int(donor_boundary) - donor_hit$donor_ref_boundary_coord[1],

    acceptor_ref_exon_number = acceptor_hit$exon_number[1],
    acceptor_ref_exon_start = acceptor_hit$exon_start[1],
    acceptor_ref_exon_end = acceptor_hit$exon_end[1],
    acceptor_ref_boundary_coord = acceptor_hit$acceptor_ref_boundary_coord[1],
    acceptor_ref_boundary_offset_nt = safe_int(acceptor_boundary) - acceptor_hit$acceptor_ref_boundary_coord[1],

    reference_exon_exon_context = paste0(
      "Exon ",
      donor_hit$exon_number[1],
      " -> Exon ",
      acceptor_hit$exon_number[1]
    ),

    reference_exon_exon_context_with_coords = paste0(
      "Exon ",
      donor_hit$exon_number[1],
      " [",
      donor_hit$exon_start[1],
      "-",
      donor_hit$exon_end[1],
      "] -> Exon ",
      acceptor_hit$exon_number[1],
      " [",
      acceptor_hit$exon_start[1],
      "-",
      acceptor_hit$exon_end[1],
      "]"
    )
  )
}

ref_context_tbl <- sf3a2_strict_only %>%
  select(gene_id, junction_key, mut_donor_boundary, mut_acceptor_boundary) %>%
  distinct() %>%
  rowwise() %>%
  do(
    bind_cols(
      tibble(
        gene_id = .$gene_id,
        junction_key = .$junction_key,
        mut_donor_boundary = .$mut_donor_boundary,
        mut_acceptor_boundary = .$mut_acceptor_boundary
      ),
      annotate_reference_context_one(
        gene_id = .$gene_id,
        donor_boundary = .$mut_donor_boundary,
        acceptor_boundary = .$mut_acceptor_boundary,
        ref_exons = ref_exons
      )
    )
  ) %>%
  ungroup()

sf3a2_strict_only <- sf3a2_strict_only %>%
  left_join(
    ref_context_tbl,
    by = c("gene_id", "junction_key", "mut_donor_boundary", "mut_acceptor_boundary")
  )


## ============================================================
## STEP 9. Add donor, acceptor, branchpoint, and splice motif windows
## ============================================================

sf3a2_with_windows <- sf3a2_strict_only %>%
  rowwise() %>%
  mutate(
    splice_motif_from_genome = extract_splice_motif_from_genome(
      chrom = chrom,
      donor_boundary = mut_donor_boundary,
      acceptor_boundary = mut_acceptor_boundary,
      strand = strand
    ),

    donor_5ss_seq_5exon_20intron = site_window_seq(
      chrom = chrom,
      site = mut_donor_boundary,
      strand = strand,
      upstream_before_site = donor_exonic_nt - 1,
      downstream_after_site = donor_intronic_nt
    ),

    acceptor_3ss_seq_60intron_10exon = site_window_seq(
      chrom = chrom,
      site = mut_acceptor_boundary,
      strand = strand,
      upstream_before_site = acceptor_intronic_nt,
      downstream_after_site = acceptor_exonic_nt - 1
    ),

    acceptor_core_seq_20intron_5exon = site_window_seq(
      chrom = chrom,
      site = mut_acceptor_boundary,
      strand = strand,
      upstream_before_site = acceptor_core_upstream_nt,
      downstream_after_site = acceptor_core_exonic_nt - 1
    ),

    branchpoint_candidate_seq_60to15nt_upstream_3ss = upstream_region_seq(
      chrom = chrom,
      site = mut_acceptor_boundary,
      strand = strand,
      from_upstream = bp_upstream_from,
      to_upstream = bp_upstream_to
    ),

    bp_candidate_A_fraction = base_fraction(
      branchpoint_candidate_seq_60to15nt_upstream_3ss,
      bases = c("A")
    ),

    bp_candidate_pyrimidine_fraction = base_fraction(
      branchpoint_candidate_seq_60to15nt_upstream_3ss,
      bases = c("C", "T")
    ),

    acceptor_core_A_fraction = base_fraction(
      acceptor_core_seq_20intron_5exon,
      bases = c("A")
    ),

    acceptor_core_pyrimidine_fraction = base_fraction(
      acceptor_core_seq_20intron_5exon,
      bases = c("C", "T")
    )
  ) %>%
  ungroup()


## ============================================================
## STEP 10. Add matched reference splice-site windows
## ============================================================

sf3a2_with_ref_windows <- sf3a2_with_windows %>%
  rowwise() %>%
  mutate(
    ref_splice_motif_from_genome = extract_splice_motif_from_genome(
      chrom = chrom,
      donor_boundary = donor_ref_boundary_coord,
      acceptor_boundary = acceptor_ref_boundary_coord,
      strand = strand
    ),

    ref_donor_5ss_seq_5exon_20intron = site_window_seq(
      chrom = chrom,
      site = donor_ref_boundary_coord,
      strand = strand,
      upstream_before_site = donor_exonic_nt - 1,
      downstream_after_site = donor_intronic_nt
    ),

    ref_acceptor_3ss_seq_60intron_10exon = site_window_seq(
      chrom = chrom,
      site = acceptor_ref_boundary_coord,
      strand = strand,
      upstream_before_site = acceptor_intronic_nt,
      downstream_after_site = acceptor_exonic_nt - 1
    ),

    ref_acceptor_core_seq_20intron_5exon = site_window_seq(
      chrom = chrom,
      site = acceptor_ref_boundary_coord,
      strand = strand,
      upstream_before_site = acceptor_core_upstream_nt,
      downstream_after_site = acceptor_core_exonic_nt - 1
    ),

    ref_branchpoint_candidate_seq_60to15nt_upstream_3ss = upstream_region_seq(
      chrom = chrom,
      site = acceptor_ref_boundary_coord,
      strand = strand,
      from_upstream = bp_upstream_from,
      to_upstream = bp_upstream_to
    ),

    ref_bp_candidate_A_fraction = base_fraction(
      ref_branchpoint_candidate_seq_60to15nt_upstream_3ss,
      bases = c("A")
    ),

    ref_bp_candidate_pyrimidine_fraction = base_fraction(
      ref_branchpoint_candidate_seq_60to15nt_upstream_3ss,
      bases = c("C", "T")
    ),

    ref_acceptor_core_A_fraction = base_fraction(
      ref_acceptor_core_seq_20intron_5exon,
      bases = c("A")
    ),

    ref_acceptor_core_pyrimidine_fraction = base_fraction(
      ref_acceptor_core_seq_20intron_5exon,
      bases = c("C", "T")
    )
  ) %>%
  ungroup() %>%
  mutate(
    delta_bp_A_fraction = bp_candidate_A_fraction - ref_bp_candidate_A_fraction,
    delta_bp_pyrimidine_fraction = bp_candidate_pyrimidine_fraction - ref_bp_candidate_pyrimidine_fraction,
    delta_acceptor_core_A_fraction = acceptor_core_A_fraction - ref_acceptor_core_A_fraction,
    delta_acceptor_core_pyrimidine_fraction = acceptor_core_pyrimidine_fraction - ref_acceptor_core_pyrimidine_fraction
  )


## ============================================================
## STEP 11. Add direct sequence comparison to nearest reference junction
## ============================================================

junction_seq_tbl <- sf3a2_with_ref_windows %>%
  select(
    gene_id,
    junction_key,
    chrom,
    strand,
    donor_ref_boundary_coord,
    acceptor_ref_boundary_coord,
    mut_donor_boundary,
    mut_acceptor_boundary
  ) %>%
  distinct() %>%
  rowwise() %>%
  do(
    bind_cols(
      tibble(
        gene_id = .$gene_id,
        junction_key = .$junction_key,
        chrom = .$chrom,
        strand = .$strand,
        donor_ref_boundary_coord = .$donor_ref_boundary_coord,
        acceptor_ref_boundary_coord = .$acceptor_ref_boundary_coord,
        mut_donor_boundary = .$mut_donor_boundary,
        mut_acceptor_boundary = .$mut_acceptor_boundary
      ),

      make_ref_and_mut_junction_seq(
        chrom = .$chrom,
        strand = .$strand,
        ref_donor = .$donor_ref_boundary_coord,
        ref_acceptor = .$acceptor_ref_boundary_coord,
        mut_donor = .$mut_donor_boundary,
        mut_acceptor = .$mut_acceptor_boundary,
        nt_each_side = ref_junction_nt_each_side
      )
    )
  ) %>%
  ungroup() %>%
  mutate(
    ref_junction_dna_codon_line = purrr::map2_chr(
      ref_junction_seq,
      ref_junction_boundary_after_nt,
      format_dna_codons_with_boundary
    ),

    sf3a2_junction_dna_codon_line = purrr::map2_chr(
      sf3a2_junction_seq,
      sf3a2_junction_boundary_after_nt,
      format_dna_codons_with_boundary
    ),

    ref_junction_aa_line = purrr::map_chr(ref_junction_seq, format_aa_codons),
    sf3a2_junction_aa_line = purrr::map_chr(sf3a2_junction_seq, format_aa_codons),

    ref_junction_has_stop = purrr::map_lgl(ref_junction_aa_line, aa_has_stop),
    sf3a2_junction_has_stop = purrr::map_lgl(sf3a2_junction_aa_line, aa_has_stop),

    sf3a2_first_stop_codon_index = purrr::map_int(
      sf3a2_junction_aa_line,
      first_stop_codon_index
    ),

    direct_sequence_changed = ref_junction_seq != sf3a2_junction_seq,
    direct_translation_changed = ref_junction_aa_line != sf3a2_junction_aa_line
  )

sf3a2_summary <- sf3a2_with_ref_windows %>%
  left_join(
    junction_seq_tbl,
    by = c(
      "gene_id",
      "junction_key",
      "chrom",
      "strand",
      "donor_ref_boundary_coord",
      "acceptor_ref_boundary_coord",
      "mut_donor_boundary",
      "mut_acceptor_boundary"
    )
  )


## ============================================================
## STEP 12. Add predicted coding-frame effects
## ============================================================

frame_tbl <- sf3a2_summary %>%
  select(
    gene_id,
    junction_key,
    strand,
    donor_ref_boundary_coord,
    acceptor_ref_boundary_coord,
    mut_donor_boundary,
    mut_acceptor_boundary
  ) %>%
  distinct() %>%
  rowwise() %>%
  do(
    bind_cols(
      tibble(
        gene_id = .$gene_id,
        junction_key = .$junction_key
      ),

      compute_phase_shift(
        ref_donor = .$donor_ref_boundary_coord,
        ref_acceptor = .$acceptor_ref_boundary_coord,
        mut_donor = .$mut_donor_boundary,
        mut_acceptor = .$mut_acceptor_boundary,
        strand = .$strand
      )
    )
  ) %>%
  ungroup() %>%
  mutate(
    predicted_frame_effect = case_when(
      is.na(computed_frameshift) ~ "undetermined",
      computed_frameshift ~ "frameshift_predicted",
      !computed_frameshift & net_exonic_change_nt == 0 ~ "same_length_or_boundary_shift_only",
      !computed_frameshift & net_exonic_change_nt != 0 ~ "in_frame_gain_or_loss",
      TRUE ~ "undetermined"
    ),

    predicted_stop_effect = case_when(
      sf3a2_first_stop_codon_index > 0 ~ paste0("STOP_detected_at_local_codon_", sf3a2_first_stop_codon_index),
      sf3a2_junction_has_stop ~ "STOP_detected",
      TRUE ~ "no_STOP_in_displayed_window"
    )
  )

sf3a2_summary <- sf3a2_summary %>%
  left_join(frame_tbl, by = c("gene_id", "junction_key"))


## ============================================================
## STEP 13. Reorder and write full output table
## ============================================================

full_cols_first <- c(
  "sample",
  "gene_id",
  "stage",
  "pbid",
  "chrom",
  "strand",

  "donor_site",
  "acceptor_site",
  "mut_donor_boundary",
  "mut_acceptor_boundary",
  "junction_key",

  "pbid_fl_support",
  "fl_support",
  "pbid_fl_cpm",
  "pbid_fl_fraction",

  "junction_category",
  "start_site_category",
  "end_site_category",
  "strict_novel",
  "sf3a2_only_strict_novel",
  "splice_site",
  "canonical",
  "splice_motif_from_genome",

  "structural_category",
  "subcategory",
  "filter_result",
  "associated_transcript",
  "gene_assignment_method",

  "ref_transcript_id",
  "reference_exon_exon_context",
  "reference_exon_exon_context_with_coords",

  "donor_ref_exon_number",
  "donor_ref_exon_start",
  "donor_ref_exon_end",
  "donor_ref_boundary_coord",
  "donor_ref_boundary_offset_nt",

  "acceptor_ref_exon_number",
  "acceptor_ref_exon_start",
  "acceptor_ref_exon_end",
  "acceptor_ref_boundary_coord",
  "acceptor_ref_boundary_offset_nt",

  "donor_5ss_seq_5exon_20intron",
  "acceptor_3ss_seq_60intron_10exon",
  "acceptor_core_seq_20intron_5exon",
  "branchpoint_candidate_seq_60to15nt_upstream_3ss",

  "ref_donor_5ss_seq_5exon_20intron",
  "ref_acceptor_3ss_seq_60intron_10exon",
  "ref_acceptor_core_seq_20intron_5exon",
  "ref_branchpoint_candidate_seq_60to15nt_upstream_3ss",

  "bp_candidate_A_fraction",
  "bp_candidate_pyrimidine_fraction",
  "acceptor_core_A_fraction",
  "acceptor_core_pyrimidine_fraction",

  "ref_bp_candidate_A_fraction",
  "ref_bp_candidate_pyrimidine_fraction",
  "ref_acceptor_core_A_fraction",
  "ref_acceptor_core_pyrimidine_fraction",

  "delta_bp_A_fraction",
  "delta_bp_pyrimidine_fraction",
  "delta_acceptor_core_A_fraction",
  "delta_acceptor_core_pyrimidine_fraction",

  "ref_junction_seq",
  "sf3a2_junction_seq",
  "ref_junction_len",
  "sf3a2_junction_len",
  "junction_len_change",

  "ref_junction_dna_codon_line",
  "sf3a2_junction_dna_codon_line",
  "ref_junction_aa_line",
  "sf3a2_junction_aa_line",

  "direct_sequence_changed",
  "direct_translation_changed",

  "upstream_exon_change_nt",
  "downstream_exon_change_nt",
  "net_exonic_change_nt",
  "computed_phase_shift",
  "computed_frameshift",
  "predicted_frame_effect",
  "predicted_stop_effect",

  "ref_junction_has_stop",
  "sf3a2_junction_has_stop",
  "sf3a2_first_stop_codon_index"
)

full_output <- sf3a2_summary %>%
  select(any_of(full_cols_first), everything()) %>%
  arrange(gene_id, stage, desc(safe_num(pbid_fl_support)), pbid, junction_key)

full_outfile <- file.path(
  out_dir,
  "SF3A2_specific_strict_novel_SP_EF_EM_splice_junction_summary_full.tsv"
)

readr::write_tsv(full_output, full_outfile)

message("Full output written to:")
message(full_outfile)


## ============================================================
## STEP 14. Write compact summary table and FASTA files
## ============================================================

compact_output <- full_output %>%
  transmute(
    gene_id,
    stage,
    pbid,

    pbid_fl_support = if ("pbid_fl_support" %in% colnames(full_output)) pbid_fl_support else NA,
    fl_support = if ("fl_support" %in% colnames(full_output)) fl_support else NA,
    pbid_fl_cpm = if ("pbid_fl_cpm" %in% colnames(full_output)) pbid_fl_cpm else NA,
    pbid_fl_fraction = if ("pbid_fl_fraction" %in% colnames(full_output)) pbid_fl_fraction else NA,

    chrom,
    strand,

    sf3a2_junction = paste0(mut_donor_boundary, " -> ", mut_acceptor_boundary),
    nearest_reference_junction = paste0(donor_ref_boundary_coord, " -> ", acceptor_ref_boundary_coord),

    reference_exon_exon_context,
    reference_exon_exon_context_with_coords,

    junction_category,
    start_site_category,
    end_site_category,
    structural_category,

    splice_site = if ("splice_site" %in% colnames(full_output)) splice_site else NA,
    canonical = if ("canonical" %in% colnames(full_output)) canonical else NA,

    splice_motif_from_genome,
    ref_splice_motif_from_genome,

    donor_5ss_seq_5exon_20intron,
    acceptor_3ss_seq_60intron_10exon,
    acceptor_core_seq_20intron_5exon,
    branchpoint_candidate_seq_60to15nt_upstream_3ss,

    ref_donor_5ss_seq_5exon_20intron,
    ref_acceptor_3ss_seq_60intron_10exon,
    ref_acceptor_core_seq_20intron_5exon,
    ref_branchpoint_candidate_seq_60to15nt_upstream_3ss,

    bp_candidate_A_fraction,
    bp_candidate_pyrimidine_fraction,
    acceptor_core_A_fraction,
    acceptor_core_pyrimidine_fraction,

    ref_bp_candidate_A_fraction,
    ref_bp_candidate_pyrimidine_fraction,
    ref_acceptor_core_A_fraction,
    ref_acceptor_core_pyrimidine_fraction,

    delta_bp_A_fraction,
    delta_bp_pyrimidine_fraction,
    delta_acceptor_core_A_fraction,
    delta_acceptor_core_pyrimidine_fraction,

    ref_junction_seq,
    sf3a2_junction_seq,
    ref_junction_len,
    sf3a2_junction_len,
    junction_len_change,

    ref_junction_dna_codon_line,
    ref_junction_aa_line,
    sf3a2_junction_dna_codon_line,
    sf3a2_junction_aa_line,

    net_exonic_change_nt,
    computed_phase_shift,
    computed_frameshift,
    predicted_frame_effect,
    predicted_stop_effect,

    direct_sequence_changed,
    direct_translation_changed,
    ref_junction_has_stop,
    sf3a2_junction_has_stop,
    sf3a2_first_stop_codon_index
  )

compact_outfile <- file.path(
  out_dir,
  "SF3A2_specific_strict_novel_SP_EF_EM_splice_junction_summary_compact.tsv"
)

readr::write_tsv(compact_output, compact_outfile)

message("Compact output written to:")
message(compact_outfile)

## FASTA outputs

ref_sequences_clean <- purrr::map_chr(full_output$ref_junction_seq, clean_dna)
sf3a2_sequences_clean <- purrr::map_chr(full_output$sf3a2_junction_seq, clean_dna)

ref_sequences_clean[is.na(ref_sequences_clean) | ref_sequences_clean == ""] <- "N"
sf3a2_sequences_clean[is.na(sf3a2_sequences_clean) | sf3a2_sequences_clean == ""] <- "N"

ref_fa <- Biostrings::DNAStringSet(ref_sequences_clean)

names(ref_fa) <- paste(
  full_output$gene_id,
  full_output$stage,
  full_output$pbid,
  full_output$reference_exon_exon_context,
  "REF",
  sep = "|"
)

sf3a2_fa <- Biostrings::DNAStringSet(sf3a2_sequences_clean)

names(sf3a2_fa) <- paste(
  full_output$gene_id,
  full_output$stage,
  full_output$pbid,
  full_output$reference_exon_exon_context,
  "SF3A2",
  sep = "|"
)

ref_fasta_out <- file.path(
  out_dir,
  "SF3A2_specific_strict_novel_matched_reference_junction_sequences.fasta"
)

sf3a2_fasta_out <- file.path(
  out_dir,
  "SF3A2_specific_strict_novel_mutant_junction_sequences.fasta"
)

Biostrings::writeXStringSet(ref_fa, filepath = ref_fasta_out)
Biostrings::writeXStringSet(sf3a2_fa, filepath = sf3a2_fasta_out)

message("Reference FASTA written to:")
message(ref_fasta_out)

message("SF3A2 FASTA written to:")
message(sf3a2_fasta_out)

message("Done.")