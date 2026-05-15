#!/usr/bin/env Rscript
# ============================================================
# DEXSeq exon-usage analysis from stage-specific pseudobulk BAMs.
# Usage:
#   Rscript long_reads/scripts/02_dexseq_pseudobulk_exon_usage.R Progenitor
# ============================================================

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(GenomicFeatures)
  library(txdbmaker)
  library(GenomicAlignments)
  library(Rsamtools)
  library(DEXSeq)
  library(SummarizedExperiment)
})

config_file <- Sys.getenv("SF3A2_CONFIG", unset = "config/config.yml")
if (!file.exists(config_file)) config_file <- "config/config_template.yml"
config <- yaml::read_yaml(config_file)

stage_tag <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(stage_tag)) {
  stop("Provide a stage tag, e.g. Progenitor, female, male, earlyrings, laterings.")
}

bam_map <- read.csv(config$paths$long_reads$stage_bam_map, stringsAsFactors = FALSE)
stage_bams <- bam_map %>% filter(stage == stage_tag)

if (nrow(stage_bams) == 0) stop("No BAMs found for stage: ", stage_tag)

bam_files <- stage_bams$bam
if (!all(file.exists(bam_files))) {
  stop("Missing BAM files: ", paste(bam_files[!file.exists(bam_files)], collapse = ", "))
}

txdb <- txdbmaker::makeTxDbFromGFF(config$paths$long_reads$gtf)

flattened <- GenomicFeatures::exonicParts(txdb, linked.to.single.gene.only = TRUE)
names(flattened) <- sprintf("%s:E%0.3d", flattened$gene_id, flattened$exonic_part)

se <- GenomicAlignments::summarizeOverlaps(
  features = flattened,
  reads = Rsamtools::BamFileList(bam_files),
  mode = "Union",
  singleEnd = TRUE,
  fragments = FALSE,
  ignore.strand = FALSE
)

colnames(se) <- paste(stage_bams$condition, paste0("rep", stage_bams$replicate), sep = "_")

sample_table <- data.frame(
  row.names = colnames(se),
  condition = factor(stage_bams$condition, levels = c("NF54", "SF3A2")),
  replicate = factor(stage_bams$replicate)
)

dxd <- DEXSeq::DEXSeqDataSet(
  countData = assays(se)$counts,
  sampleData = sample_table,
  design = ~ sample + exon + condition:exon,
  featureID = sub(".*:", "", names(flattened)),
  groupID = flattened$gene_id
)

dxd <- DEXSeq::estimateSizeFactors(dxd)
dxd <- DEXSeq::estimateDispersions(dxd)
dxd <- DEXSeq::testForDEU(dxd)
dxd <- DEXSeq::estimateExonFoldChanges(dxd, fitExpToVar = "condition")

dxr <- DEXSeq::DEXSeqResults(dxd)
res <- as.data.frame(dxr)
res$stage <- stage_tag

outdir <- file.path(config$paths$long_reads$tables, "dexseq")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

write.csv(res, file.path(outdir, paste0("DEXSeq_", stage_tag, "_all_exonic_parts.csv")), row.names = FALSE)

sig <- res %>%
  filter(!is.na(padj), padj <= 0.05)

write.csv(sig, file.path(outdir, paste0("DEXSeq_", stage_tag, "_significant_exonic_parts.csv")), row.names = FALSE)

message("Significant DEXSeq exonic parts for ", stage_tag, ": ", nrow(sig))
