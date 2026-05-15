# Pipeline overview

## Short-read analysis

The short-read pipeline is divided into preprocessing, integration, annotation, expression visualization, and pseudotime analysis.

Important correction: pseudotime is handled as **two separate workflows**:

1. **Asexual pseudotime**  
   Uses ER/LR/ET/LT/ES/LS cells.  
   Script: `short_reads/scripts/07a_slingshot_tradeseq_asexual.R`

2. **Sexual/gametocyte pseudotime**  
   Uses SP/EF/EM cells from the induced/sexual gametocyte object.  
   Script: `short_reads/scripts/07b_slingshot_tradeseq_gametocyte.R`

The corresponding heatmaps are also separated:

- `08a_pseudotime_heatmap_asexual.R`
- `08b_pseudotime_heatmap_gametocyte.R`

## Long-read analysis

The long-read workflow contains independent scripts for:

1. Stage-labeling SQANTI3 classification files
2. DEXSeq exon-usage analysis from pseudobulk BAMs
3. SQANTI3 junction comparison
4. Strict novel-junction filtering using PBID FL support
5. Cis-element BED/FASTA preparation
6. STREME/FIMO motif analysis
7. Motif and sequence-composition summary
8. Branchpoint/acceptor logo plotting
9. Target-gene cis-context extraction


## Final strict-novel target-junction summary

`long_reads/scripts/09_sf3a2_strict_novel_sp_ef_em_summary.R` is the final target-junction summary script. It focuses on SF3A2-specific strict-novel splice junctions in SP, EF, and EM stages. It reports PBID support, splice-site novelty, donor/acceptor/branchpoint-proximal sequence windows, nearest reference exon-exon context, matched REF versus SF3A2 junction sequences, and predicted local coding-frame effects.
