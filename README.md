# SF3A2 transcriptomics analysis workflows

This repository contains modular analysis scripts for the SF3A2 *Plasmodium falciparum* transcriptomics project.

The repository is organized into two major modules:

1. **Short-read single-cell RNA-seq analysis**
2. **Long-read PacBio Iso-Seq / SQANTI3 / splice-junction analysis**

The short-read module includes separate Slingshot/tradeSeq workflows for:

- **Asexual-stage pseudotime**
- **Sexual/gametocyte pseudotime**

The long-read module includes workflows for:

- Adding stage labels to SQANTI3 classification outputs
- DEXSeq-based exon-usage analysis from pseudobulk BAM files
- SQANTI3 junction comparison between NF54 and SF3A2
- Strict novel-junction filtering
- PBID-supported junction filtering
- Cis-element window preparation
- STREME/FIMO motif analysis
- Branchpoint/acceptor sequence-logo visualization
- Target-gene cis-context analysis for candidate genes such as `PF3D7_1438800` and `PF3D7_1343000`
- SF3A2-only strict-novel SP/EF/EM target-junction summary with matched REF/SF3A2 sequence comparison and predicted frame effects

## Repository structure

```text
sf3a2_transcriptomics_github/
├── config/
│   ├── config_template.yml
│   └── longread_stage_bam_map_template.csv
├── R/
│   ├── helper_functions.R
│   └── longread_helpers.R
├── short_reads/
│   └── scripts/
│       ├── 00_setup.R
│       ├── 01_qc_doublet_filtering.R
│       ├── 02_integrate_nf54_sf3a2.R
│       ├── 03_mca_label_transfer.R
│       ├── 04_cluster_annotation_and_markers.R
│       ├── 05_gametocyte_subset_analysis.R
│       ├── 06_representative_gene_expression.R
│       ├── 07a_slingshot_tradeseq_asexual.R
│       ├── 07b_slingshot_tradeseq_gametocyte.R
│       ├── 08a_pseudotime_heatmap_asexual.R
│       └── 08b_pseudotime_heatmap_gametocyte.R
├── long_reads/
│   ├── scripts/
│   │   ├── 01_add_stage_to_sqanti_classification.py
│   │   ├── 02_dexseq_pseudobulk_exon_usage.R
│   │   ├── 03_sqanti_junction_comparison.R
│   │   ├── 04_sqanti_junction_comparison_strict_novelty.R
│   │   ├── 05_prepare_cis_element_inputs.R
│   │   ├── 06_summarize_motif_results.R
│   │   ├── 07_plot_branchpoint_acceptor_logos.R
│   │   ├── 08_target_gene_cis_context.R
│   │   └── 09_sf3a2_strict_novel_sp_ef_em_summary.R
│   └── slurm/
│       └── run_streme_fimo_cis_elements.sh
├── data/
├── results/
└── docs/
```

## Before running

Copy the config template:

```bash
cp config/config_template.yml config/config.yml
```

Edit paths in `config/config.yml` so they match your local workstation or HPC environment.

## Short-read workflow

Run these in order:

```bash
Rscript short_reads/scripts/01_qc_doublet_filtering.R
Rscript short_reads/scripts/02_integrate_nf54_sf3a2.R
Rscript short_reads/scripts/03_mca_label_transfer.R
Rscript short_reads/scripts/04_cluster_annotation_and_markers.R
Rscript short_reads/scripts/05_gametocyte_subset_analysis.R
Rscript short_reads/scripts/06_representative_gene_expression.R
```

Run pseudotime separately:

```bash
Rscript short_reads/scripts/07a_slingshot_tradeseq_asexual.R
Rscript short_reads/scripts/08a_pseudotime_heatmap_asexual.R

Rscript short_reads/scripts/07b_slingshot_tradeseq_gametocyte.R
Rscript short_reads/scripts/08b_pseudotime_heatmap_gametocyte.R
```

## Long-read workflow

Run these as needed:

```bash
python long_reads/scripts/01_add_stage_to_sqanti_classification.py \
  --classification data/long_reads/sqanti3/SF3A2_classification.txt \
  --abundance data/long_reads/sqanti3/SF3A2_mapped.abundance.txt \
  --barcode-stage data/long_reads/sqanti3/SF3A2_barcode_to_stage.txt \
  --output results/long_reads/tables/SF3A2_classification_with_stage.txt

Rscript long_reads/scripts/02_dexseq_pseudobulk_exon_usage.R Progenitor
Rscript long_reads/scripts/03_sqanti_junction_comparison.R
Rscript long_reads/scripts/04_sqanti_junction_comparison_strict_novelty.R
Rscript long_reads/scripts/05_prepare_cis_element_inputs.R
bash long_reads/slurm/run_streme_fimo_cis_elements.sh
Rscript long_reads/scripts/06_summarize_motif_results.R
Rscript long_reads/scripts/07_plot_branchpoint_acceptor_logos.R
Rscript long_reads/scripts/08_target_gene_cis_context.R
Rscript long_reads/scripts/09_sf3a2_strict_novel_sp_ef_em_summary.R
```

## GitHub data policy

Do **not** commit:

- FASTQ files
- BAM/BAI files
- large `.rds` objects
- SQANTI3 full raw outputs if unpublished/private
- Cell Ranger matrices
- personal local paths

Use `config/config.yml` locally and keep it out of version control.
