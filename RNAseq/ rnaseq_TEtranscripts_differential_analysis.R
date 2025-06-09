#!/usr/bin/env Rscript
# ==============================================================================
# Script: rnaseq_TEtranscripts_differential_analysis.R
# Author: Alireza Ghahramani
# Contact: aghahram@uwo.ca
#
# Methods Overview:
#   - TE expression analysis was performed using **TEtranscripts** (v2.0.3) (Jin et al., 2015),
#     designed for accurate TE quantification.
#   - The **curated GTF file** from the TEtranscripts authors was used as the TE annotation.
#     This file excludes simple repeats and non-TE RNA species (rRNAs, scRNAs, snRNAs,
#     srpRNAs, tRNAs), ensuring analysis focused on TEs only.
#   - **Read alignment:** RNA-seq reads were aligned to GRCm38/mm10 using STAR (v2.7.3a)
#     with parameters optimized for TE quantification:
#        --winAnchorMultimapNmax 100   (allow up to 100 anchors for multimappers)
#        --outFilterMultimapNmax 100   (allow each read to map to up to 100 locations)
#   - **TEtranscripts was run in --mode multi** to appropriately assign multi-mapping reads,
#     as recommended for repetitive elements.
#
# Shell commands for alignment and TEtranscripts quantification:
#
#   # STAR alignment (example per sample; edit as needed)
#   STAR --runThreadN 8 \
#        --genomeDir /path/to/STAR_index \
#        --readFilesIn <sample_R1.fastq.gz> <sample_R2.fastq.gz> \
#        --readFilesCommand zcat \
#        --outSAMtype BAM SortedByCoordinate \
#        --winAnchorMultimapNmax 100 \
#        --outFilterMultimapNmax 100 \
#        --outFileNamePrefix <sample>
#
#   # TEtranscripts quantification (edit for your sample groupings and paths)
#   TEtranscripts --format BAM \
#       -t <KO1.bam,KO2.bam,...> \
#       -c <CTRL1.bam,CTRL2.bam,...> \
#       --GTF <gene_annotation.gtf> \
#       --TE <curated_TE_annotation.gtf> \
#       --project <output_prefix> \
#       --mode multi \
#       --stranded no \
#       --minl 50
#
# Outputs:
#   - <output_prefix>_te.counts: raw TE count matrix for all samples
#
# This R script performs:
#   - Differential expression analysis of TEs using DESeq2
#   - Volcano plot and heatmap visualization
#   - sessionInfo export for reproducibility
#
# Inputs for R script:
#   - TE count table from TEtranscripts (<output_prefix>_te.counts)
#   - Sample metadata (CSV: sample, condition, etc.)
#
# Requirements:
#   - R packages: DESeq2, ggplot2, pheatmap
#
# Usage:
#   1. Run STAR and TEtranscripts as described above.
#   2. Edit USER PARAMETERS below.
#   3. Run: source("rnaseq_TEtranscripts_differential_analysis.R")
# ==============================================================================

# -------- USER PARAMETERS (EDIT THESE) --------
te_counts_file <- "project_te.counts"          # Path to TEtranscripts counts output
metadata_file <- "sample_metadata.csv"         # Path to sample metadata (must match columns to counts)
output_dir <- "rnaseq_TEtranscripts_results"   # Output directory for results
condition_column <- "condition"                # Metadata column for group comparison

dir.create(output_dir, showWarnings = FALSE)

# -------- 1. LOAD LIBRARIES --------
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

# -------- 2. LOAD DATA --------
# Load TEtranscripts counts
te_counts <- read.table(te_counts_file, header=TRUE, row.names=1, check.names=FALSE)
te_counts <- te_counts[!grepl("^__", rownames(te_counts)), ] # Remove summary rows if present

# Load sample metadata
metadata <- read.csv(metadata_file, stringsAsFactors=FALSE)
metadata <- metadata[match(colnames(te_counts), metadata$sample), ]
rownames(metadata) <- metadata$sample

# -------- 3. DESEQ2 ANALYSIS --------
dds <- DESeqDataSetFromMatrix(
  countData = te_counts,
  colData = metadata,
  design = as.formula(paste("~", condition_column))
)
dds <- dds[rowSums(counts(dds)) > 10, ]    # Filter low-count TEs

dds <- DESeq(dds)
res <- results(dds, alpha=0.05)
res_df <- as.data.frame(res)
res_df$TE_id <- rownames(res_df)

write.csv(res_df, file=file.path(output_dir, "TE_deseq2_results.csv"))

# -------- 4. PLOTS --------
# Volcano plot
ggplot(res_df, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color=padj < 0.05), alpha=0.7) +
  theme_minimal() +
  labs(title="RNA-seq TE Volcano Plot (TEtranscripts)", x="log2 Fold Change", y="-log10 adj p-value") +
  scale_color_manual(values=c("gray", "red")) +
  ggsave(file.path(output_dir, "te_volcano_plot.png"), width=6, height=5)

# Heatmap of top 30 variable TEs
top_tes <- head(order(res_df$padj), 30)
if (length(top_tes) >= 2) {
  vsd <- varianceStabilizingTransformation(dds)
  mat <- assay(vsd)[top_tes, ]
  pheatmap(mat, scale="row", annotation_col=metadata, filename=file.path(output_dir, "top_tes_heatmap.png"))
}

# Session info
sink(file.path(output_dir, "sessionInfo.txt"))
sessionInfo()
sink()

cat("RNA-seq TE differential analysis (TEtranscripts) complete. Results in", output_dir, "\n")
