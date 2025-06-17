#!/usr/bin/env Rscript
# ========================================================================
# Script: te_featurecounts_differential_expression.R
# Author: Alireza Ghahramani
# Contact: aghahram@uwo.ca
#
# Purpose:
#   Differential expression analysis of transposable elements (TEs) using
#   featureCounts (SAF, RepeatMasker annotation) and DESeq2 in R.
#
# Overview:
#   1. Quantify TE counts from BAM files using featureCounts and a
#      RepeatMasker SAF annotation file (see command below).
#   2. Run this R script to perform differential expression analysis of TEs.
#
# Example featureCounts command (run outside of R):
#   featureCounts -p -F SAF --ignoreDup -T 8 \
#     -a <repeatmasker_annotation.saf> \
#     -o <output_counts.txt> \
#     <bam_files>
#
# Inputs:
#   - featureCounts TE output (SAF format; can be renamed as needed)
#   - Sample metadata CSV (columns: sample, condition, [optional covariates])
#
# ========================================================================

# -------- USER PARAMETERS (EDIT THESE) --------
te_count_file <- "featurecounts_TE_counts.txt"     # Path to featureCounts output (TEs, SAF format)
metadata_file <- "sample_metadata.csv"             # Path to sample metadata file
output_dir <- "te_deseq2_results"                  # Output directory for results
condition_column <- "condition"                    # Column in metadata with condition/group info

dir.create(output_dir, showWarnings = FALSE)

# -------- 1. LOAD DATA --------
# Load featureCounts TE matrix
counts_raw <- read.table(te_count_file, header=TRUE, sep="\t", comment.char="#", check.names=FALSE)
counts <- counts_raw[ , 7:ncol(counts_raw)]
rownames(counts) <- counts_raw$Geneid

# Remove summary rows (e.g., "__no_feature") if present
counts <- counts[!grepl("^__", rownames(counts)), ]

# Load sample metadata
metadata <- read.csv(metadata_file, stringsAsFactors=FALSE)
# Ensure order matches columns in counts
metadata <- metadata[match(colnames(counts), metadata$sample), ]
rownames(metadata) <- metadata$sample

# -------- 2. DESEQ2 ANALYSIS --------
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = as.formula(paste("~", condition_column))
)
dds <- dds[rowSums(counts(dds)) > 10, ]    # Filter low-count TEs

dds <- DESeq(dds)
res <- results(dds, alpha=0.05)
res_df <- as.data.frame(res)
res_df$TE_id <- rownames(res_df)

# -------- 3. EXPORT RESULTS --------
write.csv(res_df, file=file.path(output_dir, "TE_deseq2_results.csv"))

# -------- 4. PLOTS --------
# Volcano plot
ggplot(res_df, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color=padj < 0.05), alpha=0.7) +
  theme_minimal() +
  labs(title="TE Volcano Plot", x="log2 Fold Change", y="-log10 adj p-value") +
  scale_color_manual(values=c("gray", "red")) +
  ggsave(file.path(output_dir, "TE_volcano_plot.png"), width=6, height=5)

# Heatmap of top 30 TEs
top_tes <- head(order(res_df$padj), 30)
if (length(top_tes) >= 2) {
  vsd <- varianceStabilizingTransformation(dds)
  mat <- assay(vsd)[top_tes, ]
  pheatmap(mat, scale="row", annotation_col=metadata, filename=file.path(output_dir, "top_TEs_heatmap.png"))
}

# -------- 5. SESSION INFO --------
sink(file.path(output_dir, "sessionInfo.txt"))
sessionInfo()
sink()

cat("TE differential expression analysis complete. Results in", output_dir, "\n")
