#!/usr/bin/env Rscript
# =========================================================================
# Script: atac_te_diff_analysis_in_r.R
# Author: Alireza Ghahramani
# Contact: aghahram@uwo.ca
#
# Purpose:
#   Quantify transposable elements (TEs) in ATAC-seq BAMs using featureCounts
#   (with a RepeatMasker SAF annotation), and perform differential analysis
#   with DESeq2.
#
# Inputs:
#   - Directory with sample BAM files (deduped, sorted, indexed)
#   - Sample metadata CSV (columns: sample, bam, condition, etc.)
#   - RepeatMasker annotation in SAF format
#
# Outputs:
#   - TE count matrix (for record-keeping)
#   - DESeq2 results (CSV)
#   - Volcano plot (PNG)
#   - Heatmap of top variable TEs (PNG)
#   - Session info (TXT)
#
# Requirements:
#   - R packages: Rsubread, DESeq2, ggplot2, pheatmap
#
# =========================================================================

# -------- USER PARAMETERS (EDIT THESE) --------
bam_dir <- "aligned"                              # Directory with BAM files
metadata_file <- "sample_metadata.csv"            # Sample metadata (must include 'sample' and 'bam' columns)
saf_file <- "mm10_repeatmasker_for_featurecounts.txt"  # RepeatMasker SAF annotation file
output_dir <- "atac_te_deseq2_results"            # Output directory
condition_column <- "condition"                   # Metadata column for DE analysis

dir.create(output_dir, showWarnings = FALSE)

# -------- 1. LOAD LIBRARIES --------
suppressPackageStartupMessages({
  library(Rsubread)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

# -------- 2. LOAD SAMPLE METADATA --------
metadata <- read.csv(metadata_file, stringsAsFactors=FALSE)
# Ensure BAM path is correct for each sample
bam_files <- file.path(bam_dir, metadata$bam)    # metadata$bam: BAM filenames

# -------- 3. QUANTIFY TEs USING featureCounts --------
fc <- featureCounts(
  files = bam_files,
  annot.ext = saf_file,
  isPairedEnd = TRUE,
  nthreads = 8
)
counts <- fc$counts
colnames(counts) <- metadata$sample
rownames(counts) <- rownames(fc$annotation)

# Save the TE count matrix (for records)
write.csv(counts, file=file.path(output_dir, "TE_counts_matrix.csv"))

# -------- 4. DESEQ2 ANALYSIS --------
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

write.csv(res_df, file=file.path(output_dir, "te_deseq2_results.csv"))

# -------- 5. PLOTS --------
# Volcano plot
ggplot(res_df, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color=padj < 0.05), alpha=0.7) +
  theme_minimal() +
  labs(title="ATAC-seq TE Volcano Plot", x="log2 Fold Change", y="-log10 adj p-value") +
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

cat("ATAC-seq TE differential analysis complete. Results in", output_dir, "\n")
