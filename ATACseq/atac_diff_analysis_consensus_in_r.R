#!/usr/bin/env Rscript
# =========================================================================
# Script: atac_diff_analysis_consensus_in_r.R
# Author: Alireza Ghahramani
# Contact: aghahram@uwo.ca
#
# Purpose:
#   1. Load all sample MACS2 narrowPeak files
#   2. Merge/union them to create a consensus peak set (BED)
#   3. Quantify reads in all peaks for all samples (featureCounts in R)
#   4. Differential analysis of peaks using DESeq2
#
# Inputs:
#   - Directory with all sample narrowPeak files (*.narrowPeak)
#   - Directory with all sample BAM files (deduplicated, sorted, indexed)
#   - Sample metadata CSV (with sample names, condition, and filenames)
#
#
# Requirements:
#   - R packages: GenomicRanges, rtracklayer, Rsubread, DESeq2, ggplot2, pheatmap
#
# =========================================================================

# -------- USER PARAMETERS (EDIT THESE) --------
narrowpeak_dir <- "peaks"                 # Directory containing *.narrowPeak files
bam_dir <- "aligned"                      # Directory containing BAM files
metadata_file <- "sample_metadata.csv"    # Path to sample metadata CSV
output_dir <- "atac_deseq2_results"      # Output directory
condition_column <- "condition"           # Metadata column for comparison

dir.create(output_dir, showWarnings = FALSE)

# -------- 1. LOAD LIBRARIES --------
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(Rsubread)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

# -------- 2. LOAD SAMPLE METADATA --------
metadata <- read.csv(metadata_file, stringsAsFactors=FALSE)

# -------- 3. LOAD AND MERGE ALL narrowPeak FILES INTO CONSENSUS PEAKS --------
peak_files <- list.files(narrowpeak_dir, pattern="*.narrowPeak$", full.names=TRUE)
peak_list <- lapply(peak_files, function(f) {
  gr <- import(f, format="narrowPeak")
  seqlevelsStyle(gr) <- "UCSC"
  gr
})

# Merge all peaks to consensus (union)
all_peaks <- suppressWarnings(reduce(do.call(c, peak_list)))
# Optionally, filter for width > X, or for peak score

# Save consensus as BED (for records)
consensus_bed <- file.path(output_dir, "consensus_peaks.bed")
export(all_peaks, consensus_bed, format="BED")

cat("Consensus peaks:", length(all_peaks), "regions\n")

# -------- 4. QUANTIFY READS IN CONSENSUS PEAKS WITH featureCounts --------
# Prepare SAF annotation from GRanges for Rsubread::featureCounts
saf <- data.frame(
  GeneID = paste(seqnames(all_peaks), start(all_peaks), end(all_peaks), sep="_"),
  Chr = as.character(seqnames(all_peaks)),
  Start = start(all_peaks),
  End = end(all_peaks),
  Strand = "*"
)
saf_file <- file.path(output_dir, "consensus_peaks.saf")
write.table(saf, saf_file, sep="\t", row.names=FALSE, quote=FALSE)

# Get BAM files in order matching metadata
bam_files <- file.path(bam_dir, metadata$bam)   # metadata$bam should list BAM filenames

fc <- featureCounts(
  files = bam_files,
  annot.ext = saf_file,
  isPairedEnd = TRUE,
  nthreads = 8
)
peak_counts <- fc$counts
colnames(peak_counts) <- metadata$sample
rownames(peak_counts) <- saf$GeneID

# -------- 5. DESEQ2 ANALYSIS --------
dds <- DESeqDataSetFromMatrix(
  countData = peak_counts,
  colData = metadata,
  design = as.formula(paste("~", condition_column))
)
dds <- dds[rowSums(counts(dds)) > 10, ]    # Filter low-count peaks

dds <- DESeq(dds)
res <- results(dds, alpha=0.05)
res_df <- as.data.frame(res)
res_df$peak_id <- rownames(res_df)

write.csv(res_df, file=file.path(output_dir, "consensus_peak_deseq2_results.csv"))

# -------- 6. PLOTS --------
ggplot(res_df, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color=padj < 0.05), alpha=0.7) +
  theme_minimal() +
  labs(title="Consensus Peak Volcano Plot", x="log2 Fold Change", y="-log10 adj p-value") +
  scale_color_manual(values=c("gray", "red")) +
  ggsave(file.path(output_dir, "consensus_peak_volcano_plot.png"), width=6, height=5)

top_peaks <- head(order(res_df$padj), 30)
if (length(top_peaks) >= 2) {
  vsd <- varianceStabilizingTransformation(dds)
  mat <- assay(vsd)[top_peaks, ]
  pheatmap(mat, scale="row", annotation_col=metadata, filename=file.path(output_dir, "top_consensus_peaks_heatmap.png"))
}

sink(file.path(output_dir, "sessionInfo.txt"))
sessionInfo()
sink()

cat("ATAC-seq consensus peak differential analysis complete. Results in", output_dir, "\n")
