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
reference_level <- "control"                   # baseline level of condition_column (denominator)
treatment_level <- "KO"                        # level compared against the baseline (numerator)

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

# Set the baseline explicitly. Without this DESeq2 orders levels alphabetically,
# which silently decides the sign of every log2 fold change.
dds[[condition_column]] <- relevel(factor(dds[[condition_column]]), ref = reference_level)

dds <- DESeq(dds)
res <- results(dds, contrast=c(condition_column, treatment_level, reference_level), alpha=0.05)
res_df <- as.data.frame(res)
res_df$peak_id <- rownames(res_df)

write.csv(res_df, file=file.path(output_dir, "consensus_peak_deseq2_results.csv"))

# -------- 6. PLOTS --------
# padj is NA wherever DESeq2 independent filtering removed the feature; drop those
# rows explicitly so they are not silently discarded by the plot.
res_df_plot <- res_df[!is.na(res_df$padj), ]
p_consensus_peak_volcano_plot <- ggplot(res_df_plot, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha=0.7) +
  theme_minimal() +
  labs(title="Consensus Peak Volcano Plot", x="log2 Fold Change", y="-log10 adj p-value") +
  scale_color_manual(values=c("FALSE"="gray", "TRUE"="red"), name="padj < 0.05")
ggsave(file.path(output_dir, "consensus_peak_volcano_plot.png"), p_consensus_peak_volcano_plot, width=6, height=5)

top_peaks <- head(order(res_df$padj), 30)
if (length(top_peaks) >= 2) {
  vsd <- varianceStabilizingTransformation(dds)
  mat <- assay(vsd)[top_peaks, ]
  pheatmap(mat, scale="row", annotation_col=metadata[, condition_column, drop=FALSE], filename=file.path(output_dir, "top_consensus_peaks_heatmap.png"))
}

sink(file.path(output_dir, "sessionInfo.txt"))
sessionInfo()
sink()

cat("ATAC-seq consensus peak differential analysis complete. Results in", output_dir, "\n")
