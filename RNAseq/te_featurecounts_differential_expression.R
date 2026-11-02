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
#   This is the per-locus arm of the TE analysis. It complements the subfamily-level
#   TEtranscripts arm in `rnaseq_TEtranscripts_differential_analysis.R`; the manuscript
#   reports a TE result as supported when both arms agree.
#
#   The libraries are reverse-stranded (dUTP), and TEs are frequently transcribed
#   antisense to, or read through from, an overlapping host gene. Counting both strands
#   together therefore conflates autonomous TE transcription with host transcription. The
#   published workflow avoids this by splitting each BAM by transcriptional strand first
#   and testing the two strands separately:
#
#   1. Split each BAM into sense and antisense with `split_bam_and_bigwig.sh`.
#   2. Count reads over the UCSC RepeatMasker annotation for each strand separately,
#      using featureCounts (Rsubread) with a SAF whose GeneID encodes the locus:
#         repName#repFamily#repClass#chr#start#end
#   3. Run this script once on the sense matrix and once on the antisense matrix.
#   4. Intersect the significant loci from the two runs (GenomicRanges::intersect) and
#      carry forward the elements recovered on both.
#
#   Step 2 was run twice, counting multi-mapping reads and excluding them
#   (`countMultiMappingReads = TRUE` / `FALSE` in Rsubread::featureCounts). The reported
#   LTR result is present in both; report which variant a given figure used.
#
# Equivalent featureCounts command line for one strand:
#   featureCounts -p -F SAF -T 8 \
#     -a <mm10_repeatmasker.saf> \
#     -o <output_counts.txt> \
#     <sample1.fwd.bam> <sample2.fwd.bam> ...
#
#   Do not pass `-s`: the BAMs have already been split by strand, so each split BAM is
#   counted unstranded. Note also that `--ignoreDup` was NOT used; RNA-seq duplicates are
#   largely expression-driven and dropping them biases high-expression loci.
#
# Inputs:
#   - featureCounts output for ONE strand (SAF-based; rename per strand as needed)
#   - Sample metadata CSV (columns: sample, condition, [optional covariates])
#
# ========================================================================

# -------- USER PARAMETERS (EDIT THESE) --------
te_count_file <- "featurecounts_TE_counts_fwd.txt" # featureCounts output for ONE strand (see header)
metadata_file <- "sample_metadata.csv"             # Path to sample metadata file
output_dir <- "te_deseq2_results_fwd"              # Output directory for results (one per strand)
condition_column <- "condition"                    # Column in metadata with condition/group info
reference_level <- "control"                   # baseline level of condition_column (denominator)
treatment_level <- "KO"                        # level compared against the baseline (numerator)

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

# Set the baseline explicitly. Without this DESeq2 orders levels alphabetically,
# which silently decides the sign of every log2 fold change.
dds[[condition_column]] <- relevel(factor(dds[[condition_column]]), ref = reference_level)

dds <- DESeq(dds)
res <- results(dds, contrast=c(condition_column, treatment_level, reference_level), alpha=0.05)
res_df <- as.data.frame(res)
res_df$TE_id <- rownames(res_df)

# -------- 3. EXPORT RESULTS --------
# Split the SAF GeneID back into its parts so the result can be turned into a GRanges
# and intersected with the other strand.
id_parts <- do.call(rbind, strsplit(rownames(res_df), "#", fixed=TRUE))
if (ncol(id_parts) == 6) {
  colnames(id_parts) <- c("repName", "repFamily", "repClass", "chr", "start", "end")
  res_df <- cbind(res_df, as.data.frame(id_parts, stringsAsFactors=FALSE))
  res_df$start <- as.numeric(res_df$start)
  res_df$end   <- as.numeric(res_df$end)
} else {
  warning("GeneID is not repName#repFamily#repClass#chr#start#end; locus columns not added.")
}

write.csv(res_df, file=file.path(output_dir, "TE_deseq2_results.csv"))

# -------- 4. PLOTS --------
# Volcano plot
# padj is NA wherever DESeq2 independent filtering removed the feature; drop those
# rows explicitly so they are not silently discarded by the plot.
res_df_plot <- res_df[!is.na(res_df$padj), ]
p_te_volcano_plot <- ggplot(res_df_plot, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha=0.7) +
  theme_minimal() +
  labs(title="TE Volcano Plot", x="log2 Fold Change", y="-log10 adj p-value") +
  scale_color_manual(values=c("FALSE"="gray", "TRUE"="red"), name="padj < 0.05")
ggsave(file.path(output_dir, "TE_volcano_plot.png"), p_te_volcano_plot, width=6, height=5)

# Heatmap of top 30 TEs
top_tes <- head(order(res_df$padj), 30)
if (length(top_tes) >= 2) {
  vsd <- varianceStabilizingTransformation(dds)
  mat <- assay(vsd)[top_tes, ]
  pheatmap(mat, scale="row", annotation_col=metadata[, condition_column, drop=FALSE], filename=file.path(output_dir, "top_TEs_heatmap.png"))
}

# -------- 5. SESSION INFO --------
sink(file.path(output_dir, "sessionInfo.txt"))
sessionInfo()
sink()

cat("TE differential expression analysis complete. Results in", output_dir, "\n")
