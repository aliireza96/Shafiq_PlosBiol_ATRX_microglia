#!/usr/bin/env Rscript
# ==============================================================================
# Script: rnaseq_TEtranscripts_differential_analysis.R
# Author: Alireza Ghahramani
# Contact: aghahram@uwo.ca
#
# Methods Overview:
#   - TE expression was quantified with **TEtranscripts** (Jin et al., 2015) in
#     `--mode multi`, which probabilistically redistributes multi-mapping reads across
#     TE subfamilies rather than discarding them.
#   - TE annotation: the **curated rmsk TE GTF** distributed by the TEtranscripts authors
#     (GRCm38_Ensembl_rmsk_TE.gtf). It excludes simple repeats and non-TE RNA species
#     (rRNA, scRNA, snRNA, srpRNA, tRNA), so the analysis is restricted to TEs.
#   - **Libraries are reverse-stranded (dUTP).** TEtranscripts must therefore be run with
#     `--stranded reverse`. Running it with `--stranded no` on these libraries counts
#     antisense and intronic reads toward TEs and inflates TE signal.
#   - **Alignment:** BAMs come from the pipeline in `rna_seq_preprocessing_pipeline.sh`
#     (Trim Galore v0.6.6, HISAT2 v2.0.4, GRCm38.p6), as described in the manuscript
#     Methods. Note that `--mode multi` is only meaningful if the aligner actually
#     reported multi-mapping alignments; if you re-run this on your own data, check your
#     aligner's multimapper settings before trusting subfamily-level counts for young,
#     high-copy elements (L1Md, IAPEz, MMERVK10C).
#
# TEtranscripts command used (edit paths and sample groupings):
#
#   TEtranscripts --format BAM --mode multi \
#       -t <KO1.bam> <KO2.bam> <KO3.bam> \
#       -c <CTRL1.bam> <CTRL2.bam> <CTRL3.bam> \
#       --GTF <gene_annotation.gtf> \
#       --TE GRCm38_Ensembl_rmsk_TE.gtf \
#       --project <output_prefix> \
#       --stranded reverse \
#       --sortByPos
#
# Outputs written by TEtranscripts:
#   - <output_prefix>.cntTable                     raw counts, GENES AND TEs TOGETHER
#   - <output_prefix>_gene_TE_analysis.txt         its own DESeq2 result
#   - <output_prefix>_sigdiff_gene_TE.txt          significant subset
#   - <output_prefix>_DESeq2.R                     the DESeq2 script it generates
#
# This R script is an alternative to the DESeq2 script TEtranscripts emits. It performs:
#   - Differential expression of TEs using DESeq2, with an explicit contrast
#   - Volcano plot and heatmap visualization
#   - sessionInfo export for reproducibility
#
# Inputs for R script:
#   - <output_prefix>.cntTable from TEtranscripts
#   - Sample metadata (CSV with a `sample` column matching the cntTable column names)
#
# Requirements:
#   - R packages: DESeq2, ggplot2, pheatmap
#
# Usage:
#   1. Run TEtranscripts as above.
#   2. Edit USER PARAMETERS below.
#   3. Run: Rscript rnaseq_TEtranscripts_differential_analysis.R
# ==============================================================================

# -------- USER PARAMETERS (EDIT THESE) --------
te_counts_file <- "project.cntTable"           # TEtranscripts output (genes + TEs in one table)
metadata_file <- "sample_metadata.csv"         # Path to sample metadata (must match columns to counts)
output_dir <- "rnaseq_TEtranscripts_results"   # Output directory for results
condition_column <- "condition"                # Metadata column for group comparison
reference_level <- "control"                   # baseline level of condition_column (denominator)
treatment_level <- "KO"                        # level compared against the baseline (numerator)

dir.create(output_dir, showWarnings = FALSE)

# -------- 1. LOAD LIBRARIES --------
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

# -------- 2. LOAD DATA --------
# The cntTable holds genes AND TEs. Both are kept in the DESeq2 object so that size
# factors are estimated from the whole library; TE rows are separated afterwards, at the
# results stage. Dropping genes here would normalise TEs against TEs only.
counts_all <- read.table(te_counts_file, header=TRUE, row.names=1, check.names=FALSE)
counts_all <- counts_all[!grepl("^__", rownames(counts_all)), ] # htseq-style summary rows

# TEtranscripts names TE rows "subfamily:family:class" and gene rows by gene ID.
is_te <- grepl(":", rownames(counts_all), fixed = TRUE)
message(sprintf("Loaded %d rows: %d TE, %d gene", nrow(counts_all), sum(is_te), sum(!is_te)))

# Load sample metadata
metadata <- read.csv(metadata_file, stringsAsFactors=FALSE)
metadata <- metadata[match(colnames(counts_all), metadata$sample), ]
rownames(metadata) <- metadata$sample
stopifnot(!any(is.na(metadata$sample)))

# -------- 3. DESEQ2 ANALYSIS --------
dds <- DESeqDataSetFromMatrix(
  countData = counts_all,
  colData = metadata,
  design = as.formula(paste("~", condition_column))
)
dds <- dds[rowSums(counts(dds)) > 10, ]    # Filter low-count TEs

# Set the baseline explicitly. Without this DESeq2 orders levels alphabetically,
# which silently decides the sign of every log2 fold change.
dds[[condition_column]] <- relevel(factor(dds[[condition_column]]), ref = reference_level)

dds <- DESeq(dds)
res <- results(dds, contrast=c(condition_column, treatment_level, reference_level), alpha=0.05)
res_all <- as.data.frame(res)
res_all$feature_id <- rownames(res_all)

# Now restrict to TEs. Everything above was fitted on genes + TEs together.
res_df <- res_all[grepl(":", rownames(res_all), fixed = TRUE), ]
res_df$TE_id <- rownames(res_df)
# "subfamily:family:class" -> subfamily, for plot labels
res_df$TE_subfamily <- sub(":.*", "", res_df$TE_id)

write.csv(res_all, file=file.path(output_dir, "gene_and_TE_deseq2_results.csv"))
write.csv(res_df,  file=file.path(output_dir, "TE_deseq2_results.csv"))

# -------- 4. PLOTS --------
# Volcano plot
# padj is NA wherever DESeq2 independent filtering removed the feature; drop those
# rows explicitly so they are not silently discarded by the plot.
res_df_plot <- res_df[!is.na(res_df$padj), ]
p_te_volcano_plot <- ggplot(res_df_plot, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha=0.7) +
  theme_minimal() +
  labs(title="RNA-seq TE Volcano Plot (TEtranscripts)", x="log2 Fold Change", y="-log10 adj p-value") +
  scale_color_manual(values=c("FALSE"="gray", "TRUE"="red"), name="padj < 0.05")
ggsave(file.path(output_dir, "te_volcano_plot.png"), p_te_volcano_plot, width=6, height=5)

# Heatmap of the 30 most significant TEs (by padj, not by variance)
# Rank TEs by padj and index the VST matrix BY NAME. res_df is a subset of the full
# fit, so positional indexing into assay(vsd) would pull the wrong rows.
top_tes <- head(res_df$TE_id[order(res_df$padj)], 30)
if (length(top_tes) >= 2) {
  vsd <- varianceStabilizingTransformation(dds)
  mat <- assay(vsd)[top_tes, ]
  pheatmap(mat, scale="row", annotation_col=metadata[, condition_column, drop=FALSE], filename=file.path(output_dir, "top_tes_heatmap.png"))
}

# Session info
sink(file.path(output_dir, "sessionInfo.txt"))
sessionInfo()
sink()

cat("RNA-seq TE differential analysis (TEtranscripts) complete. Results in", output_dir, "\n")
