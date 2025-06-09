#!/usr/bin/env Rscript
# ==============================================================================
# Script: rnaseq_transcript_and_gene_level_differential_analysis.R
# Author: Alireza Ghahramani
# Contact: aghahram@uwo.ca
#
# Purpose:
#   - Transcript-level RNA-seq DESeq2 analysis (StringTie quantification)
#   - Gene-level differential analysis by aggregating transcript p-values
#     with the Lancaster method, and summarizing log2FC per gene
#
# Inputs:
#   - Transcript abundance files (StringTie output, e.g., *.abund.tab)
#   - Sample metadata (CSV: sample, condition, etc.)
#   - GTF annotation file (for transcript→gene mapping)
#
# Requirements:
#   - R packages: tximport, DESeq2, aggregation, GenomicFeatures, clusterProfiler,
#     org.Mm.eg.db, ggplot2, pheatmap
# ==============================================================================

# -------- USER PARAMETERS (EDIT THESE) --------
abundance_dir <- "analysis_output/stringtie"    # Folder with StringTie .abund.tab files
metadata_file <- "sample_metadata.csv"          # CSV file: sample, condition, etc.
gtf_file <- "/path/to/annotation.gtf"           # GTF file for tx2gene mapping
output_dir <- "rnaseq_results"
condition_column <- "condition"                 # Column name in metadata for group/condition

dir.create(output_dir, showWarnings = FALSE)

# -------- 1. LOAD LIBRARIES --------
suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(GenomicFeatures)
  library(aggregation)
  library(ggplot2)
  library(pheatmap)
})

# -------- 2. PREPARE FILES AND METADATA --------
metadata <- read.csv(metadata_file, stringsAsFactors=FALSE)
rownames(metadata) <- metadata$sample

# Build tx2gene mapping from GTF
txdb <- makeTxDbFromGFF(gtf_file)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")

# Collect all StringTie abundance files for the samples in metadata
files <- file.path(abundance_dir, paste0(metadata$sample, ".abund.tab"))
names(files) <- metadata$sample

# -------- 3. TXIMPORT: READ TRANSCRIPT-LEVEL QUANTIFICATIONS --------
txi <- tximport(
  files,
  type = "stringtie",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE,
  txOut = TRUE           # <-- key for transcript-level
)

# -------- 4. TRANSCRIPT-LEVEL DESEQ2 ANALYSIS --------
dds <- DESeqDataSetFromTximport(
  txi,
  colData = metadata,
  design = as.formula(paste("~", condition_column)),
  txOut = TRUE
)
dds <- dds[rowSums(counts(dds)) > 10, ]  # Filter out low-expression transcripts

dds <- DESeq(dds)
res <- results(dds, alpha=0.05)
res_df <- as.data.frame(res)
res_df$transcript_id <- rownames(res_df)

# Map transcripts to genes using tx2gene mapping
res_df$gene_id <- tx2gene$GENEID[match(res_df$transcript_id, tx2gene$TXNAME)]

# Export transcript-level results
write.csv(res_df, file=file.path(output_dir, "transcript_level_deseq2_results.csv"), row.names=FALSE)

# -------- 5. LANCASTER AGGREGATION TO GENE LEVEL --------
# Prepare gene p-value list
keep <- !is.na(res_df$gene_id) & !is.na(res_df$pvalue)
gene_pvals <- split(res_df$pvalue[keep], res_df$gene_id[keep])

# Lancaster method (from "aggregation" package)
agg_results <- aggregate_pvalues(gene_pvals, method = "lancaster")
agg_results$padj <- p.adjust(agg_results$pvalue, method = "BH")

# Gene-level log2FC summary: Use log2FC of transcript with smallest p-value for each gene (as in your script)
log2fc_summary <- tapply(
  X = res_df$log2FoldChange[keep],
  INDEX = res_df$gene_id[keep],
  FUN = function(fc) {
    idx <- which.min(res_df$pvalue[keep][res_df$gene_id[keep] == names(fc)[1]])
    fc[idx]
  }
)
gene_ids <- names(log2fc_summary)
gene_log2fc_df <- data.frame(gene_id = gene_ids, log2FoldChange = as.numeric(log2fc_summary), stringsAsFactors = FALSE)

# Merge Lancaster p-values and log2FC summary
final_gene_results <- merge(agg_results, gene_log2fc_df, by = "gene_id")

# Export gene-level results
write.csv(final_gene_results, file=file.path(output_dir, "gene_level_DE_lancaster.csv"), row.names=FALSE)

# -------- 6. PLOTS --------
# Volcano plot: transcript level
ggplot(res_df, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha=0.7) +
  theme_minimal() +
  labs(title="Transcript-level Volcano plot", x="log2 Fold Change", y="-log10 adj p-value") +
  scale_color_manual(values=c("gray", "red")) +
  ggsave(file.path(output_dir, "volcano_plot_transcript.png"), width=6, height=5)

# Volcano plot: gene level
ggplot(final_gene_results, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha=0.7) +
  theme_minimal() +
  labs(title="Gene-level Volcano plot (Lancaster aggregation)", x="log2 Fold Change", y="-log10 adj p-value") +
  scale_color_manual(values=c("gray", "red")) +
  ggsave(file.path(output_dir, "volcano_plot_gene.png"), width=6, height=5)

# Heatmap of top 30 variable transcripts (optional)
top_trans <- head(order(res_df$padj), 30)
if (length(top_trans) >= 2) {
  vsd <- varianceStabilizingTransformation(dds)
  mat <- assay(vsd)[top_trans, ]
  pheatmap(mat, scale="row", annotation_col=metadata, filename=file.path(output_dir, "top_transcripts_heatmap.png"))
}

# -------- 7. SESSION INFO --------
sink(file.path(output_dir, "sessionInfo.txt"))
sessionInfo()
sink()

cat("RNA-seq transcript-level and gene-level (Lancaster) analysis complete. Results in", output_dir, "\n")
