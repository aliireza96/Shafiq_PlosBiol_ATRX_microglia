# Shafiq_PlosBiol_ATRX_microglia
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15679297.svg)](https://doi.org/10.5281/zenodo.15679297)


# ATRX-null Microglia Multi-omics Analysis

This repository contains the scripts and pipelines used in our study of ATRX deletion in microglia. The analysis integrates RNA-seq, ATAC-seq, and transposable element (TE) profiling, with dedicated workflows for LTR classification and data visualization.

## Organization

- **/RNAseq/**: RNA-seq preprocessing, gene and transcript differential expression, TE quantification by two independent routes, and LTR solo/full-length classification.
- **/ATACseq/**: ATAC-seq preprocessing, peak calling, and differential analysis of consensus peaks and of TEs.

## TE analysis: two independent routes

TE results in the manuscript rest on agreement between two quantification strategies, not on either alone. They answer different questions and fail in different ways.

| | Per-locus | Per-subfamily |
|---|---|---|
| Script | `RNAseq/te_featurecounts_differential_expression.R` | `RNAseq/rnaseq_TEtranscripts_differential_analysis.R` |
| Tool | featureCounts (Rsubread), UCSC RepeatMasker SAF | TEtranscripts, `--mode multi` |
| Resolution | individual TE loci | TE subfamily |
| Multi-mapping reads | counted, and separately excluded, as a sensitivity check | redistributed probabilistically |
| Strand | BAMs split into sense/antisense first (`RNAseq/split_bam_and_bigwig.sh`), each tested separately, significant loci intersected | `--stranded reverse` |

The per-locus route shows which elements carry the signal. The subfamily route gains power on repetitive sequence where unique placement is impossible. The reported LTR upregulation is present in both.

A reviewer asked that proviruses be distinguished from solo LTRs, since the two support different claims about ERV-driven viral mimicry. That classification is `RNAseq/ltr_classification_solo_full_length.R`.

## Library chemistry

The RNA-seq libraries are **reverse-stranded (dUTP)**. This determines several parameters, so check it against your own data before reusing these scripts: `--stranded reverse` for TEtranscripts, the flag combinations in `split_bam_and_bigwig.sh` (sense = 99/147, antisense = 83/163), and `-s 2` wherever featureCounts is run on unsplit BAMs.

## Reusing these scripts on your own data

The R scripts take a counts file and a sample metadata CSV, and expose their settings in a `USER PARAMETERS` block at the top. Two things to set deliberately:

- `reference_level` and `treatment_level`. These define the contrast. If you leave DESeq2 to order the levels itself it does so alphabetically, which silently decides the sign of every log2 fold change.
- Your aligner's multi-mapping settings. TEtranscripts `--mode multi` can only redistribute alignments that the aligner actually reported. Aligners that cap reported loci at a low default discard exactly the young, high-copy elements (L1Md, IAPEz, MMERVK10C) that TE studies care about, and the loss is not uniform across subfamilies, so it cannot be corrected downstream.

## Changes since deposition

The scripts were updated after publication. The analysis and the published results are unchanged.

- TEtranscripts runs with `--stranded reverse`, matching the dUTP reverse-stranded libraries.
- The TEtranscripts input is `<prefix>.cntTable`. That table holds genes and TEs together, so both are kept for size-factor estimation and TEs are separated at the results stage.
- The TEtranscripts header records the alignment used in the manuscript Methods: Trim Galore v0.6.6 and HISAT2 v2.0.4 against GRCm38.p6, per `RNAseq/rna_seq_preprocessing_pipeline.sh`.
- `--minl 50` dropped, `--sortByPos` added.
- All five DESeq2 scripts take an explicit contrast through `reference_level` and `treatment_level`, and relevel before fitting.
- `ggsave()` is called on its own rather than chained onto the plot object.
- `pheatmap` receives the grouping column alone as `annotation_col`.
- Volcano plots filter out features whose `padj` is NA from independent filtering.
- The per-locus script documents the strand-splitting step from the Methods and points at `split_bam_and_bigwig.sh`.
- `rnaseq_TEtranscripts_differential_analysis.R` no longer has a leading space in its filename.

The Zenodo archive (DOI above) is a snapshot taken at deposition and does not carry these changes. Cite this repository for the current scripts.

**Questions or issues?**
Open an Issue on GitHub or email Alireza Ghahramani at aghahram@uwo.ca.
