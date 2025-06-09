#!/bin/bash

###############################################################################
# ATAC-seq Analysis Pipeline
# Author: Alireza Ghahramani
# Contact: aghahram@uwo.ca
#
# Description:
#   This script processes bulk ATAC-seq fastq files through adapter trimming,
#   alignment, deduplication, peak calling, coverage track generation, and
#   transcription factor (TF) footprinting/motif analysis using TOBIAS.
#
# Main Steps:
#   1. Adapter/quality trimming with Trim Galore (+FastQC)
#   2. Alignment with Bowtie2
#   3. Sorting, indexing, and deduplication of BAM files
#   4. Peak calling with MACS2 (using standard ATAC-seq settings)
#   5. Coverage track generation with deepTools (bamCoverage)
#   6. TF footprinting and motif analysis with TOBIAS (ATACorrect, FootprintScores, BindDetect)
#
# Requirements:
#   - Software: Trim Galore, FastQC, Bowtie2, SAMtools, MACS2, deepTools (bamCoverage), TOBIAS
#   - Bowtie2 genome index and reference genome FASTA (for TOBIAS)
#   - Consensus/merged peak BED file (for TOBIAS)
#   - Motif database (e.g., JASPAR PFM file for TOBIAS BindDetect)
#
# Usage:
#   1. Set all variables under "USER PARAMETERS" below.
#   2. Place raw paired-end fastq files as *_R1.fastq.gz and *_R2.fastq.gz in $RAW_DIR.
#   3. Run: bash atac_seq_pipeline.sh
#
# Output:
#   - Output directories: trimmed fastq, aligned BAMs, peaks, bigwig tracks, TOBIAS results, logs
#
###############################################################################

# ---------------- USER PARAMETERS (EDIT BELOW) ----------------

BOWTIE2_INDEX="/path/to/bowtie2_index/genome"        # Path to Bowtie2 genome index prefix
GENOME_FASTA="/path/to/genome.fa"                    # Path to reference genome FASTA (for TOBIAS)
ANNOTATION_GTF="/path/to/annotation.gtf"             # Optional, if using for downstream annotation
THREADS=8                                            # Number of threads for parallel steps

RAW_DIR="raw_fastq"                                  # Folder with input FASTQ files
OUT_DIR="atac_output"                                # Main output directory

TRIM_DIR="${OUT_DIR}/trimmed"
ALIGN_DIR="${OUT_DIR}/aligned"
PEAK_DIR="${OUT_DIR}/peaks"
BIGWIG_DIR="${OUT_DIR}/bigwig"
TOBIAS_DIR="${OUT_DIR}/tobias"
LOG_DIR="${OUT_DIR}/logs"

CONSENSUS_PEAKS="/path/to/consensus_peaks.bed"       # Consensus/merged peak BED (for TOBIAS, e.g., union of all samples)
MOTIF_PFM="/path/to/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt" # JASPAR motifs for TOBIAS

# Create all output directories if they do not exist
mkdir -p $TRIM_DIR $ALIGN_DIR $PEAK_DIR $BIGWIG_DIR $TOBIAS_DIR $LOG_DIR

# ------------------ MAIN PIPELINE ------------------

for fq1 in ${RAW_DIR}/*_R1.fastq.gz; do
    fq2="${fq1/_R1.fastq.gz/_R2.fastq.gz}"
    sample=$(basename $fq1 _R1.fastq.gz)
    echo "---- Processing sample: $sample ----" | tee -a $LOG_DIR/pipeline.log

    # 1. Adapter/quality trimming with Trim Galore (includes FastQC)
    trim_galore --paired --fastqc --output_dir $TRIM_DIR $fq1 $fq2 > $LOG_DIR/${sample}_trim_galore.log 2>&1
    trimmed1="${TRIM_DIR}/${sample}_R1_val_1.fq.gz"
    trimmed2="${TRIM_DIR}/${sample}_R2_val_2.fq.gz"

    # 2. Alignment with Bowtie2 (paired-end)
    bowtie2 -p $THREADS -x $BOWTIE2_INDEX -1 $trimmed1 -2 $trimmed2 | \
        samtools view -bS - | \
        samtools sort -@ $THREADS -o $ALIGN_DIR/${sample}.sorted.bam
    samtools index $ALIGN_DIR/${sample}.sorted.bam

    # 3. Deduplication (removal of PCR duplicates, important for ATAC-seq)
    samtools rmdup $ALIGN_DIR/${sample}.sorted.bam $ALIGN_DIR/${sample}.rmdup.bam
    samtools index $ALIGN_DIR/${sample}.rmdup.bam

    # 4. Peak calling with MACS2 (standard ATAC-seq settings)
    macs2 callpeak -t $ALIGN_DIR/${sample}.rmdup.bam \
        -f BAMPE -g mm -n $sample \
        --outdir $PEAK_DIR --nomodel --shift -100 --extsize 200 -q 0.01 > $LOG_DIR/${sample}_macs2.log 2>&1

    # 5. Generate normalized BigWig coverage tracks (deepTools)
    bamCoverage -b $ALIGN_DIR/${sample}.rmdup.bam \
        -o $BIGWIG_DIR/${sample}.bw \
        --normalizeUsing RPGC \
        --effectiveGenomeSize 2652783500 \
        --binSize 10 \
        --ignoreDuplicates

    echo "Sample $sample complete." | tee -a $LOG_DIR/pipeline.log
done

# ------------------ TOBIAS ANALYSIS (TF FOOTPRINTING) ------------------

# The following steps assume:
# - You have the reference genome FASTA, a consensus/merged peak BED, and sample BAM files.
# - All commands are run after activating the TOBIAS environment.

for bam in $ALIGN_DIR/*.rmdup.bam; do
    sample=$(basename $bam .rmdup.bam)
    echo "---- TOBIAS Analysis for: $sample ----" | tee -a $LOG_DIR/pipeline.log

    # 1. Correction for Tn5 insertion bias (ATACorrect)
    TOBIAS ATACorrect \
        --bam $bam \
        --genome $GENOME_FASTA \
        --peaks $CONSENSUS_PEAKS \
        --outdir $TOBIAS_DIR/${sample}_atacorrect/ \
        --cores $THREADS \
        > $LOG_DIR/${sample}_tobias_atacorrect.log 2>&1

    # 2. Calculate footprint scores (FootprintScores)
    TOBIAS FootprintScores \
        --signal $TOBIAS_DIR/${sample}_atacorrect/corrected.bw \
        --regions $CONSENSUS_PEAKS \
        --output $TOBIAS_DIR/${sample}_footprints.bw \
        --cores $THREADS \
        > $LOG_DIR/${sample}_tobias_footprintscores.log 2>&1

    # 3. Detect TF binding sites using motifs (BindDetect)
    TOBIAS BindDetect \
        --motifs $MOTIF_PFM \
        --signals $TOBIAS_DIR/${sample}_footprints.bw \
        --regions $CONSENSUS_PEAKS \
        --genome $GENOME_FASTA \
        --outdir $TOBIAS_DIR/${sample}_bindetect/ \
        > $LOG_DIR/${sample}_tobias_bindetect.log 2>&1

    echo "TOBIAS TF footprinting complete for $sample." | tee -a $LOG_DIR/pipeline.log
done

echo "ATAC-seq pipeline (including TOBIAS) completed for all samples."
