#!/bin/bash

# RNA-seq processing pipeline (multi-sample, robust, reproducible)

# -------- USER CONFIG --------
HISAT2_INDEX="/path/to/hisat2_index/genome"
ANNOTATION_GTF="/path/to/annotation.gtf"
THREADS=8
RAW_DIR="raw_fastq"
OUT_DIR="analysis_output"
TRIM_DIR="${OUT_DIR}/trimmed"
ALIGN_DIR="${OUT_DIR}/aligned"
STRINGTIE_DIR="${OUT_DIR}/stringtie"
COUNT_DIR="${OUT_DIR}/featureCounts"
FASTQC_DIR="${OUT_DIR}/fastqc"
LOG_DIR="${OUT_DIR}/logs"

# Create directories
mkdir -p $TRIM_DIR $ALIGN_DIR $STRINGTIE_DIR $COUNT_DIR $FASTQC_DIR $LOG_DIR

# -------- MAIN PIPELINE --------

for fq1 in ${RAW_DIR}/*_R1.fastq.gz; do
    fq2="${fq1/_R1.fastq.gz/_R2.fastq.gz}"
    sample=$(basename $fq1 _R1.fastq.gz)
    echo "Processing sample: $sample" | tee -a $LOG_DIR/pipeline.log

    # 1. Trim Galore & FastQC
    trim_galore --paired --fastqc --output_dir $TRIM_DIR $fq1 $fq2 > $LOG_DIR/${sample}_trim_galore.log 2>&1
    trimmed1="${TRIM_DIR}/${sample}_R1_val_1.fq.gz"
    trimmed2="${TRIM_DIR}/${sample}_R2_val_2.fq.gz"

    # Move FastQC output to separate directory
    mv $TRIM_DIR/*_fastqc* $FASTQC_DIR/ 2>/dev/null

    # 2. HISAT2 alignment
    hisat2 -p $THREADS -x $HISAT2_INDEX -1 $trimmed1 -2 $trimmed2 -S $ALIGN_DIR/${sample}.sam > $LOG_DIR/${sample}_hisat2.log 2>&1

    # 3. SAMtools: SAM -> BAM, sort, index
    samtools view -bS $ALIGN_DIR/${sample}.sam | samtools sort -o $ALIGN_DIR/${sample}.sorted.bam
    samtools index $ALIGN_DIR/${sample}.sorted.bam

    # Optional: remove large SAM to save space
    rm $ALIGN_DIR/${sample}.sam

    # 4. StringTie for transcript quantification
    stringtie $ALIGN_DIR/${sample}.sorted.bam -p $THREADS -G $ANNOTATION_GTF -o $STRINGTIE_DIR/${sample}.gtf -A $STRINGTIE_DIR/${sample}.abund.tab > $LOG_DIR/${sample}_stringtie.log 2>&1

    # 5. FeatureCounts for gene-level quantification
    featureCounts -p -T $THREADS -a $ANNOTATION_GTF -o $COUNT_DIR/${sample}.counts.txt $ALIGN_DIR/${sample}.sorted.bam > $LOG_DIR/${sample}_featureCounts.log 2>&1

    echo "Sample $sample complete." | tee -a $LOG_DIR/pipeline.log
done

echo "All samples processed."
