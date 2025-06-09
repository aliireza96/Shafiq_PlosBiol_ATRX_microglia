#!/usr/bin/env bash
###############################################################################
# Script: split_bam_and_bigwig.sh
#
# Purpose:
#   For each paired-end, stranded RNA-seq BAM:
#     1. Split by *transcriptional sense and antisense* using proper flag combinations
#        (dUTP/reverse-stranded protocol standard, see below)
#     2. Generate indexed BAMs for sense (fwd) and antisense (rev)
#     3. Create normalized BigWig (BW) tracks for each strand using deepTools
#
# Usage:
#   bash split_bam_and_bigwig.sh <input.bam> <output_directory>
#
# Requirements:
#   - samtools (v1.9+)
#   - deepTools (bamCoverage)
#
# Detailed BAM splitting logic (stranded, paired-end, dUTP "reverse" protocol):
#   FORWARD (sense) BAM:
#     - Read 1, properly paired, mapped to reverse strand:     flag 147
#     - Read 2, properly paired, mapped to forward strand:     flag 99
#   REVERSE (antisense) BAM:
#     - Read 1, properly paired, mapped to forward strand:     flag 83
#     - Read 2, properly paired, mapped to reverse strand:     flag 163
#
# Output:
#   - For each <input.bam>, creates:
#       <name>.fwd.bam, <name>.rev.bam, and their .bai indexes
#       <name>.fwd.bw, <name>.rev.bw   (BigWig coverage tracks)
###############################################################################

set -euo pipefail

BAM=$1
TARGET_D=$2

FILE=$(basename "$BAM")
NAME=${FILE%.*}
BAMF1=${TARGET_D}/${NAME}.fwd1.bam
BAMF2=${TARGET_D}/${NAME}.fwd2.bam
BAMF=${TARGET_D}/${NAME}.fwd.bam
BAMR1=${TARGET_D}/${NAME}.rev1.bam
BAMR2=${TARGET_D}/${NAME}.rev2.bam
BAMR=${TARGET_D}/${NAME}.rev.bam

# Create output directory if not exist
mkdir -p "$TARGET_D"

echo "Splitting $BAM by transcriptional strand..."

# FORWARD/SENSE: (flag 147 and 99)
samtools view -@ 8 -bh -f 147 "$BAM" > "$BAMF1"
samtools index -@ 8 "$BAMF1"
samtools view -@ 8 -bh -f 99 "$BAM" > "$BAMF2"
samtools index -@ 8 "$BAMF2"
samtools merge -@ 8 -f "$BAMF" "$BAMF1" "$BAMF2"
samtools index -@ 8 "$BAMF"

# REVERSE/ANTISENSE: (flag 83 and 163)
samtools view -@ 8 -bh -f 83 "$BAM" > "$BAMR1"
samtools index -@ 8 "$BAMR1"
samtools view -@ 8 -bh -f 163 "$BAM" > "$BAMR2"
samtools index -@ 8 "$BAMR2"
samtools merge -@ 8 -f "$BAMR" "$BAMR1" "$BAMR2"
samtools index -@ 8 "$BAMR"

# Optionally remove temporary split BAMs (uncomment if you want cleanup)
# rm "$BAMF1" "$BAMF2" "$BAMR1" "$BAMR2"

echo "Generating BigWig tracks for UCSC/IGV..."

# Create normalized (CPM) BigWig for sense and antisense
bamCoverage -b "$BAMF" -o "${TARGET_D}/${NAME}.fwd.bw" --binSize 10 --normalizeUsing CPM
bamCoverage -b "$BAMR" -o "${TARGET_D}/${NAME}.rev.bw" --binSize 10 --normalizeUsing CPM

echo "Done. Outputs in $TARGET_D:"
echo "  - ${NAME}.fwd.bam, ${NAME}.fwd.bw (sense/forward)"
echo "  - ${NAME}.rev.bam, ${NAME}.rev.bw (antisense/reverse)"
