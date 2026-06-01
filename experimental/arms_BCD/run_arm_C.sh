#!/bin/bash
# Arm C: Aggressive STAR parameters (maximum coverage)
set -euo pipefail
source ~/miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

OUTDIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/pilot/arm_C_star_aggressive"
FASTQ_DIR="/media/yusanghyeon/30B4E366B4E32D52/01_Projects/korea_PHMG/RNA_fastq"
STAR_INDEX="/home/yusanghyeon/RAT_project/PHMG_IT/reference/star_index"

echo "=== Arm C: Aggressive STAR alignment ==="

STAR \
  --genomeDir "$STAR_INDEX" \
  --readFilesIn "$FASTQ_DIR/C1_1.fastq.gz" "$FASTQ_DIR/C1_2.fastq.gz" \
  --readFilesCommand zcat \
  --runThreadN 14 \
  --twopassMode Basic \
  --alignEndsType Local \
  --alignSoftClipAtReferenceEnds Yes \
  --outFilterMultimapNmax 50 \
  --outMultimapperOrder Random \
  --outFilterMismatchNoverLmax 0.1 \
  --outFilterMismatchNmax 10 \
  --outFilterScoreMinOverLread 0.33 \
  --outFilterMatchNminOverLread 0.33 \
  --winAnchorMultimapNmax 100 \
  --alignIntronMin 20 \
  --alignIntronMax 1000000 \
  --alignMatesGapMax 1000000 \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMstrandField intronMotif \
  --outSAMattributes NH HI AS nM XS MD \
  --outSAMunmapped Within \
  --outSAMattrRGline ID:C1 SM:C1 PL:ILLUMINA LB:lib1 \
  --outFileNamePrefix "$OUTDIR/C1_" \
  --limitBAMsortRAM 50000000000

samtools index "$OUTDIR/C1_Aligned.sortedByCoord.out.bam"

echo "=== Arm C: Done ==="
