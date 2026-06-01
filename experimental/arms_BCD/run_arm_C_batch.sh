#!/bin/bash
# Arm C (aggressive STAR) alignment + coverage evaluation for multiple samples
set -euo pipefail
source ~/miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

FASTQ_DIR="/media/yusanghyeon/30B4E366B4E32D52/01_Projects/korea_PHMG/RNA_fastq"
STAR_INDEX="/home/yusanghyeon/RAT_project/PHMG_IT/reference/star_index"
GT_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/ground_truth"
RESULT_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/pilot"
EVAL_SCRIPT="/home/yusanghyeon/RAT_project/PHMG_IT/scripts/evaluate_coverage.sh"

SAMPLES=("C2" "C3" "P1" "P2" "P3")
TOTAL=${#SAMPLES[@]}

for i in "${!SAMPLES[@]}"; do
    S="${SAMPLES[$i]}"
    N=$((i+1))
    OUTDIR="${RESULT_DIR}/arm_C_${S}"
    mkdir -p "$OUTDIR"

    echo ""
    echo "================================================================"
    echo "[$N/$TOTAL] Sample: $S — STAR aggressive alignment"
    echo "================================================================"

    STAR \
      --genomeDir "$STAR_INDEX" \
      --readFilesIn "$FASTQ_DIR/${S}_1.fastq.gz" "$FASTQ_DIR/${S}_2.fastq.gz" \
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
      --outSAMattrRGline ID:${S} SM:${S} PL:ILLUMINA LB:lib1 \
      --outFileNamePrefix "${OUTDIR}/${S}_" \
      --limitBAMsortRAM 50000000000

    echo "[$N/$TOTAL] Indexing BAM..."
    samtools index "${OUTDIR}/${S}_Aligned.sortedByCoord.out.bam"

    echo "[$N/$TOTAL] Evaluating coverage..."
    bash "$EVAL_SCRIPT" \
      "${OUTDIR}/${S}_Aligned.sortedByCoord.out.bam" \
      "${GT_DIR}/${S}_dna_snps.vcf.gz" \
      "$OUTDIR" \
      "arm_C_${S}"

    # Cleanup STAR temp directories
    rm -rf "${OUTDIR}/${S}__STARgenome" "${OUTDIR}/${S}__STARpass1" 2>/dev/null

    echo "[$N/$TOTAL] $S complete!"
done

echo ""
echo "================================================================"
echo "ALL DONE! Generating summary..."
echo "================================================================"

# Collect summaries
echo "Sample|Mapped|MapRate|Cov_DP1_Mb|Cov_DP5_Mb|Cov_DP5_Pct|SNP_DP1|SNP_DP5|SNP_DP5_Pct|Total_DNA" > "${RESULT_DIR}/batch_summary.tsv"
# Add C1 (already done)
cat "${RESULT_DIR}/arm_C_star_aggressive/arm_C_summary.tsv" >> "${RESULT_DIR}/batch_summary.tsv"
for S in "${SAMPLES[@]}"; do
    cat "${RESULT_DIR}/arm_C_${S}/arm_C_${S}_summary.tsv" >> "${RESULT_DIR}/batch_summary.tsv"
done
echo "Summary saved to: ${RESULT_DIR}/batch_summary.tsv"
cat "${RESULT_DIR}/batch_summary.tsv"
