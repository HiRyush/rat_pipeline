#!/bin/bash
# Arm C alignment + coverage evaluation for remaining 9 samples (C4, C5, P4-P10)
# Step 1: Generate ground truth DNA SNP VCFs
# Step 2: STAR aggressive alignment + coverage evaluation
set -euo pipefail
source ~/miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

FASTQ_DIR="/media/yusanghyeon/30B4E366B4E32D52/01_Projects/korea_PHMG/RNA_fastq"
DNA_VCF_DIR="/media/yusanghyeon/30B4E366B4E32D52/01_Projects/korea_PHMG/DNA_vcf"
STAR_INDEX="/home/yusanghyeon/RAT_project/PHMG_IT/reference/star_index"
GT_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/ground_truth"
RESULT_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/pilot"
EVAL_SCRIPT="/home/yusanghyeon/RAT_project/PHMG_IT/scripts/evaluate_coverage.sh"

# Sample mapping: RNA_SAMPLE -> DNA_SAMPLE
declare -A DNA_MAP
DNA_MAP[C4]="P5S54W-6"
DNA_MAP[C5]="P5S54W-8"
DNA_MAP[P4]="P5H48W-12"
DNA_MAP[P5]="P5H49W-9"
DNA_MAP[P6]="P5H49W-18"
DNA_MAP[P7]="P5H52W-8"
DNA_MAP[P8]="P5H54W-2"
DNA_MAP[P9]="P5H54W-19"
DNA_MAP[P10]="P5M44W-11"

SAMPLES=("C4" "C5" "P4" "P5" "P6" "P7" "P8" "P9" "P10")
TOTAL=${#SAMPLES[@]}

echo "================================================================"
echo "Step 1: Generate ground truth DNA SNP VCFs"
echo "================================================================"

for S in "${SAMPLES[@]}"; do
    DNA_SAMPLE="${DNA_MAP[$S]}"
    DNA_VCF="${DNA_VCF_DIR}/${DNA_SAMPLE}_sorted.genome.vcf.gz"
    GT_OUT="${GT_DIR}/${S}_dna_snps.vcf.gz"

    if [ -f "$GT_OUT" ]; then
        echo "  $S: ground truth already exists, skipping"
        continue
    fi

    echo "  $S <- ${DNA_SAMPLE}: extracting SNPs..."
    bcftools view -v snps \
        -i 'QUAL>0 && (FILTER="PASS" || FILTER=".")' \
        "$DNA_VCF" \
        -Oz -o "$GT_OUT"
    bcftools index "$GT_OUT"
    echo "  $S: done"
done

echo ""
echo "================================================================"
echo "Step 2: STAR aggressive alignment + coverage evaluation"
echo "================================================================"

for i in "${!SAMPLES[@]}"; do
    S="${SAMPLES[$i]}"
    N=$((i+1))
    OUTDIR="${RESULT_DIR}/arm_C_${S}"
    mkdir -p "$OUTDIR"

    # Skip if BAM already exists
    BAM="${OUTDIR}/${S}_Aligned.sortedByCoord.out.bam"
    if [ -f "$BAM" ] && [ -f "${BAM}.bai" ]; then
        echo "[$N/$TOTAL] $S: BAM already exists, skipping alignment"
    else
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
        samtools index "$BAM"
    fi

    # Evaluate coverage
    GT_VCF="${GT_DIR}/${S}_dna_snps.vcf.gz"
    SUMMARY="${OUTDIR}/arm_C_${S}_summary.tsv"
    if [ -f "$SUMMARY" ]; then
        echo "[$N/$TOTAL] $S: coverage already evaluated, skipping"
    else
        echo "[$N/$TOTAL] Evaluating coverage..."
        bash "$EVAL_SCRIPT" "$BAM" "$GT_VCF" "$OUTDIR" "arm_C_${S}"
    fi

    # Cleanup STAR temp directories
    rm -rf "${OUTDIR}/${S}__STARgenome" "${OUTDIR}/${S}__STARpass1" 2>/dev/null

    echo "[$N/$TOTAL] $S complete!"
done

echo ""
echo "================================================================"
echo "ALL DONE! Generating full 15-sample summary..."
echo "================================================================"

# Collect all 15 samples summary
SUMMARY_FILE="${RESULT_DIR}/batch_summary_all15.tsv"
echo "Sample|Mapped|MapRate|Cov_DP1_Mb|Cov_DP5_Mb|Cov_DP5_Pct|SNP_DP1|SNP_DP5|SNP_DP5_Pct|Total_DNA" > "$SUMMARY_FILE"

# C1 (original arm_C_star_aggressive)
cat "${RESULT_DIR}/arm_C_star_aggressive/arm_C_summary.tsv" >> "$SUMMARY_FILE"

# C2-C5, P1-P10
for S in C2 C3 C4 C5 P1 P2 P3 P4 P5 P6 P7 P8 P9 P10; do
    cat "${RESULT_DIR}/arm_C_${S}/arm_C_${S}_summary.tsv" >> "$SUMMARY_FILE"
done

echo "Summary saved to: $SUMMARY_FILE"
echo ""
cat "$SUMMARY_FILE" | column -t -s'|'
