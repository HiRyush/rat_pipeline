#!/bin/bash
# Multi-sample callable region analysis
# "15개 중 K개 샘플에서 DP≥1"이면 callable로 정의
# 각 K 기준별 coverage와 DNA SNP coverage 측정
set -euo pipefail
source ~/miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

PILOT_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/pilot"
GT_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/ground_truth"
OUT_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/multi_sample"
GENOME_SIZE=2648062885
REF="/home/yusanghyeon/RAT_project/PHMG_IT/reference/rn7.fa"

mkdir -p "$OUT_DIR"

echo "================================================================"
echo "Multi-sample Callable Region Analysis"
echo "================================================================"

# Step 1: Create bedtools genome file from .fai
GENOME_FILE="$OUT_DIR/rn7.genome"
awk '{print $1"\t"$2}' "${REF}.fai" > "$GENOME_FILE"

# Step 2: Convert each sample's DP≥1 BED to per-base coverage (0/1) using genomecov
# Then sum across samples using bedtools unionbedg
echo "[1/3] Preparing per-sample BED files..."

BED_LIST=()
SAMPLE_NAMES=()
for S in C1 C2 C3 C4 C5 P1 P2 P3 P4 P5 P6 P7 P8 P9 P10; do
    if [ "$S" = "C1" ]; then
        BED="${PILOT_DIR}/arm_C_star_aggressive/arm_C_dp1.bed"
    else
        BED="${PILOT_DIR}/arm_C_${S}/arm_C_${S}_dp1.bed"
    fi

    # Convert to sorted, genome-clipped bedGraph (value=1 where covered)
    SORTED_BG="$OUT_DIR/${S}_dp1_sorted.bg"
    awk '{print $1"\t"$2"\t"$3"\t1"}' "$BED" | sort -k1,1 -k2,2n > "$SORTED_BG"
    BED_LIST+=("$SORTED_BG")
    SAMPLE_NAMES+=("$S")
    echo "  $S ready"
done

# Step 3: Use bedtools unionbedg to merge all 15 samples
echo ""
echo "[2/3] Computing unionbedg across 15 samples..."
UNION_BG="$OUT_DIR/union_15samples.bg"
bedtools unionbedg -i "${BED_LIST[@]}" -g "$GENOME_FILE" -filler 0 \
    -header -names "${SAMPLE_NAMES[@]}" > "$UNION_BG"

echo "  Union BedGraph created: $(wc -l < "$UNION_BG") regions"

# Step 4: For each K threshold, compute callable regions
echo ""
echo "[3/3] Computing callable regions at each K threshold..."
echo ""

SUMMARY_FILE="$OUT_DIR/multi_sample_callable_summary.tsv"
echo "K_threshold|Callable_Mb|Callable_Pct|Description" > "$SUMMARY_FILE"

for K in 1 2 3 5 8 10 12 15; do
    CALLABLE_BED="$OUT_DIR/callable_K${K}_of_15.bed"

    # Skip header line, sum columns 4-18 (15 samples), filter >= K
    awk -v k="$K" 'NR>1 {
        sum=0; for(i=4;i<=NF;i++) sum+=$i;
        if(sum>=k) print $1"\t"$2"\t"$3
    }' "$UNION_BG" | bedtools merge -i - > "$CALLABLE_BED"

    CALLABLE_BP=$(awk '{sum += $3 - $2} END {print sum+0}' "$CALLABLE_BED")
    CALLABLE_MB=$(echo "scale=1; $CALLABLE_BP / 1000000" | bc)
    CALLABLE_PCT=$(echo "scale=1; $CALLABLE_BP * 100 / $GENOME_SIZE" | bc)

    echo "K>=${K} (${K}/15 samples DP≥1): ${CALLABLE_MB} Mb (${CALLABLE_PCT}%)"
    echo "${K}|${CALLABLE_MB}|${CALLABLE_PCT}|${K}/15 samples DP>=1" >> "$SUMMARY_FILE"
done

# Step 5: DNA SNP coverage at each K threshold (using all 15 ground truth VCFs)
echo ""
echo "================================================================"
echo "DNA SNP coverage at each K threshold"
echo "================================================================"
echo ""

SNP_SUMMARY="$OUT_DIR/multi_sample_snp_coverage.tsv"
echo "K_threshold|Sample|SNPs_in_callable|Total_SNPs|Pct" > "$SNP_SUMMARY"

# Use 6 samples with ground truth for validation
GT_SAMPLES=("C1" "C2" "C3" "P1" "P2" "P3")
# Also include new ones
for S in C4 C5 P4 P5 P6 P7 P8 P9 P10; do
    if [ -f "${GT_DIR}/${S}_dna_snps.vcf.gz" ]; then
        GT_SAMPLES+=("$S")
    fi
done

echo "Samples with ground truth: ${GT_SAMPLES[*]}"
echo ""

for K in 1 3 5 8 10 12 15; do
    CALLABLE_BED="$OUT_DIR/callable_K${K}_of_15.bed"

    TOTAL_SNP_SUM=0
    COVERED_SNP_SUM=0

    for S in "${GT_SAMPLES[@]}"; do
        GT_VCF="${GT_DIR}/${S}_dna_snps.vcf.gz"
        TOTAL_SNP=$(bcftools view -H "$GT_VCF" | wc -l)
        COVERED_SNP=$(bcftools view -R "$CALLABLE_BED" -H "$GT_VCF" | wc -l)
        SNP_PCT=$(echo "scale=1; $COVERED_SNP * 100 / $TOTAL_SNP" | bc)

        echo "${K}|${S}|${COVERED_SNP}|${TOTAL_SNP}|${SNP_PCT}" >> "$SNP_SUMMARY"

        TOTAL_SNP_SUM=$((TOTAL_SNP_SUM + TOTAL_SNP))
        COVERED_SNP_SUM=$((COVERED_SNP_SUM + COVERED_SNP))
    done

    AVG_PCT=$(echo "scale=1; $COVERED_SNP_SUM * 100 / $TOTAL_SNP_SUM" | bc)
    echo "K>=${K}: DNA SNP coverage = ${AVG_PCT}% (avg across ${#GT_SAMPLES[@]} samples)"
done

# Cleanup intermediate files
echo ""
echo "Cleaning up intermediate files..."
rm -f "$OUT_DIR"/*_dp1_sorted.bg

echo ""
echo "================================================================"
echo "DONE! Results in: $OUT_DIR/"
echo "  - multi_sample_callable_summary.tsv"
echo "  - multi_sample_snp_coverage.tsv"
echo "  - callable_K{1,2,3,5,8,10,12,15}_of_15.bed"
echo "================================================================"
