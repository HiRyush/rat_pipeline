#!/bin/bash
# evaluate_coverage.sh — BAM의 coverage를 평가하고 DNA SNP 커버 수를 계산
#
# Usage: bash evaluate_coverage.sh <BAM> <DNA_SNP_VCF> <OUTPUT_DIR> [LABEL]

set -euo pipefail

BAM="$1"
DNA_VCF="$2"
OUT_DIR="$3"
LABEL="${4:-$(basename "$BAM" .bam)}"
REF="/home/yusanghyeon/RAT_project/PHMG_IT/reference/rn7.fa"
GENOME_SIZE=2648062885  # rn7 genome size in bp

mkdir -p "$OUT_DIR"

echo "=== Coverage Evaluation: $LABEL ==="
echo "BAM: $BAM"
echo "DNA VCF: $DNA_VCF"
echo ""

# 1. Flagstat
echo "[1/4] samtools flagstat..."
samtools flagstat "$BAM" > "$OUT_DIR/${LABEL}_flagstat.txt"
TOTAL_READS=$(grep "in total" "$OUT_DIR/${LABEL}_flagstat.txt" | awk '{print $1}')
MAPPED_READS=$(grep "mapped (" "$OUT_DIR/${LABEL}_flagstat.txt" | head -1 | awk '{print $1}')
MAP_RATE=$(grep "mapped (" "$OUT_DIR/${LABEL}_flagstat.txt" | head -1 | grep -oP '\([\d.]+' | tr -d '(')
echo "  Total reads: $TOTAL_READS"
echo "  Mapped reads: $MAPPED_READS ($MAP_RATE%)"

# 2. Coverage BED at multiple DP thresholds
echo "[2/4] bedtools genomecov..."
GENOMECOV_BG="$OUT_DIR/${LABEL}_genomecov.bg"
bedtools genomecov -ibam "$BAM" -bg > "$GENOMECOV_BG"

echo ""
echo "DP_threshold | Covered_Mb | Genome_pct"
echo "-------------|------------|----------"

for DP in 1 3 5 10; do
    BED="$OUT_DIR/${LABEL}_dp${DP}.bed"
    awk -v dp="$DP" '$4 >= dp' "$GENOMECOV_BG" | bedtools merge -i - > "$BED"
    COVERED_BP=$(awk '{sum += $3 - $2} END {print sum}' "$BED")
    COVERED_MB=$(echo "scale=1; $COVERED_BP / 1000000" | bc)
    GENOME_PCT=$(echo "scale=1; $COVERED_BP * 100 / $GENOME_SIZE" | bc)
    echo "DP>=${DP}        | ${COVERED_MB} Mb   | ${GENOME_PCT}%"
done

# 3. DNA SNP count in covered regions
echo ""
echo "[3/4] DNA SNPs in covered regions..."
echo ""
echo "DP_threshold | DNA_SNPs_in_region | Total_DNA_SNPs | Pct"
echo "-------------|--------------------|-----------|---------"

TOTAL_DNA=$(bcftools view -H "$DNA_VCF" | wc -l)

for DP in 1 3 5 10; do
    BED="$OUT_DIR/${LABEL}_dp${DP}.bed"
    SNP_COUNT=$(bcftools view -R "$BED" -H "$DNA_VCF" | wc -l)
    SNP_PCT=$(echo "scale=1; $SNP_COUNT * 100 / $TOTAL_DNA" | bc)
    echo "DP>=${DP}        | ${SNP_COUNT}          | ${TOTAL_DNA}  | ${SNP_PCT}%"
done

# 4. Summary line for easy parsing
echo ""
echo "[4/4] Summary JSON..."
# DP>=5 metrics for the summary
BED5="$OUT_DIR/${LABEL}_dp5.bed"
COV5_BP=$(awk '{sum += $3 - $2} END {print sum}' "$BED5")
COV5_MB=$(echo "scale=1; $COV5_BP / 1000000" | bc)
COV5_PCT=$(echo "scale=1; $COV5_BP * 100 / $GENOME_SIZE" | bc)
SNP5=$(bcftools view -R "$BED5" -H "$DNA_VCF" | wc -l)
SNP5_PCT=$(echo "scale=1; $SNP5 * 100 / $TOTAL_DNA" | bc)

# DP>=1 metrics
BED1="$OUT_DIR/${LABEL}_dp1.bed"
COV1_BP=$(awk '{sum += $3 - $2} END {print sum}' "$BED1")
COV1_MB=$(echo "scale=1; $COV1_BP / 1000000" | bc)
SNP1=$(bcftools view -R "$BED1" -H "$DNA_VCF" | wc -l)

echo "${LABEL}|${MAPPED_READS}|${MAP_RATE}|${COV1_MB}|${COV5_MB}|${COV5_PCT}|${SNP1}|${SNP5}|${SNP5_PCT}|${TOTAL_DNA}" > "$OUT_DIR/${LABEL}_summary.tsv"

echo "  Label: $LABEL"
echo "  Mapped: $MAPPED_READS ($MAP_RATE%)"
echo "  Coverage DP>=1: ${COV1_MB} Mb"
echo "  Coverage DP>=5: ${COV5_MB} Mb (${COV5_PCT}%)"
echo "  DNA SNPs DP>=1: ${SNP1}"
echo "  DNA SNPs DP>=5: ${SNP5} (${SNP5_PCT}%)"
echo ""
echo "=== Done: $LABEL ==="

# Cleanup large intermediate
rm -f "$GENOMECOV_BG"
