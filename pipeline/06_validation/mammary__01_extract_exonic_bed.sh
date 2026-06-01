#!/bin/bash
set -euo pipefail

# Step 1: rn6 GTF에서 exonic interval BED 추출
# Output: results/mutect2/intervals/exonic_rn6.bed
# 목적: MuTect2 calling을 exome 영역으로 제한 (WES validation에 맞춤)

WORK="/home/yusanghyeon/RAT_project/mammary_cancer"
GTF="${WORK}/reference/rn6/Rnor_6.0.104.gtf"
OUT_BED="${WORK}/results/mutect2/intervals/exonic_rn6.bed"
OUT_INTERVALS="${WORK}/results/mutect2/intervals/exonic_rn6.interval_list"

source /home/yusanghyeon/RAT_project/miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

echo "[$(date +%H:%M:%S)] Extracting exonic intervals from GTF..."

# GTF에서 exon feature만 추출 → BED (chr, start-1, end) → merge
awk -F'\t' 'BEGIN{OFS="\t"} $3=="exon" && $1 ~ /^([0-9]+|X|Y|MT)$/ {print $1, $4-1, $5}' "$GTF" \
  | sort -k1,1V -k2,2n -k3,3n \
  | bedtools merge -i - \
  > "$OUT_BED"

echo "Exonic BED:"
wc -l "$OUT_BED"
echo "Total exonic bp:"
awk '{sum += $3-$2} END {print sum}' "$OUT_BED"

# GATK용 interval_list 변환 (BED → Picard interval_list with header)
gatk BedToIntervalList \
  -I "$OUT_BED" \
  -O "$OUT_INTERVALS" \
  -SD "${WORK}/reference/rn6/Rnor_6.0.dict" 2>/dev/null

echo "[$(date +%H:%M:%S)] Done. Output:"
echo "  BED: $OUT_BED"
echo "  interval_list: $OUT_INTERVALS"
