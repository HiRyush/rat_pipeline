#!/bin/bash
set -euo pipefail

# Step 9: Sample identity verification via genotype concordance
# 각 RNA sample의 PASS 변이를 모든 WES와 pairwise 비교
# 매칭된 짝(같은 개체)이라면 변이 overlap이 다른 짝보다 압도적으로 높아야 함

WORK="/home/yusanghyeon/RAT_project/mammary_cancer"
FILT_DIR="${WORK}/results/mutect2/filtered"
GT_DIR="${WORK}/ground_truth"
VAL_DIR="${WORK}/results/observation_first/validation"
OUT_DIR="${WORK}/results/observation_first/identity_check"
mkdir -p "$OUT_DIR"

source /home/yusanghyeon/RAT_project/miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

# Sample → claimed WES 매핑 (from sample_mapping_verified.tsv)
declare -A CLAIMED
CLAIMED[SRR33625148]="GSM8994506_20m_4"
CLAIMED[SRR33625149]="GSM8994505_20m_3"
CLAIMED[SRR33625150]="GSM8994504_20m_2"
CLAIMED[SRR33625151]="GSM8994503_20m_1"
CLAIMED[SRR33625152]="GSM8994502_12m_3"
CLAIMED[SRR33625153]="GSM8994501_12m_2"
CLAIMED[SRR33625154]="GSM8994500_12m_1"
CLAIMED[SRR33625165]="GSM8994488_1m_7"
CLAIMED[SRR33625166]="GSM8994487_1m_6"
CLAIMED[SRR33625167]="GSM8994486_1m_5"
CLAIMED[SRR33625168]="GSM8994484_1m_3"
CLAIMED[SRR33625169]="GSM8994483_1m_2"
CLAIMED[SRR33625170]="GSM8994482_1m_1"

SAMPLES=(SRR33625148 SRR33625149 SRR33625150 SRR33625151 SRR33625152 SRR33625153 SRR33625154 \
         SRR33625165 SRR33625166 SRR33625167 SRR33625168 SRR33625169 SRR33625170)

# Step 1: RNA sample 별 position+allele set 추출 (CHROM\tPOS\tREF\tALT)
echo "[$(date +%H:%M:%S)] Step 1: RNA variant set 추출..."
for srr in "${SAMPLES[@]}"; do
  PASS="${FILT_DIR}/${srr}.pass.vcf.gz"
  bcftools view -H "$PASS" 2>/dev/null \
    | awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $4, $5}' \
    | sort -u > "${OUT_DIR}/rna_${srr}.set"
done

# Step 2: WES 별 position+allele set 추출 (gzip 처리)
echo "[$(date +%H:%M:%S)] Step 2: WES variant set 추출..."
for vcf in "$GT_DIR"/*_tumor.pass.vcf; do
  base=$(basename "$vcf" _tumor.pass.vcf)
  awk -F'\t' '!/^#/{print $1"\t"$2"\t"$4"\t"$5}' "$vcf" \
    | sort -u > "${OUT_DIR}/wes_${base}.set"
done

# Step 3: pairwise overlap 매트릭스
echo "[$(date +%H:%M:%S)] Step 3: Pairwise overlap 계산..."
MATRIX="${OUT_DIR}/identity_matrix.tsv"

# Header
{
  printf "RNA_sample"
  for vcf in "$GT_DIR"/*_tumor.pass.vcf; do
    base=$(basename "$vcf" _tumor.pass.vcf)
    printf "\t%s" "$base"
  done
  printf "\tCLAIMED_match\tBEST_match\tMATCH?\n"
} > "$MATRIX"

# 각 RNA에 대해 모든 WES와 overlap 계산
for srr in "${SAMPLES[@]}"; do
  RNA_SET="${OUT_DIR}/rna_${srr}.set"
  RNA_N=$(wc -l < "$RNA_SET")
  printf "%s" "$srr" >> "$MATRIX"

  best_match=""
  best_count=0

  for vcf in "$GT_DIR"/*_tumor.pass.vcf; do
    base=$(basename "$vcf" _tumor.pass.vcf)
    WES_SET="${OUT_DIR}/wes_${base}.set"
    # Overlap count
    overlap=$(comm -12 "$RNA_SET" "$WES_SET" | wc -l)
    printf "\t%d" "$overlap" >> "$MATRIX"

    if [ "$overlap" -gt "$best_count" ]; then
      best_count=$overlap
      best_match=$base
    fi
  done

  claimed="${CLAIMED[$srr]}"
  match_result="?"
  if [ "$best_match" = "$claimed" ]; then
    match_result="✅"
  else
    match_result="❌ (claimed=${claimed}, best=${best_match})"
  fi

  printf "\t%s\t%s\t%s\n" "$claimed" "$best_match" "$match_result" >> "$MATRIX"
done

echo ""
echo "[$(date +%H:%M:%S)] Step 4: Summary 출력..."
echo ""
echo "=== Sample Identity Check Results ==="
echo ""
echo "각 RNA sample의 claimed match vs 실제 best match:"
echo ""
awk -F'\t' 'NR==1{
  # find columns
  for(i=1; i<=NF; i++) {
    if ($i == "CLAIMED_match") claimed_col=i
    if ($i == "BEST_match") best_col=i
    if ($i == "MATCH?") match_col=i
  }
  next
}
{
  printf "  %-15s  claimed=%-20s  best=%-20s  %s\n", $1, $claimed_col, $best_col, $match_col
}' "$MATRIX"

# Top 3 matches per sample
echo ""
echo "=== 각 RNA의 Top 3 overlap (relative) ==="
HEADERS=$(head -1 "$MATRIX" | tr '\t' '\n')
NCOL=$(head -1 "$MATRIX" | awk -F'\t' '{print NF}')

awk -F'\t' -v ncol=$NCOL 'NR==1{
  for(i=1; i<=NF; i++) hdr[i]=$i
  next
}
{
  rna=$1
  # Last 3 cols are summary, so iterate 2..(NF-3)
  delete vals
  delete labels
  n=0
  for(i=2; i<=NF-3; i++) {
    n++
    vals[n]=$i
    labels[n]=hdr[i]
  }
  # Sort
  for(i=1; i<=n; i++) {
    for(j=i+1; j<=n; j++) {
      if (vals[j] > vals[i]) {
        t=vals[i]; vals[i]=vals[j]; vals[j]=t
        t=labels[i]; labels[i]=labels[j]; labels[j]=t
      }
    }
  }
  printf "  %s: ", rna
  for(i=1; i<=3 && i<=n; i++) {
    printf "%s=%d  ", labels[i], vals[i]
  }
  printf "\n"
}' "$MATRIX"

echo ""
echo "Matrix saved: $MATRIX"
