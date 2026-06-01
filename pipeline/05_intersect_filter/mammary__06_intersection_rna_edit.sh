#!/bin/bash
set -euo pipefail

# Step 6: Imputed differential intersection + RNA editing filter
# PHMG_IT Step 4 + Step 5a 로직과 동일

WORK="/home/yusanghyeon/RAT_project/mammary_cancer"
FILT_DIR="${WORK}/results/mutect2/filtered"
DIFF_DIR="${WORK}/results/mutect2/differential"
IMPUTED_DIR="${WORK}/imputation/imputed"
OUT_DIR="${WORK}/results/observation_first"
mkdir -p "$OUT_DIR"

source /home/yusanghyeon/RAT_project/miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

EARLY=(SRR33625170 SRR33625169 SRR33625168 SRR33625167 SRR33625166 SRR33625165)
LATE=(SRR33625154 SRR33625153 SRR33625152 SRR33625151 SRR33625150 SRR33625149 SRR33625148)
CHRS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 X)

echo "[$(date +%H:%M:%S)] Step 1: Merge imputed VCFs (per group)..."

# Per-sample, per-chr imputed VCFs를 sample별로 concat (모든 chr)
merge_sample_imputed() {
  local SRR=$1
  local OUT="${OUT_DIR}/imputed_${SRR}.vcf.gz"
  if [ -f "$OUT" ]; then return; fi
  local FILES=()
  for chr in "${CHRS[@]}"; do
    local F="${IMPUTED_DIR}/${SRR}_${chr}.vcf.gz"
    [ -f "$F" ] && FILES+=("$F")
  done
  bcftools concat "${FILES[@]}" -Oz -o "$OUT" 2>/dev/null
  bcftools index -t "$OUT"
}

for SRR in "${EARLY[@]}" "${LATE[@]}"; do
  merge_sample_imputed "$SRR" &
done
wait

echo "[$(date +%H:%M:%S)] Step 2: Merge groups (Early imputed, Late imputed)..."

# Group-level merge: ALT genotype 있는 sample 수로 group differential 정의
EARLY_IMP="${OUT_DIR}/early_imputed_merged.vcf.gz"
LATE_IMP="${OUT_DIR}/late_imputed_merged.vcf.gz"

EARLY_IMP_FILES=()
for SRR in "${EARLY[@]}"; do
  EARLY_IMP_FILES+=("${OUT_DIR}/imputed_${SRR}.vcf.gz")
done
bcftools merge "${EARLY_IMP_FILES[@]}" -Oz -o "$EARLY_IMP" 2>/dev/null
bcftools index -t "$EARLY_IMP"

LATE_IMP_FILES=()
for SRR in "${LATE[@]}"; do
  LATE_IMP_FILES+=("${OUT_DIR}/imputed_${SRR}.vcf.gz")
done
bcftools merge "${LATE_IMP_FILES[@]}" -Oz -o "$LATE_IMP" 2>/dev/null
bcftools index -t "$LATE_IMP"

echo "[$(date +%H:%M:%S)] Step 3: Imputed differential (Late ALT-rich, Early ALT-poor)..."

# Late group에서 ALT allele 있는 sample 수 >= 50%, Early group에서는 <= 10%
LATE_IMP_DIFF="${OUT_DIR}/imputed_late_only.vcf.gz"
bcftools view "$LATE_IMP" -i 'COUNT(GT[*]="alt") >= 4' -Oz -o "${OUT_DIR}/late_imp_majority.vcf.gz" 2>/dev/null
bcftools index -t "${OUT_DIR}/late_imp_majority.vcf.gz"
bcftools view "$EARLY_IMP" -i 'COUNT(GT[*]="alt") <= 1' -Oz -o "${OUT_DIR}/early_imp_minority.vcf.gz" 2>/dev/null
bcftools index -t "${OUT_DIR}/early_imp_minority.vcf.gz"

# Late-rich AND Early-poor 교집합
bcftools isec -p "${OUT_DIR}/imp_isec" -n=2 -w1 -Oz \
  "${OUT_DIR}/late_imp_majority.vcf.gz" \
  "${OUT_DIR}/early_imp_minority.vcf.gz"
mv "${OUT_DIR}/imp_isec/0000.vcf.gz" "$LATE_IMP_DIFF"
mv "${OUT_DIR}/imp_isec/0000.vcf.gz.tbi" "${LATE_IMP_DIFF}.tbi" 2>/dev/null || bcftools index -t "$LATE_IMP_DIFF"
rm -rf "${OUT_DIR}/imp_isec"

IMP_DIFF_COUNT=$(bcftools view -H "$LATE_IMP_DIFF" | wc -l)
echo "  Imputed Late-only differential: ${IMP_DIFF_COUNT}"

echo ""
echo "[$(date +%H:%M:%S)] Step 4: Intersection (MuTect2 differential ∩ Imputed differential)..."

MUTECT_DIFF="${DIFF_DIR}/late_only_differential.vcf.gz"
INTERSECTION="${OUT_DIR}/step4_intersection.vcf.gz"

bcftools isec -p "${OUT_DIR}/final_isec" -n=2 -w1 -Oz \
  "$MUTECT_DIFF" \
  "$LATE_IMP_DIFF"
mv "${OUT_DIR}/final_isec/0000.vcf.gz" "$INTERSECTION"
mv "${OUT_DIR}/final_isec/0000.vcf.gz.tbi" "${INTERSECTION}.tbi" 2>/dev/null || bcftools index -t "$INTERSECTION"
rm -rf "${OUT_DIR}/final_isec"

INT_COUNT=$(bcftools view -H "$INTERSECTION" | wc -l)
echo "  Step 4 Intersection: ${INT_COUNT}"

echo ""
echo "[$(date +%H:%M:%S)] Step 5: RNA editing filter (A>G, T>C 제거)..."

RNA_EDIT_FILTERED="${OUT_DIR}/step5_rna_editing_removed.vcf.gz"
bcftools view "$INTERSECTION" -e '(REF="A" && ALT="G") || (REF="T" && ALT="C")' -Oz -o "$RNA_EDIT_FILTERED"
bcftools index -t "$RNA_EDIT_FILTERED"

FINAL_COUNT=$(bcftools view -H "$RNA_EDIT_FILTERED" | wc -l)
echo "  After RNA editing filter: ${FINAL_COUNT}"

echo ""
echo "=== Summary ==="
echo "  MuTect2 Late-only differential:  $(bcftools view -H "$MUTECT_DIFF" | wc -l)"
echo "  Imputed Late-only differential:   ${IMP_DIFF_COUNT}"
echo "  Step 4 Intersection:              ${INT_COUNT}"
echo "  Step 5 RNA editing removed:       ${FINAL_COUNT}"
echo ""
echo "Output: $RNA_EDIT_FILTERED"
