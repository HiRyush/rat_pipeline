#!/bin/bash
set -euo pipefail

# Step 5: Early (M1) vs Late (M12+M20) differential 추출
# PHMG_IT의 Phase 9 differential과 동일한 로직:
#   - Late group에서 보이지만 Early group에서 안 보이는 variant = "late-only"
#   - 이 set이 carcinogen-progression-induced somatic candidate

WORK="/home/yusanghyeon/RAT_project/mammary_cancer"
FILT_DIR="${WORK}/results/mutect2/filtered"
DIFF_DIR="${WORK}/results/mutect2/differential"
mkdir -p "$DIFF_DIR"

source /home/yusanghyeon/RAT_project/miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

# Sample 그룹
EARLY=(SRR33625170 SRR33625169 SRR33625168 SRR33625167 SRR33625166 SRR33625165)
LATE=(SRR33625154 SRR33625153 SRR33625152 SRR33625151 SRR33625150 SRR33625149 SRR33625148)

echo "[$(date +%H:%M:%S)] Step 1: Union variant sets..."

# Early union
EARLY_UNION="${DIFF_DIR}/early_union.vcf.gz"
EARLY_FILES=()
for SRR in "${EARLY[@]}"; do
  EARLY_FILES+=("${FILT_DIR}/${SRR}.pass.vcf.gz")
done
bcftools merge --force-samples "${EARLY_FILES[@]}" 2>/dev/null \
  | bcftools view -i 'N_PASS(GT[*]="alt") >= 1' -Oz -o "$EARLY_UNION"
bcftools index -t "$EARLY_UNION"
EARLY_COUNT=$(bcftools view -H "$EARLY_UNION" | wc -l)
echo "  Early union (${#EARLY[@]} samples): ${EARLY_COUNT} variants"

# Late union
LATE_UNION="${DIFF_DIR}/late_union.vcf.gz"
LATE_FILES=()
for SRR in "${LATE[@]}"; do
  LATE_FILES+=("${FILT_DIR}/${SRR}.pass.vcf.gz")
done
bcftools merge --force-samples "${LATE_FILES[@]}" 2>/dev/null \
  | bcftools view -i 'N_PASS(GT[*]="alt") >= 1' -Oz -o "$LATE_UNION"
bcftools index -t "$LATE_UNION"
LATE_COUNT=$(bcftools view -H "$LATE_UNION" | wc -l)
echo "  Late union (${#LATE[@]} samples): ${LATE_COUNT} variants"

echo ""
echo "[$(date +%H:%M:%S)] Step 2: Late-only differential (Late ∖ Early)..."

# Late에는 있고 Early에는 없는 위치 추출
LATE_ONLY="${DIFF_DIR}/late_only_differential.vcf.gz"
bcftools isec -p "${DIFF_DIR}/isec_tmp" -n=1 -w1 -Oz "$LATE_UNION" "$EARLY_UNION"
mv "${DIFF_DIR}/isec_tmp/0000.vcf.gz" "$LATE_ONLY"
mv "${DIFF_DIR}/isec_tmp/0000.vcf.gz.tbi" "${LATE_ONLY}.tbi" 2>/dev/null || bcftools index -t "$LATE_ONLY"
rm -rf "${DIFF_DIR}/isec_tmp"

DIFF_COUNT=$(bcftools view -H "$LATE_ONLY" | wc -l)
echo "  Late-only differential: ${DIFF_COUNT} variants"

echo ""
echo "[$(date +%H:%M:%S)] DONE"
echo "Output: $LATE_ONLY"
