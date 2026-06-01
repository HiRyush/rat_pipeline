#!/bin/bash
set -euo pipefail

# Step 7: WES ground truth validation
# Mammary는 PHMG_IT과 달리:
#   - WES (exome only) — 우리도 exonic만 calling했으므로 영역 일치
#   - Per-tumor matched WES — 더 정밀한 per-sample 검증 가능
#
# 우리 final candidate set의 PPV 산출:
#   - True positive: WES에서도 confirmed (treatment-induced non-germline)
#   - False positive: WES에 없음 (likely RNA artifact or germline contamination)

WORK="/home/yusanghyeon/RAT_project/mammary_cancer"
GT_DIR="${WORK}/ground_truth"
OUT_DIR="${WORK}/results/observation_first"
VAL_DIR="${OUT_DIR}/validation"
mkdir -p "$VAL_DIR"

source /home/yusanghyeon/RAT_project/miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

# Late group의 WES VCF list (검증 기준)
# tumor_name → GSM 매핑 사용
declare -A LATE_WES
LATE_WES[M12_1]="GSM8994500_12m_1_tumor.pass.vcf"
LATE_WES[M12_2]="GSM8994501_12m_2_tumor.pass.vcf"
LATE_WES[M12_3]="GSM8994502_12m_3_tumor.pass.vcf"
LATE_WES[M20_1]="GSM8994503_20m_1_tumor.pass.vcf"
LATE_WES[M20_2]="GSM8994504_20m_2_tumor.pass.vcf"
LATE_WES[M20_3]="GSM8994505_20m_3_tumor.pass.vcf"
LATE_WES[M20_4]="GSM8994506_20m_4_tumor.pass.vcf"

echo "[$(date +%H:%M:%S)] Step 1: Late group WES union 생성..."

LATE_WES_UNION="${VAL_DIR}/late_wes_union.vcf.gz"
WES_LIST=()
for tumor in "${!LATE_WES[@]}"; do
  WES="${GT_DIR}/${LATE_WES[$tumor]}"
  if [ -f "$WES" ]; then
    # bgzip + index 안 됐을 가능성 → 처리
    GZ="${VAL_DIR}/$(basename "$WES").gz"
    if [ ! -f "$GZ" ]; then
      bgzip -c "$WES" > "$GZ"
      bcftools index -t "$GZ"
    fi
    WES_LIST+=("$GZ")
  fi
done

# WES union (any Late tumor's WES PASS)
bcftools merge --force-samples "${WES_LIST[@]}" 2>/dev/null \
  | bcftools view -i 'N_PASS(GT[*]="alt") >= 1' -Oz -o "$LATE_WES_UNION"
bcftools index -t "$LATE_WES_UNION"

WES_COUNT=$(bcftools view -H "$LATE_WES_UNION" | wc -l)
echo "  Late WES union (${#WES_LIST[@]} samples): ${WES_COUNT} variants"

# Early WES union도 (germline leak 추정용)
echo ""
echo "[$(date +%H:%M:%S)] Step 2: Early group WES union (germline leak 평가용)..."

declare -A EARLY_WES
EARLY_WES[M1_1]="GSM8994482_1m_1_tumor.pass.vcf"
EARLY_WES[M1_2]="GSM8994483_1m_2_tumor.pass.vcf"
EARLY_WES[M1_3]="GSM8994484_1m_3_tumor.pass.vcf"
EARLY_WES[M1_5]="GSM8994486_1m_5_tumor.pass.vcf"
EARLY_WES[M1_6]="GSM8994487_1m_6_tumor.pass.vcf"
EARLY_WES[M1_7]="GSM8994488_1m_7_tumor.pass.vcf"

EARLY_WES_UNION="${VAL_DIR}/early_wes_union.vcf.gz"
EWES_LIST=()
for tumor in "${!EARLY_WES[@]}"; do
  WES="${GT_DIR}/${EARLY_WES[$tumor]}"
  if [ -f "$WES" ]; then
    GZ="${VAL_DIR}/$(basename "$WES").gz"
    if [ ! -f "$GZ" ]; then
      bgzip -c "$WES" > "$GZ"
      bcftools index -t "$GZ"
    fi
    EWES_LIST+=("$GZ")
  fi
done

bcftools merge --force-samples "${EWES_LIST[@]}" 2>/dev/null \
  | bcftools view -i 'N_PASS(GT[*]="alt") >= 1' -Oz -o "$EARLY_WES_UNION"
bcftools index -t "$EARLY_WES_UNION"

EWES_COUNT=$(bcftools view -H "$EARLY_WES_UNION" | wc -l)
echo "  Early WES union (${#EWES_LIST[@]} samples): ${EWES_COUNT} variants"

# Validation
echo ""
echo "[$(date +%H:%M:%S)] Step 3: Candidate validation..."

CANDIDATES="${OUT_DIR}/step5_rna_editing_removed.vcf.gz"
CAND_COUNT=$(bcftools view -H "$CANDIDATES" | wc -l)
echo "  Total candidates: ${CAND_COUNT}"

# True positive: candidate ∩ Late WES (confirmed in late tumor DNA)
TP_VCF="${VAL_DIR}/true_positive.vcf.gz"
bcftools isec -p "${VAL_DIR}/tp_isec" -n=2 -w1 -Oz "$CANDIDATES" "$LATE_WES_UNION"
mv "${VAL_DIR}/tp_isec/0000.vcf.gz" "$TP_VCF"
mv "${VAL_DIR}/tp_isec/0000.vcf.gz.tbi" "${TP_VCF}.tbi" 2>/dev/null || bcftools index -t "$TP_VCF"
rm -rf "${VAL_DIR}/tp_isec"
TP=$(bcftools view -H "$TP_VCF" | wc -l)

# Germline leak: candidate ∩ Early WES (Early tumor에도 있으면 germline)
GL_VCF="${VAL_DIR}/germline_leak.vcf.gz"
bcftools isec -p "${VAL_DIR}/gl_isec" -n=2 -w1 -Oz "$CANDIDATES" "$EARLY_WES_UNION"
mv "${VAL_DIR}/gl_isec/0000.vcf.gz" "$GL_VCF"
mv "${VAL_DIR}/gl_isec/0000.vcf.gz.tbi" "${GL_VCF}.tbi" 2>/dev/null || bcftools index -t "$GL_VCF"
rm -rf "${VAL_DIR}/gl_isec"
GL=$(bcftools view -H "$GL_VCF" | wc -l)

# RNA artifact: candidate에 있지만 WES (Late+Early) 어디에도 없음
RNA_ART_VCF="${VAL_DIR}/rna_artifact.vcf.gz"
ANY_WES="${VAL_DIR}/any_wes.vcf.gz"
bcftools merge "$LATE_WES_UNION" "$EARLY_WES_UNION" --force-samples 2>/dev/null | bcftools view -Oz -o "$ANY_WES"
bcftools index -t "$ANY_WES"
bcftools isec -p "${VAL_DIR}/art_isec" -n=1 -w1 -Oz "$CANDIDATES" "$ANY_WES"
mv "${VAL_DIR}/art_isec/0000.vcf.gz" "$RNA_ART_VCF"
mv "${VAL_DIR}/art_isec/0000.vcf.gz.tbi" "${RNA_ART_VCF}.tbi" 2>/dev/null || bcftools index -t "$RNA_ART_VCF"
rm -rf "${VAL_DIR}/art_isec"
RNA_ART=$(bcftools view -H "$RNA_ART_VCF" | wc -l)

# True somatic: Late WES에 있지만 Early WES에 없음 (treatment-enriched)
TRUE_SOMATIC_VCF="${VAL_DIR}/true_somatic.vcf.gz"
bcftools isec -p "${VAL_DIR}/ts_isec" -n=1 -w1 -Oz "$TP_VCF" "$EARLY_WES_UNION"
mv "${VAL_DIR}/ts_isec/0000.vcf.gz" "$TRUE_SOMATIC_VCF"
mv "${VAL_DIR}/ts_isec/0000.vcf.gz.tbi" "${TRUE_SOMATIC_VCF}.tbi" 2>/dev/null || bcftools index -t "$TRUE_SOMATIC_VCF"
rm -rf "${VAL_DIR}/ts_isec"
TS=$(bcftools view -H "$TRUE_SOMATIC_VCF" | wc -l)

echo ""
echo "=== Mammary Validation Results ==="
echo "  Total candidates:           ${CAND_COUNT}"
echo "  WES-confirmed (Late):       ${TP}"
echo "    True somatic (Late ∖ Early): ${TS}"
echo "    Germline leak (Late ∩ Early): $((TP - TS))"
echo "  Germline leak (Early WES):  ${GL}"
echo "  RNA artifact (no WES):      ${RNA_ART}"
echo ""
PPV_TP=$(awk -v tp=$TP -v c=$CAND_COUNT 'BEGIN{printf "%.1f", (c>0)?tp/c*100:0}')
PPV_TS=$(awk -v ts=$TS -v c=$CAND_COUNT 'BEGIN{printf "%.1f", (c>0)?ts/c*100:0}')
RNA_RATE=$(awk -v r=$RNA_ART -v c=$CAND_COUNT 'BEGIN{printf "%.1f", (c>0)?r/c*100:0}')
GL_RATE=$(awk -v g=$GL -v c=$CAND_COUNT 'BEGIN{printf "%.1f", (c>0)?g/c*100:0}')

echo "  PPV (any WES confirm):    ${PPV_TP}%"
echo "  PPV (true somatic):       ${PPV_TS}%"
echo "  RNA artifact rate:        ${RNA_RATE}%"
echo "  Germline leak rate:       ${GL_RATE}%"

# 결과 저장
cat > "${VAL_DIR}/summary.tsv" <<EOF
metric	value
total_candidates	${CAND_COUNT}
wes_confirmed	${TP}
true_somatic	${TS}
germline_leak	${GL}
rna_artifact	${RNA_ART}
ppv_any_wes	${PPV_TP}
ppv_true_somatic	${PPV_TS}
rna_artifact_rate	${RNA_RATE}
germline_leak_rate	${GL_RATE}
EOF

echo ""
echo "Summary saved: ${VAL_DIR}/summary.tsv"
