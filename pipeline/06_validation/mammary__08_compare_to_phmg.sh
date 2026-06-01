#!/bin/bash
set -euo pipefail

# Step 8: PHMG_IT vs Mammary 성능 비교 (method generalizability)
# Both datasets에서 동일 method가 일관된 성능을 보이는지 평가

WORK="/home/yusanghyeon/RAT_project/mammary_cancer"
OUT="${WORK}/results/observation_first/comparison_to_phmg.md"

MAMMARY_SUM="${WORK}/results/observation_first/validation/summary.tsv"

# PHMG_IT 결과 (CLAUDE.md 기록)
PHMG_TOTAL=961
PHMG_TS=657
PHMG_PPV=68.4
PHMG_RNA=3.0
PHMG_GL=28.6
PHMG_COV=72.6

# Mammary 결과 추출
M_TOTAL=$(awk -F'\t' '$1=="total_candidates" {print $2}' "$MAMMARY_SUM")
M_TS=$(awk -F'\t' '$1=="true_somatic" {print $2}' "$MAMMARY_SUM")
M_PPV=$(awk -F'\t' '$1=="ppv_true_somatic" {print $2}' "$MAMMARY_SUM")
M_RNA=$(awk -F'\t' '$1=="rna_artifact_rate" {print $2}' "$MAMMARY_SUM")
M_GL=$(awk -F'\t' '$1=="germline_leak_rate" {print $2}' "$MAMMARY_SUM")
M_COV_AVG=$(awk 'NR>1 {sum+=$7; n++} END {if(n>0) printf "%.1f", sum/n}' "${WORK}/imputation/evaluation/imputation_results.tsv")

cat > "$OUT" <<EOF
# Method Generalizability: PHMG_IT vs Mammary Cancer

**Date:** $(date +%Y-%m-%d)
**Method:** Observation-first + Imputation-as-Filter (RNA-seq only)

## 비교 표

| Metric | PHMG_IT | Mammary | 일관성 |
|---|:---:|:---:|:---:|
| Design | Control vs Treatment (5C+10T) | Early M1 vs Late M12+M20 (6+7) | 다른 group structure |
| Reference | rn7 | rn6 | 다름 |
| HRDP panel | 75 strains | 48 strains | 다름 |
| Ground truth | WGS (genome-wide) | WES (exome only) | 다름 |
| Imputation coverage | ${PHMG_COV}% | ${M_COV_AVG}% (sample 평균) | 다름 (HRDP panel size 차이) |
| **Total candidates** | ${PHMG_TOTAL} | ${M_TOTAL} | — |
| **True somatic** | ${PHMG_TS} | ${M_TS} | — |
| **PPV** | ${PHMG_PPV}% | ${M_PPV}% | $(awk -v a=$PHMG_PPV -v b=$M_PPV 'BEGIN{d=a-b; if(d<0)d=-d; if(d<10) print "✅ 일관"; else if(d<20) print "△ 차이 있음"; else print "❌ 불일관"}') |
| **RNA artifact rate** | ${PHMG_RNA}% | ${M_RNA}% | $(awk -v a=$PHMG_RNA -v b=$M_RNA 'BEGIN{d=a-b; if(d<0)d=-d; if(d<5) print "✅ 일관"; else print "△ 차이"}') |
| **Germline leak rate** | ${PHMG_GL}% | ${M_GL}% | $(awk -v a=$PHMG_GL -v b=$M_GL 'BEGIN{d=a-b; if(d<0)d=-d; if(d<10) print "✅ 일관"; else print "△ 차이"}') |

## 해석

### Method generalizability 판정
- PPV 차이 < 10%p: Method가 두 dataset에서 일관된 성능 → **generalizable**
- PPV 차이 10-20%p: Dataset-dependent variation 있음 → method 한계 명시 필요
- PPV 차이 > 20%p: Method가 dataset-specific 가능성 → 재설계 검토

### Confounding factors (해석 시 고려)
1. Reference 차이 (rn6 vs rn7) — coordinate/annotation 호환성
2. HRDP panel 크기 (48 vs 75 strains) — imputation 정확도에 영향
3. Ground truth 차이 (WES vs WGS) — exome 영역에서만 평가 가능
4. Group design 차이 (control-treatment vs early-late) — 생물학적 의미 다름

### Method 자체의 강점/약점
- ✓ RNA editing filter 효과 일관성
- ? Imputation-as-filter의 효과 (HRDP panel 크기 의존성)
- ? Group differential 로직의 carcinogenesis time-course에서의 작동

## 다음 단계
EOF

# 일관성에 따른 다음 단계 추천
if awk -v a=$PHMG_PPV -v b=$M_PPV 'BEGIN{d=a-b; if(d<0)d=-d; exit !(d<10)}'; then
  cat >> "$OUT" <<EOF
- ✅ Method generalizable 입증 → Phase B (sensitivity 특성화) 진행
- Phase A (PoN ablation)로 두 dataset에서 일관된 개선 효과 검증
EOF
else
  cat >> "$OUT" <<EOF
- ⚠️ Dataset-dependent variation 감지 → 차이 원인 분석 필요
- HRDP panel 크기 영향 검토 / Group design 영향 검토
EOF
fi

echo "Comparison written: $OUT"
echo ""
cat "$OUT"
