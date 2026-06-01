#!/bin/bash
# ==============================================================================
# Filter Optimization — 필터 수준별 Sensitivity/Precision 비교
# ==============================================================================
# 1) joint_raw.vcf.gz에서 다양한 필터로 SNP 추출
# 2) 각 샘플별 DP≥5 영역 한정 ground truth 생성
# 3) 전체 genome 대비 + DP≥5 한정 Sensitivity/Precision 계산
# ==============================================================================
set -euo pipefail
source ~/miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

REF="/home/yusanghyeon/RAT_project/PHMG_IT/reference/rn7.fa"
PILOT_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/pilot"
GT_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/ground_truth"
JOINT_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/joint_calling"
RAW_VCF="$JOINT_DIR/joint_raw.vcf.gz"
OUT_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/filter_optimization"
THREADS=14

SAMPLES=(C1 C2 C3 C4 C5 P1 P2 P3 P4 P5 P6 P7 P8 P9 P10)

mkdir -p "$OUT_DIR"

# ── Step 1: 필터 수준별 VCF 생성 ─────────────────────────────────────────
echo "================================================================"
echo "Step 1: Generating filtered VCFs at multiple thresholds"
echo "================================================================"

declare -A FILTER_EXPR
FILTER_EXPR[F1_current]='QUAL>=20 && INFO/DP>=15'
FILTER_EXPR[F2_moderate]='QUAL>=10 && INFO/DP>=10'
FILTER_EXPR[F3_lenient]='QUAL>=5 && INFO/DP>=5'
FILTER_EXPR[F4_minimal]='QUAL>=1'

FILTER_ORDER=(F1_current F2_moderate F3_lenient F4_minimal)

for FNAME in "${FILTER_ORDER[@]}"; do
    EXPR="${FILTER_EXPR[$FNAME]}"
    FILT_VCF="$OUT_DIR/${FNAME}_snps.vcf.gz"
    echo "  $FNAME: $EXPR"

    bcftools view -v snps "$RAW_VCF" \
    | bcftools filter -i "$EXPR" \
        -Oz -o "$FILT_VCF"
    bcftools index "$FILT_VCF"

    COUNT=$(bcftools view -H "$FILT_VCF" | wc -l)
    echo "    → $COUNT SNPs"
done

# ── Step 2: 샘플별 DP≥5 한정 ground truth 생성 ──────────────────────────
echo ""
echo "================================================================"
echo "Step 2: Generating DP>=5 restricted ground truth per sample"
echo "================================================================"

GT_RESTRICTED_DIR="$OUT_DIR/ground_truth_dp5"
mkdir -p "$GT_RESTRICTED_DIR"

for S in "${SAMPLES[@]}"; do
    GT_VCF="${GT_DIR}/${S}_dna_snps.vcf.gz"

    # Get DP≥5 BED
    if [ "$S" = "C1" ]; then
        DP5_BED="${PILOT_DIR}/arm_C_star_aggressive/arm_C_dp5.bed"
    else
        DP5_BED="${PILOT_DIR}/arm_C_${S}/arm_C_${S}_dp5.bed"
    fi

    GT_DP5="$GT_RESTRICTED_DIR/${S}_dna_snps_dp5.vcf.gz"

    bcftools view -R "$DP5_BED" "$GT_VCF" -Oz -o "$GT_DP5"
    bcftools index "$GT_DP5"

    COUNT=$(bcftools view -H "$GT_DP5" | wc -l)
    TOTAL=$(bcftools view -H "$GT_VCF" | wc -l)
    echo "  $S: $COUNT / $TOTAL DNA SNPs in DP≥5 regions ($(echo "scale=1; $COUNT * 100 / $TOTAL" | bc)%)"
done

# ── Step 3: 필터별 × 샘플별 평가 ─────────────────────────────────────────
echo ""
echo "================================================================"
echo "Step 3: Per-filter, per-sample evaluation"
echo "================================================================"

SUMMARY="$OUT_DIR/filter_optimization_summary.tsv"
echo -e "Filter\tSample\tScope\tRNA_calls\tDNA_SNPs\tTP\tFP\tFN\tSensitivity\tPrecision\tF1" > "$SUMMARY"

for FNAME in "${FILTER_ORDER[@]}"; do
    FILT_VCF="$OUT_DIR/${FNAME}_snps.vcf.gz"
    ISEC_BASE="$OUT_DIR/isec_${FNAME}"
    mkdir -p "$ISEC_BASE"

    echo ""
    echo "── Filter: $FNAME ──"

    for S in "${SAMPLES[@]}"; do
        GT_VCF="${GT_DIR}/${S}_dna_snps.vcf.gz"
        GT_DP5="$GT_RESTRICTED_DIR/${S}_dna_snps_dp5.vcf.gz"

        # Extract sample genotypes
        SAMPLE_VCF="$ISEC_BASE/${S}_calls.vcf.gz"
        bcftools view -s "$S" "$FILT_VCF" \
        | bcftools view -i 'GT="alt"' \
            -Oz -o "$SAMPLE_VCF"
        bcftools index "$SAMPLE_VCF"

        # --- Evaluation A: Whole genome ---
        ISEC_DIR="$ISEC_BASE/isec_${S}_whole"
        rm -rf "$ISEC_DIR"
        bcftools isec -p "$ISEC_DIR" "$SAMPLE_VCF" "$GT_VCF" -Oz

        TP=$(bcftools view -H "$ISEC_DIR/0002.vcf.gz" | wc -l)
        FP=$(bcftools view -H "$ISEC_DIR/0000.vcf.gz" | wc -l)
        FN=$(bcftools view -H "$ISEC_DIR/0001.vcf.gz" | wc -l)
        RNA_CALLS=$((TP + FP))
        DNA_SNPS=$(bcftools view -H "$GT_VCF" | wc -l)

        if [ "$TP" -gt 0 ]; then
            SENS=$(echo "scale=4; $TP / ($TP + $FN)" | bc)
            PREC=$(echo "scale=4; $TP / ($TP + $FP)" | bc)
            F1_SCORE=$(echo "scale=4; 2 * $SENS * $PREC / ($SENS + $PREC)" | bc)
        else
            SENS="0"; PREC="0"; F1_SCORE="0"
        fi
        echo -e "${FNAME}\t${S}\twhole_genome\t${RNA_CALLS}\t${DNA_SNPS}\t${TP}\t${FP}\t${FN}\t${SENS}\t${PREC}\t${F1_SCORE}" >> "$SUMMARY"

        # --- Evaluation B: DP≥5 restricted ---
        ISEC_DIR="$ISEC_BASE/isec_${S}_dp5"
        rm -rf "$ISEC_DIR"
        bcftools isec -p "$ISEC_DIR" "$SAMPLE_VCF" "$GT_DP5" -Oz

        TP5=$(bcftools view -H "$ISEC_DIR/0002.vcf.gz" | wc -l)
        FP5=$(bcftools view -H "$ISEC_DIR/0000.vcf.gz" | wc -l)
        FN5=$(bcftools view -H "$ISEC_DIR/0001.vcf.gz" | wc -l)
        RNA5=$((TP5 + FP5))
        DNA5=$(bcftools view -H "$GT_DP5" | wc -l)

        if [ "$TP5" -gt 0 ]; then
            SENS5=$(echo "scale=4; $TP5 / ($TP5 + $FN5)" | bc)
            PREC5=$(echo "scale=4; $TP5 / ($TP5 + $FP5)" | bc)
            F1_5=$(echo "scale=4; 2 * $SENS5 * $PREC5 / ($SENS5 + $PREC5)" | bc)
        else
            SENS5="0"; PREC5="0"; F1_5="0"
        fi
        echo -e "${FNAME}\t${S}\tdp5_restricted\t${RNA5}\t${DNA5}\t${TP5}\t${FP5}\t${FN5}\t${SENS5}\t${PREC5}\t${F1_5}" >> "$SUMMARY"

        echo "  $S: whole(Sens=${SENS},Prec=${PREC}) dp5(Sens=${SENS5},Prec=${PREC5})"
    done
done

# ── Step 4: Summary ──────────────────────────────────────────────────────
echo ""
echo "================================================================"
echo "Step 4: Summary by filter level"
echo "================================================================"

for FNAME in "${FILTER_ORDER[@]}"; do
    echo ""
    echo "── $FNAME ──"
    awk -F'\t' -v f="$FNAME" '
    $1==f && $3=="whole_genome" { n1++; s1+=$9; p1+=$10; f1+=$11; tp1+=$6; fp1+=$7; fn1+=$8 }
    $1==f && $3=="dp5_restricted" { n2++; s2+=$9; p2+=$10; f2+=$11; tp2+=$6; fp2+=$7; fn2+=$8 }
    END {
        if(n1>0) {
            os1=tp1/(tp1+fn1); op1=tp1/(tp1+fp1); of1=2*os1*op1/(os1+op1)
            printf "  Whole genome:   Avg Sens=%.4f  Avg Prec=%.4f  Avg F1=%.4f  | Overall Sens=%.4f  Prec=%.4f  F1=%.4f\n", s1/n1, p1/n1, f1/n1, os1, op1, of1
        }
        if(n2>0) {
            os2=tp2/(tp2+fn2); op2=tp2/(tp2+fp2); of2=2*os2*op2/(os2+op2)
            printf "  DP>=5 restricted: Avg Sens=%.4f  Avg Prec=%.4f  Avg F1=%.4f  | Overall Sens=%.4f  Prec=%.4f  F1=%.4f\n", s2/n2, p2/n2, f2/n2, os2, op2, of2
        }
    }' "$SUMMARY"
done

echo ""
echo "================================================================"
echo "DONE! Full results: $SUMMARY"
echo "================================================================"
