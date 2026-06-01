#!/bin/bash
# ==============================================================================
# Per-sample variant calling + Differential variant analysis
# ==============================================================================
# 1) 각 샘플별 개별 bcftools mpileup+call → per-sample VCF
# 2) 동료 파이프라인으로 Patient vs Control differential analysis
# ==============================================================================
set -euo pipefail
source ~/miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

REF="/home/yusanghyeon/RAT_project/PHMG_IT/reference/rn7.fa"
PILOT_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/pilot"
CALLABLE="/home/yusanghyeon/RAT_project/PHMG_IT/results/multi_sample/callable_K10_of_15.bed"
COWORKER="/home/yusanghyeon/RAT_project/PHMG_IT/coworker"
OUT_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/differential"
THREADS=14

CONTROLS=(C1 C2 C3 C4 C5)
PATIENTS=(P1 P2 P3 P4 P5 P6 P7 P8 P9 P10)
ALL_SAMPLES=(${CONTROLS[@]} ${PATIENTS[@]})

mkdir -p "$OUT_DIR/per_sample_vcf"

# ── Step 1: Per-sample variant calling ───────────────────────────────────
echo "================================================================"
echo "Step 1: Per-sample variant calling (bcftools mpileup + call)"
echo "================================================================"

run_sample() {
    local S="$1"
    local REF="$2"
    local CALLABLE="$3"
    local PILOT_DIR="$4"
    local OUT_DIR="$5"

    if [ "$S" = "C1" ]; then
        BAM="${PILOT_DIR}/arm_C_star_aggressive/C1_Aligned.sortedByCoord.out.bam"
    else
        BAM="${PILOT_DIR}/arm_C_${S}/${S}_Aligned.sortedByCoord.out.bam"
    fi

    VCF="$OUT_DIR/per_sample_vcf/${S}_calls.vcf.gz"

    bcftools mpileup \
        -f "$REF" \
        "$BAM" \
        -R "$CALLABLE" \
        -q 10 \
        -Q 13 \
        -a FORMAT/DP,FORMAT/AD,INFO/AD \
        -d 10000 \
        -Ou \
    | bcftools call \
        -mv \
        -Oz -o "$VCF"

    bcftools index "$VCF"

    COUNT=$(bcftools view -H "$VCF" | wc -l)
    echo "  $S: $COUNT variants"
}

export -f run_sample

RUNNING=0
for S in "${ALL_SAMPLES[@]}"; do
    run_sample "$S" "$REF" "$CALLABLE" "$PILOT_DIR" "$OUT_DIR" &
    RUNNING=$((RUNNING + 1))
    if [ "$RUNNING" -ge "$THREADS" ]; then
        wait -n
        RUNNING=$((RUNNING - 1))
    fi
done
wait
echo "  All samples completed."

# ── Step 2: Run differential variant analysis ────────────────────────────
echo ""
echo "================================================================"
echo "Step 2: Differential variant analysis (Patient vs Control)"
echo "================================================================"

# Build target and reference VCF arguments
TARGET_ARGS=""
for S in "${PATIENTS[@]}"; do
    TARGET_ARGS="$TARGET_ARGS $OUT_DIR/per_sample_vcf/${S}_calls.vcf.gz"
done

REF_ARGS=""
for S in "${CONTROLS[@]}"; do
    REF_ARGS="$REF_ARGS $OUT_DIR/per_sample_vcf/${S}_calls.vcf.gz"
done

cd "$COWORKER"

python run_pipeline.py \
    --mode independent \
    --target $TARGET_ARGS \
    --reference $REF_ARGS \
    --min-coverage 5 \
    --min-qual 5 \
    --min-alt-ratio 0.2 \
    --max-ref-alt-ratio 0.1 \
    --min-recurrence 2 \
    --min-delta-alt-ratio 0.15 \
    --output "$OUT_DIR/phmg_vs_control"

echo ""
echo "================================================================"
echo "DONE! Results in: $OUT_DIR/"
echo "================================================================"
