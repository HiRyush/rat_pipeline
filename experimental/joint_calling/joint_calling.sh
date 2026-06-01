#!/bin/bash
# ==============================================================================
# Joint Variant Calling — 15 samples with bcftools mpileup
# ==============================================================================
# 15개 BAM을 동시 투입하여 position별 reads를 pooling
# callable region (K>=10) 내에서 calling하여 속도 + 의미 있는 영역에 집중
# ==============================================================================
set -euo pipefail
source ~/miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

REF="/home/yusanghyeon/RAT_project/PHMG_IT/reference/rn7.fa"
PILOT_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/pilot"
CALLABLE="/home/yusanghyeon/RAT_project/PHMG_IT/results/multi_sample/callable_K10_of_15.bed"
GT_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/ground_truth"
OUT_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/joint_calling"
THREADS=14

mkdir -p "$OUT_DIR"

# ── Step 1: BAM list ──────────────────────────────────────────────────────
echo "================================================================"
echo "Step 1: Preparing BAM list"
echo "================================================================"

BAM_LIST="$OUT_DIR/bam_list.txt"
> "$BAM_LIST"

for S in C1 C2 C3 C4 C5 P1 P2 P3 P4 P5 P6 P7 P8 P9 P10; do
    if [ "$S" = "C1" ]; then
        BAM="${PILOT_DIR}/arm_C_star_aggressive/C1_Aligned.sortedByCoord.out.bam"
    else
        BAM="${PILOT_DIR}/arm_C_${S}/${S}_Aligned.sortedByCoord.out.bam"
    fi
    echo "$BAM" >> "$BAM_LIST"
    echo "  $S: $BAM"
done
echo "  Total: $(wc -l < "$BAM_LIST") BAMs"

# ── Step 2: Joint calling with bcftools ───────────────────────────────────
echo ""
echo "================================================================"
echo "Step 2: bcftools mpileup + call (joint, 15 samples)"
echo "================================================================"
echo "  Target region: K>=10 callable ($(wc -l < "$CALLABLE") intervals)"
echo ""

RAW_VCF="$OUT_DIR/joint_raw.vcf.gz"

bcftools mpileup \
    -f "$REF" \
    -b "$BAM_LIST" \
    -R "$CALLABLE" \
    --threads "$THREADS" \
    -q 10 \
    -Q 13 \
    -a FORMAT/DP,FORMAT/AD,INFO/AD \
    -d 10000 \
    -Ou \
| bcftools call \
    -mv \
    --threads "$THREADS" \
    -Oz -o "$RAW_VCF"

bcftools index "$RAW_VCF"

TOTAL_RAW=$(bcftools view -H "$RAW_VCF" | wc -l)
echo "  Raw variants: $TOTAL_RAW"

# ── Step 3: Filter — SNPs only, QUAL>=20 ─────────────────────────────────
echo ""
echo "================================================================"
echo "Step 3: Filter variants"
echo "================================================================"

FILT_VCF="$OUT_DIR/joint_snps_filtered.vcf.gz"

bcftools view -v snps "$RAW_VCF" \
| bcftools filter \
    -i 'QUAL>=20 && INFO/DP>=15' \
    -Oz -o "$FILT_VCF"

bcftools index "$FILT_VCF"

TOTAL_FILT=$(bcftools view -H "$FILT_VCF" | wc -l)
echo "  Filtered SNPs (QUAL>=20, DP>=15): $TOTAL_FILT"

# ── Step 4: Per-sample VCF extraction ─────────────────────────────────────
echo ""
echo "================================================================"
echo "Step 4: Extract per-sample genotypes for ground truth comparison"
echo "================================================================"

for S in C1 C2 C3 C4 C5 P1 P2 P3 P4 P5 P6 P7 P8 P9 P10; do
    GT_VCF="${GT_DIR}/${S}_dna_snps.vcf.gz"
    [ -f "$GT_VCF" ] || continue

    SAMPLE_VCF="$OUT_DIR/${S}_joint_calls.vcf.gz"

    # Extract this sample's non-ref genotypes
    bcftools view -s "$S" "$FILT_VCF" \
    | bcftools view -i 'GT="alt"' \
        -Oz -o "$SAMPLE_VCF"
    bcftools index "$SAMPLE_VCF"
done

# ── Step 5: Ground truth comparison ───────────────────────────────────────
echo ""
echo "================================================================"
echo "Step 5: Ground truth comparison (DNA WGS vs RNA joint calling)"
echo "================================================================"
echo ""

EVAL_FILE="$OUT_DIR/joint_calling_evaluation.tsv"
echo -e "Sample\tRNA_calls\tDNA_SNPs\tTrue_Pos\tFalse_Pos\tFalse_Neg\tSensitivity\tPrecision\tF1" > "$EVAL_FILE"

for S in C1 C2 C3 C4 C5 P1 P2 P3 P4 P5 P6 P7 P8 P9 P10; do
    GT_VCF="${GT_DIR}/${S}_dna_snps.vcf.gz"
    [ -f "$GT_VCF" ] || continue

    SAMPLE_VCF="$OUT_DIR/${S}_joint_calls.vcf.gz"
    ISEC_DIR="$OUT_DIR/isec_${S}"
    rm -rf "$ISEC_DIR"

    # bcftools isec: find intersections
    bcftools isec \
        -p "$ISEC_DIR" \
        "$SAMPLE_VCF" \
        "$GT_VCF" \
        -Oz

    # 0000.vcf.gz = RNA only (false positives relative to DNA)
    # 0001.vcf.gz = DNA only (false negatives)
    # 0002.vcf.gz = shared, from RNA (true positives)
    # 0003.vcf.gz = shared, from DNA

    FP=$(bcftools view -H "$ISEC_DIR/0000.vcf.gz" | wc -l)
    FN=$(bcftools view -H "$ISEC_DIR/0001.vcf.gz" | wc -l)
    TP=$(bcftools view -H "$ISEC_DIR/0002.vcf.gz" | wc -l)

    RNA_CALLS=$((TP + FP))
    DNA_SNPS=$(bcftools view -H "$GT_VCF" | wc -l)

    if [ "$RNA_CALLS" -gt 0 ] && [ "$TP" -gt 0 ]; then
        SENS=$(echo "scale=4; $TP / ($TP + $FN)" | bc)
        PREC=$(echo "scale=4; $TP / ($TP + $FP)" | bc)
        F1=$(echo "scale=4; 2 * $SENS * $PREC / ($SENS + $PREC)" | bc)
    else
        SENS="0"; PREC="0"; F1="0"
    fi

    echo -e "${S}\t${RNA_CALLS}\t${DNA_SNPS}\t${TP}\t${FP}\t${FN}\t${SENS}\t${PREC}\t${F1}"
    echo -e "${S}\t${RNA_CALLS}\t${DNA_SNPS}\t${TP}\t${FP}\t${FN}\t${SENS}\t${PREC}\t${F1}" >> "$EVAL_FILE"
done

echo ""
echo "================================================================"
echo "Step 6: Summary statistics"
echo "================================================================"

# Compute overall averages
awk -F'\t' 'NR>1 {
    n++; sens+=$7; prec+=$8; f1+=$9; tp+=$4; fp+=$5; fn+=$6
} END {
    print "Samples evaluated: " n
    print "Average Sensitivity: " sens/n
    print "Average Precision:   " prec/n
    print "Average F1:          " f1/n
    print ""
    print "Total TP: " tp
    print "Total FP: " fp
    print "Total FN: " fn
    overall_sens = tp/(tp+fn)
    overall_prec = tp/(tp+fp)
    overall_f1 = 2*overall_sens*overall_prec/(overall_sens+overall_prec)
    print "Overall Sensitivity: " overall_sens
    print "Overall Precision:   " overall_prec
    print "Overall F1:          " overall_f1
}' "$EVAL_FILE"

echo ""
echo "================================================================"
echo "DONE! Results in: $OUT_DIR/"
echo "  - joint_raw.vcf.gz          (raw joint calls)"
echo "  - joint_snps_filtered.vcf.gz (filtered SNPs)"
echo "  - joint_calling_evaluation.tsv (per-sample accuracy)"
echo "================================================================"
