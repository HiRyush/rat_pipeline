#!/bin/bash
# ==============================================================================
# FreeBayes evaluation — chr1 제외, 나머지 21개 염색체로 먼저 결과 확인
# ==============================================================================
set -euo pipefail
source ~/miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

REF="/home/yusanghyeon/RAT_project/PHMG_IT/reference/rn7.fa"
GT_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/ground_truth"
GT_DP5_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/filter_optimization/ground_truth_dp5"
OUT_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/freebayes"

SAMPLES=(C1 C2 C3 C4 C5 P1 P2 P3 P4 P5 P6 P7 P8 P9 P10)

# ── Step 1: Merge all completed chromosome VCFs (excluding chr1) ─────────
echo "================================================================"
echo "Step 1: Merge completed chromosome VCFs (excluding chr1)"
echo "================================================================"

RAW_VCF="$OUT_DIR/freebayes_no_chr1_raw.vcf.gz"
VCF_LIST=""

for CHR in chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chrX chrY scaffolds; do
    F="$OUT_DIR/per_chrom/${CHR}_raw.vcf"
    if [ -s "$F" ]; then
        GZ="${F}.gz"
        if [ ! -f "$GZ" ]; then
            bgzip -c "$F" > "$GZ"
            bcftools index "$GZ"
        fi
        VCF_LIST="$VCF_LIST $GZ"
    fi
done

bcftools concat $VCF_LIST -a -Oz -o "$RAW_VCF"
bcftools index "$RAW_VCF"

TOTAL_RAW=$(bcftools view -H "$RAW_VCF" | wc -l)
echo "  Raw variants (no chr1): $TOTAL_RAW"

# ── Step 2: Filter levels ────────────────────────────────────────────────
echo ""
echo "================================================================"
echo "Step 2: Filter SNPs at multiple thresholds"
echo "================================================================"

declare -A FILTER_EXPR
FILTER_EXPR[F1_current]='QUAL>=20 && INFO/DP>=15'
FILTER_EXPR[F3_lenient]='QUAL>=5 && INFO/DP>=5'
FILTER_EXPR[F4_minimal]='QUAL>=1'

FILTER_ORDER=(F1_current F3_lenient F4_minimal)

for FNAME in "${FILTER_ORDER[@]}"; do
    EXPR="${FILTER_EXPR[$FNAME]}"
    FILT_VCF="$OUT_DIR/no_chr1_${FNAME}_snps.vcf.gz"

    bcftools view -v snps "$RAW_VCF" \
    | bcftools filter -i "$EXPR" \
        -Oz -o "$FILT_VCF"
    bcftools index "$FILT_VCF"

    COUNT=$(bcftools view -H "$FILT_VCF" | wc -l)
    echo "  $FNAME ($EXPR): $COUNT SNPs"
done

# ── Step 3: Per-filter, per-sample evaluation ────────────────────────────
echo ""
echo "================================================================"
echo "Step 3: Per-filter, per-sample evaluation"
echo "================================================================"

SUMMARY="$OUT_DIR/freebayes_no_chr1_evaluation.tsv"
echo -e "Filter\tSample\tScope\tRNA_calls\tDNA_SNPs\tTP\tFP\tFN\tSensitivity\tPrecision\tF1" > "$SUMMARY"

for FNAME in "${FILTER_ORDER[@]}"; do
    FILT_VCF="$OUT_DIR/no_chr1_${FNAME}_snps.vcf.gz"
    ISEC_BASE="$OUT_DIR/isec_no_chr1_${FNAME}"
    mkdir -p "$ISEC_BASE"

    echo ""
    echo "── Filter: $FNAME ──"

    for S in "${SAMPLES[@]}"; do
        GT_VCF="${GT_DIR}/${S}_dna_snps.vcf.gz"
        GT_DP5="${GT_DP5_DIR}/${S}_dna_snps_dp5.vcf.gz"

        # Exclude chr1 from ground truth for fair comparison
        GT_NO1="$ISEC_BASE/${S}_gt_no_chr1.vcf.gz"
        GT_DP5_NO1="$ISEC_BASE/${S}_gt_dp5_no_chr1.vcf.gz"
        bcftools view -t ^chr1 "$GT_VCF" -Oz -o "$GT_NO1"
        bcftools index "$GT_NO1"
        bcftools view -t ^chr1 "$GT_DP5" -Oz -o "$GT_DP5_NO1"
        bcftools index "$GT_DP5_NO1"

        # Extract sample genotypes
        SAMPLE_VCF="$ISEC_BASE/${S}_calls.vcf.gz"
        bcftools view -s "$S" "$FILT_VCF" \
        | bcftools view -i 'GT="alt"' \
            -Oz -o "$SAMPLE_VCF"
        bcftools index "$SAMPLE_VCF"

        # --- Whole genome (no chr1) ---
        ISEC_DIR="$ISEC_BASE/isec_${S}_whole"
        rm -rf "$ISEC_DIR"
        bcftools isec -p "$ISEC_DIR" "$SAMPLE_VCF" "$GT_NO1" -Oz

        TP=$(bcftools view -H "$ISEC_DIR/0002.vcf.gz" | wc -l)
        FP=$(bcftools view -H "$ISEC_DIR/0000.vcf.gz" | wc -l)
        FN=$(bcftools view -H "$ISEC_DIR/0001.vcf.gz" | wc -l)
        RNA_CALLS=$((TP + FP))
        DNA_SNPS=$(bcftools view -H "$GT_NO1" | wc -l)

        if [ "$TP" -gt 0 ]; then
            SENS=$(echo "scale=4; $TP / ($TP + $FN)" | bc)
            PREC=$(echo "scale=4; $TP / ($TP + $FP)" | bc)
            F1_SCORE=$(echo "scale=4; 2 * $SENS * $PREC / ($SENS + $PREC)" | bc)
        else
            SENS="0"; PREC="0"; F1_SCORE="0"
        fi
        echo -e "${FNAME}\t${S}\twhole_genome\t${RNA_CALLS}\t${DNA_SNPS}\t${TP}\t${FP}\t${FN}\t${SENS}\t${PREC}\t${F1_SCORE}" >> "$SUMMARY"

        # --- DP≥5 restricted (no chr1) ---
        ISEC_DIR="$ISEC_BASE/isec_${S}_dp5"
        rm -rf "$ISEC_DIR"
        bcftools isec -p "$ISEC_DIR" "$SAMPLE_VCF" "$GT_DP5_NO1" -Oz

        TP5=$(bcftools view -H "$ISEC_DIR/0002.vcf.gz" | wc -l)
        FP5=$(bcftools view -H "$ISEC_DIR/0000.vcf.gz" | wc -l)
        FN5=$(bcftools view -H "$ISEC_DIR/0001.vcf.gz" | wc -l)

        if [ "$TP5" -gt 0 ]; then
            SENS5=$(echo "scale=4; $TP5 / ($TP5 + $FN5)" | bc)
            PREC5=$(echo "scale=4; $TP5 / ($TP5 + $FP5)" | bc)
            F1_5=$(echo "scale=4; 2 * $SENS5 * $PREC5 / ($SENS5 + $PREC5)" | bc)
        else
            SENS5="0"; PREC5="0"; F1_5="0"
        fi
        echo -e "${FNAME}\t${S}\tdp5_restricted\t${RNA_CALLS}\t${DNA_SNPS}\t${TP5}\t${FP5}\t${FN5}\t${SENS5}\t${PREC5}\t${F1_5}" >> "$SUMMARY"

        echo "  $S: whole(Sens=${SENS},Prec=${PREC}) dp5(Sens=${SENS5},Prec=${PREC5})"
    done
done

# ── Step 4: Summary ──────────────────────────────────────────────────────
echo ""
echo "================================================================"
echo "Step 4: Summary by filter level (excluding chr1)"
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
            printf "  Whole genome:     Avg Sens=%.4f  Avg Prec=%.4f  Avg F1=%.4f  | Overall Sens=%.4f  Prec=%.4f  F1=%.4f\n", s1/n1, p1/n1, f1/n1, os1, op1, of1
        }
        if(n2>0) {
            os2=tp2/(tp2+fn2); op2=tp2/(tp2+fp2); of2=2*os2*op2/(os2+op2)
            printf "  DP>=5 restricted: Avg Sens=%.4f  Avg Prec=%.4f  Avg F1=%.4f  | Overall Sens=%.4f  Prec=%.4f  F1=%.4f\n", s2/n2, p2/n2, f2/n2, os2, op2, of2
        }
    }' "$SUMMARY"
done

echo ""
echo "================================================================"
echo "DONE! Results: $SUMMARY"
echo "================================================================"
