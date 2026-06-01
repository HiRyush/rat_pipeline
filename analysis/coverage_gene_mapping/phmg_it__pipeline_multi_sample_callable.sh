#!/bin/bash
###############################################################################
# Multi-sample Callable Region Pipeline for RNA-seq Variant Calling
#
# Purpose:
#   단일 샘플 RNA-seq의 DP≥5 coverage는 ~40-45%로 한계가 있으나,
#   여러 샘플에서 반복 관찰되는 DP≥1 영역을 "callable"로 정의하면
#   70%+ coverage를 달성할 수 있다.
#
# Method:
#   1) STAR aggressive alignment → sorted BAM (per sample)
#   2) bedtools genomecov → DP≥1 BED (per sample)
#   3) bedtools unionbedg → 15개 샘플 합산
#   4) "N개 중 K개 샘플에서 DP≥1" 기준으로 callable region 정의
#   5) DNA WGS ground truth 대비 SNP coverage 검증
#
# Species: Rat (Rattus norvegicus), Reference: UCSC rn7 (mRatBN7.2)
# Samples: 15 (C1-C5 Control, P1-P10 PHMG-treated)
#
# Dependencies: STAR 2.7+, samtools, bedtools, bcftools
#
# Usage:
#   bash pipeline_multi_sample_callable.sh
#
# Results (K=10, i.e. 10/15 samples DP≥1):
#   - Callable region: 1,990 Mb (75.1% genome)
#   - DNA SNP coverage: 72.4%
###############################################################################
set -euo pipefail

# ── Configuration ──────────────────────────────────────────────────────────
FASTQ_DIR="/media/yusanghyeon/30B4E366B4E32D52/01_Projects/korea_PHMG/RNA_fastq"
DNA_VCF_DIR="/media/yusanghyeon/30B4E366B4E32D52/01_Projects/korea_PHMG/DNA_vcf"
REF="/home/yusanghyeon/RAT_project/PHMG_IT/reference/rn7.fa"
STAR_INDEX="/home/yusanghyeon/RAT_project/PHMG_IT/reference/star_index"
RESULT_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/pilot"
GT_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/ground_truth"
OUT_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/multi_sample"
THREADS=14
GENOME_SIZE=2648062885  # rn7 genome size (bp)

SAMPLES=("C1" "C2" "C3" "C4" "C5" "P1" "P2" "P3" "P4" "P5" "P6" "P7" "P8" "P9" "P10")
K_THRESHOLDS=(1 2 3 5 8 10 12 15)

# DNA sample ↔ RNA sample mapping (for ground truth extraction)
declare -A DNA_MAP=(
    [C1]="P5S40W-2"   [C2]="P5S54W-2"  [C3]="P5S54W-3"
    [C4]="P5S54W-6"   [C5]="P5S54W-8"
    [P1]="P5H35W-11"  [P2]="P5H40W-2"  [P3]="P5H40W-16"
    [P4]="P5H48W-12"  [P5]="P5H49W-9"  [P6]="P5H49W-18"
    [P7]="P5H52W-8"   [P8]="P5H54W-2"  [P9]="P5H54W-19"
    [P10]="P5M44W-11"
)

mkdir -p "$RESULT_DIR" "$GT_DIR" "$OUT_DIR"

###############################################################################
# STEP 1: Ground Truth Preparation (DNA gVCF → SNP-only VCF)
###############################################################################
echo "================================================================"
echo "STEP 1: Extract SNPs from DNA gVCF (ground truth)"
echo "================================================================"

for S in "${SAMPLES[@]}"; do
    GT_OUT="${GT_DIR}/${S}_dna_snps.vcf.gz"
    if [ -f "$GT_OUT" ]; then
        echo "  $S: already exists, skip"
        continue
    fi
    DNA_VCF="${DNA_VCF_DIR}/${DNA_MAP[$S]}_sorted.genome.vcf.gz"
    echo "  $S (${DNA_MAP[$S]}): extracting SNPs..."
    bcftools view -v snps \
        -i 'QUAL>0 && (FILTER="PASS" || FILTER=".")' \
        "$DNA_VCF" -Oz -o "$GT_OUT"
    bcftools index "$GT_OUT"
done

###############################################################################
# STEP 2: STAR Aggressive Alignment (per sample)
###############################################################################
echo ""
echo "================================================================"
echo "STEP 2: STAR aggressive alignment"
echo "================================================================"
echo ""
echo "Parameters:"
echo "  --alignEndsType Local  (soft-clip 허용)"
echo "  --outFilterMultimapNmax 50"
echo "  --outFilterScoreMinOverLread 0.33"
echo "  --winAnchorMultimapNmax 100"
echo "  --twopassMode Basic"
echo "  * SplitNCigarReads 생략 (coverage 54% 손실 방지)"
echo ""

TOTAL=${#SAMPLES[@]}
for i in "${!SAMPLES[@]}"; do
    S="${SAMPLES[$i]}"
    N=$((i+1))

    # C1 uses a different output directory name (legacy)
    if [ "$S" = "C1" ]; then
        OUTDIR="${RESULT_DIR}/arm_C_star_aggressive"
    else
        OUTDIR="${RESULT_DIR}/arm_C_${S}"
    fi
    mkdir -p "$OUTDIR"

    BAM="${OUTDIR}/${S}_Aligned.sortedByCoord.out.bam"
    if [ -f "$BAM" ] && [ -f "${BAM}.bai" ]; then
        echo "[$N/$TOTAL] $S: BAM exists, skip"
        continue
    fi

    echo "[$N/$TOTAL] $S: running STAR..."
    STAR \
      --genomeDir "$STAR_INDEX" \
      --readFilesIn "${FASTQ_DIR}/${S}_1.fastq.gz" "${FASTQ_DIR}/${S}_2.fastq.gz" \
      --readFilesCommand zcat \
      --runThreadN "$THREADS" \
      --twopassMode Basic \
      --alignEndsType Local \
      --alignSoftClipAtReferenceEnds Yes \
      --outFilterMultimapNmax 50 \
      --outMultimapperOrder Random \
      --outFilterMismatchNoverLmax 0.1 \
      --outFilterMismatchNmax 10 \
      --outFilterScoreMinOverLread 0.33 \
      --outFilterMatchNminOverLread 0.33 \
      --winAnchorMultimapNmax 100 \
      --alignIntronMin 20 \
      --alignIntronMax 1000000 \
      --alignMatesGapMax 1000000 \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMstrandField intronMotif \
      --outSAMattributes NH HI AS nM XS MD \
      --outSAMunmapped Within \
      --outSAMattrRGline ID:${S} SM:${S} PL:ILLUMINA LB:lib1 \
      --outFileNamePrefix "${OUTDIR}/${S}_" \
      --limitBAMsortRAM 50000000000

    samtools index "$BAM"
    rm -rf "${OUTDIR}/${S}__STARgenome" "${OUTDIR}/${S}__STARpass1" 2>/dev/null
    echo "[$N/$TOTAL] $S: done"
done

###############################################################################
# STEP 3: Per-sample Coverage (DP≥1 BED)
###############################################################################
echo ""
echo "================================================================"
echo "STEP 3: Per-sample DP≥1 coverage BED"
echo "================================================================"

for S in "${SAMPLES[@]}"; do
    if [ "$S" = "C1" ]; then
        OUTDIR="${RESULT_DIR}/arm_C_star_aggressive"
        LABEL="arm_C"
    else
        OUTDIR="${RESULT_DIR}/arm_C_${S}"
        LABEL="arm_C_${S}"
    fi

    BED="${OUTDIR}/${LABEL}_dp1.bed"
    if [ -f "$BED" ]; then
        echo "  $S: DP≥1 BED exists, skip"
        continue
    fi

    BAM="${OUTDIR}/${S}_Aligned.sortedByCoord.out.bam"
    echo "  $S: computing genomecov..."
    bedtools genomecov -ibam "$BAM" -bg \
        | awk '$4 >= 1' \
        | bedtools merge -i - > "$BED"
    echo "  $S: done"
done

###############################################################################
# STEP 4: Multi-sample Callable Region
#   bedtools unionbedg로 15개 샘플의 DP≥1 coverage를 합산한 뒤,
#   K개 이상 샘플에서 관찰된 영역을 callable로 정의
###############################################################################
echo ""
echo "================================================================"
echo "STEP 4: Multi-sample callable region (unionbedg)"
echo "================================================================"

# Genome file for bedtools
GENOME_FILE="$OUT_DIR/rn7.genome"
awk '{print $1"\t"$2}' "${REF}.fai" > "$GENOME_FILE"

# Prepare per-sample bedGraph (value=1 where DP≥1)
BED_LIST=()
SAMPLE_LABELS=()
for S in "${SAMPLES[@]}"; do
    if [ "$S" = "C1" ]; then
        BED="${RESULT_DIR}/arm_C_star_aggressive/arm_C_dp1.bed"
    else
        BED="${RESULT_DIR}/arm_C_${S}/arm_C_${S}_dp1.bed"
    fi
    SORTED_BG="$OUT_DIR/${S}_dp1_sorted.bg"
    awk '{print $1"\t"$2"\t"$3"\t1"}' "$BED" | sort -k1,1 -k2,2n > "$SORTED_BG"
    BED_LIST+=("$SORTED_BG")
    SAMPLE_LABELS+=("$S")
done

# Union across all samples
UNION_BG="$OUT_DIR/union_15samples.bg"
echo "  Running unionbedg (15 samples)..."
bedtools unionbedg \
    -i "${BED_LIST[@]}" \
    -g "$GENOME_FILE" \
    -filler 0 \
    -header -names "${SAMPLE_LABELS[@]}" > "$UNION_BG"
echo "  Union regions: $(wc -l < "$UNION_BG")"

# Generate callable BED at each K threshold
CALLABLE_SUMMARY="$OUT_DIR/multi_sample_callable_summary.tsv"
echo -e "K\tCallable_Mb\tGenome_Pct" > "$CALLABLE_SUMMARY"

echo ""
echo "  K threshold | Callable (Mb) | Genome %"
echo "  ------------|---------------|--------"

for K in "${K_THRESHOLDS[@]}"; do
    CALLABLE_BED="$OUT_DIR/callable_K${K}_of_15.bed"
    awk -v k="$K" 'NR>1 {
        sum=0; for(i=4;i<=NF;i++) sum+=$i
        if(sum>=k) print $1"\t"$2"\t"$3
    }' "$UNION_BG" | bedtools merge -i - > "$CALLABLE_BED"

    BP=$(awk '{s+=$3-$2} END{print s+0}' "$CALLABLE_BED")
    MB=$(echo "scale=1; $BP/1000000" | bc)
    PCT=$(echo "scale=1; $BP*100/$GENOME_SIZE" | bc)

    printf "  K>=%-2d       | %'10s Mb | %s%%\n" "$K" "$MB" "$PCT"
    echo -e "${K}\t${MB}\t${PCT}" >> "$CALLABLE_SUMMARY"
done

# Cleanup temp bedGraph files
rm -f "$OUT_DIR"/*_dp1_sorted.bg

###############################################################################
# STEP 5: Validate with DNA WGS Ground Truth
#   각 K 기준의 callable region 내 DNA SNP 비율 측정
###############################################################################
echo ""
echo "================================================================"
echo "STEP 5: DNA SNP coverage validation"
echo "================================================================"

SNP_SUMMARY="$OUT_DIR/multi_sample_snp_coverage.tsv"
echo -e "K\tSample\tSNPs_in_callable\tTotal_SNPs\tPct" > "$SNP_SUMMARY"

for K in "${K_THRESHOLDS[@]}"; do
    CALLABLE_BED="$OUT_DIR/callable_K${K}_of_15.bed"
    TOTAL_SUM=0
    COVERED_SUM=0

    for S in "${SAMPLES[@]}"; do
        GT_VCF="${GT_DIR}/${S}_dna_snps.vcf.gz"
        [ -f "$GT_VCF" ] || continue

        TOTAL=$(bcftools view -H "$GT_VCF" | wc -l)
        COVERED=$(bcftools view -R "$CALLABLE_BED" -H "$GT_VCF" | wc -l)
        PCT=$(echo "scale=1; $COVERED*100/$TOTAL" | bc)

        echo -e "${K}\t${S}\t${COVERED}\t${TOTAL}\t${PCT}" >> "$SNP_SUMMARY"
        TOTAL_SUM=$((TOTAL_SUM + TOTAL))
        COVERED_SUM=$((COVERED_SUM + COVERED))
    done

    AVG=$(echo "scale=1; $COVERED_SUM*100/$TOTAL_SUM" | bc)
    echo "  K>=${K}: DNA SNP coverage = ${AVG}% (avg across ${#SAMPLES[@]} samples)"
done

###############################################################################
# Summary
###############################################################################
echo ""
echo "================================================================"
echo "DONE"
echo "================================================================"
echo ""
echo "Output files:"
echo "  $CALLABLE_SUMMARY"
echo "  $SNP_SUMMARY"
echo "  ${OUT_DIR}/callable_K{1,2,3,5,8,10,12,15}_of_15.bed"
echo ""
echo "Recommended threshold: K>=10 (10/15 samples DP>=1)"
echo "  → Callable: 1,990 Mb (75.1% genome)"
echo "  → DNA SNP coverage: 72.4%"
