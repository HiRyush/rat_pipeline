# ==============================================================================
# STAR Alignment — Aggressive Parameters for Coverage Maximization
# ==============================================================================
# Based on: star_alignment_unsorted.sh
#
# Key changes from original:
#   1. --outFilterMultimapNmax     20 → 50   (multimap 허용 확대)
#   2. --outFilterScoreMinOverLread 0.5 → 0.33 (score 임계값 완화)
#   3. --outFilterMatchNminOverLread 0.5 → 0.33 (match 임계값 완화)
#   4. --winAnchorMultimapNmax     50 → 100  (anchor multimap 허용 확대)
#   5. BAM SortedByCoordinate 직접 출력 (별도 samtools sort 불필요)
#   6. Read Group 정보 STAR에서 직접 추가 (별도 AddOrReplaceReadGroups 불필요)
#
# Rationale:
#   - Read 수 자체는 +1.7%만 증가하나, 기존 read의 alignment 활용도가 크게 향상
#   - DP≥5 coverage: 41.3% → 45.6% (+4.3%p, C1 기준)
#   - DNA SNP DP≥5: 38.9% → 43.0% (+4.1%p)
#
# Important:
#   - 후속 단계에서 SplitNCigarReads를 생략함
#   - SplitNCigarReads 적용 시 junction overhang hard-clip으로
#     DP≥5 coverage가 54% 감소 (1,094 → 499 Mb)
#   - 이 BAM을 그대로 variant calling에 사용
#
# Species: Rat (Rattus norvegicus), Reference: UCSC rn7 (mRatBN7.2)
# Dependencies: STAR 2.7+, samtools
# ==============================================================================
set -euo pipefail

# ── Configuration ──────────────────────────────────────────────────────────
GENOMEDIR="reference/star_index"
FASTQDIR="fastq"
OUTDIR="result"
THREADS=14

# Sample list (adjust as needed)
SAMPLES=("C1" "C2" "C3" "C4" "C5" "P1" "P2" "P3" "P4" "P5" "P6" "P7" "P8" "P9" "P10")

mkdir -p ${OUTDIR}

for SAMPLE in "${SAMPLES[@]}"; do

    echo "=========================================="
    echo "STAR alignment (aggressive): ${SAMPLE}"
    echo "=========================================="

    BAM="${OUTDIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
    if [ -f "$BAM" ] && [ -f "${BAM}.bai" ]; then
        echo "  BAM already exists, skipping: ${SAMPLE}"
        continue
    fi

    # ── STAR Alignment ─────────────────────────────────────────────────
    #
    # Parameters unchanged from original:
    #   --twopassMode Basic
    #   --alignEndsType Local
    #   --alignSoftClipAtReferenceEnds Yes
    #   --outFilterMismatchNoverLmax 0.1
    #   --outFilterMismatchNmax 10
    #   --outMultimapperOrder Random
    #   --alignIntronMin 20
    #   --alignIntronMax 1000000
    #   --alignMatesGapMax 1000000
    #
    # Parameters CHANGED (← original value):
    #   --outFilterMultimapNmax      50   ← 20
    #   --outFilterScoreMinOverLread 0.33 ← 0.5
    #   --outFilterMatchNminOverLread 0.33 ← 0.5
    #   --winAnchorMultimapNmax      100  ← 50 (default)
    #
    STAR \
        --genomeDir ${GENOMEDIR} \
        --readFilesIn ${FASTQDIR}/${SAMPLE}_1.fastq.gz ${FASTQDIR}/${SAMPLE}_2.fastq.gz \
        --readFilesCommand zcat \
        --runThreadN ${THREADS} \
        --twopassMode Basic \
        --alignEndsType Local \
        --alignSoftClipAtReferenceEnds Yes \
        --outFilterMismatchNoverLmax 0.1 \
        --outFilterMismatchNmax 10 \
        --outFilterMultimapNmax 50 \
        --outMultimapperOrder Random \
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
        --outSAMattrRGline ID:${SAMPLE} SM:${SAMPLE} PL:ILLUMINA LB:lib1 \
        --outFileNamePrefix ${OUTDIR}/${SAMPLE}_ \
        --limitBAMsortRAM 50000000000

    # ── Index ──────────────────────────────────────────────────────────
    echo "  Indexing BAM..."
    samtools index "$BAM"

    # ── Cleanup ────────────────────────────────────────────────────────
    rm -rf ${OUTDIR}/${SAMPLE}__STARgenome ${OUTDIR}/${SAMPLE}__STARpass1 2>/dev/null

    echo "Completed: ${SAMPLE}"
    echo ""
done

echo "=========================================="
echo "All samples aligned!"
echo "=========================================="
echo ""
echo "Parameter changes vs original (star_alignment_unsorted.sh):"
echo "  outFilterMultimapNmax:      20 → 50"
echo "  outFilterScoreMinOverLread: 0.5 → 0.33"
echo "  outFilterMatchNminOverLread: 0.5 → 0.33"
echo "  winAnchorMultimapNmax:      50 → 100"
echo ""
echo "Coverage improvement (C1, single sample):"
echo "  DP>=5: 41.3% → 45.6% (+4.3%p)"
echo "  DNA SNP DP>=5: 38.9% → 43.0% (+4.1%p)"
echo "=========================================="
