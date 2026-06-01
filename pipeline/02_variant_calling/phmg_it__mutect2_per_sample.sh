#!/bin/bash
set -euo pipefail

# MuTect2 per-sample variant calling for PHMG_IT
#
# Pipeline:
#   arm_C BAM → MarkDuplicates → MuTect2 (tumor-only) → FilterMutectCalls → PASS VCF
#
# Why MuTect2 over bcftools (Phase 7):
#   - Somatic caller: detects low allele frequency variants bcftools misses
#   - RNA-seq allelic imbalance handled as somatic signal, not noise
#   - coworker vcf_parser.py는 GATK AD/DP FORMAT 지원 → 호환
#
# SplitNCigarReads: 생략 (PHMG_IT Phase 1 검증: -54% coverage 손실)
#
# Usage: ./mutect2_per_sample.sh <SAMPLE> [THREADS]
#   ex)  ./mutect2_per_sample.sh C1 16

SAMPLE=$1
THREADS=${2:-16}

WORK="/home/yusanghyeon/RAT_project/PHMG_IT"
REF="${WORK}/reference/rn7.fa"
RESULT_BASE="${WORK}/results/mutect2"
TMPDIR="${RESULT_BASE}/tmp/${SAMPLE}"

mkdir -p "${TMPDIR}" \
         "${RESULT_BASE}/markdup" \
         "${RESULT_BASE}/raw" \
         "${RESULT_BASE}/filtered" \
         "${RESULT_BASE}/pass"

GATK="${HOME}/RAT_project/miniforge3/envs/rnaseq/bin/gatk"
SAMTOOLS="${HOME}/RAT_project/miniforge3/envs/rnaseq/bin/samtools"
BCFTOOLS="${HOME}/RAT_project/miniforge3/envs/rnaseq/bin/bcftools"

# ── BAM 경로 결정 ────────────────────────────────────────────────────────────
if [ "${SAMPLE}" = "C1" ]; then
    INPUT_BAM="${WORK}/results/pilot/arm_C_star_aggressive/C1_Aligned.sortedByCoord.out.bam"
else
    INPUT_BAM="${WORK}/results/pilot/arm_C_${SAMPLE}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
fi

DEDUP_BAM="${RESULT_BASE}/markdup/${SAMPLE}.dedup.bam"
RAW_VCF="${RESULT_BASE}/raw/${SAMPLE}.mutect2.vcf.gz"
FILTERED_VCF="${RESULT_BASE}/filtered/${SAMPLE}.filtered.vcf.gz"
PASS_VCF="${RESULT_BASE}/pass/${SAMPLE}.pass.vcf.gz"

if [ ! -f "${INPUT_BAM}" ]; then
    echo "ERROR: BAM not found: ${INPUT_BAM}"
    exit 1
fi

echo "============================================"
echo "[$(date '+%H:%M:%S')] MuTect2: ${SAMPLE}"
echo "============================================"

# ── Step 1: MarkDuplicates ───────────────────────────────────────────────────
if [ ! -f "${DEDUP_BAM}" ]; then
    echo "[$(date '+%H:%M:%S')] Step 1: MarkDuplicates"
    ${GATK} MarkDuplicates \
        -I "${INPUT_BAM}" \
        -O "${DEDUP_BAM}" \
        -M "${TMPDIR}/${SAMPLE}.dup_metrics.txt" \
        --TMP_DIR "${TMPDIR}" \
        --VALIDATION_STRINGENCY SILENT
    ${SAMTOOLS} index "${DEDUP_BAM}"
    echo "[$(date '+%H:%M:%S')] Step 1 done"
else
    echo "[$(date '+%H:%M:%S')] Step 1: MarkDuplicates already done, skipping"
fi

# ── Step 2: MuTect2 (tumor-only) ────────────────────────────────────────────
# tumor-only mode: PON 없이 실행
# 이후 coworker analyzer가 C1-C5 VCF를 reference로 differential filtering 수행
#
# RNA-seq 관련 설정:
#   --dont-use-soft-clipped-bases false (default): soft-clipped base 허용 (Arm C 일관성)
#   --callable-depth 1: 최소 depth 1 이상 영역에서 calling (coverage 최대화)
if [ ! -f "${RAW_VCF}" ]; then
    echo "[$(date '+%H:%M:%S')] Step 2: MuTect2"
    ${GATK} Mutect2 \
        -R "${REF}" \
        -I "${DEDUP_BAM}" \
        -tumor "${SAMPLE}" \
        --callable-depth 1 \
        --native-pair-hmm-threads "${THREADS}" \
        --tmp-dir "${TMPDIR}" \
        -O "${RAW_VCF}"
    echo "[$(date '+%H:%M:%S')] Step 2 done"
else
    echo "[$(date '+%H:%M:%S')] Step 2: MuTect2 already done, skipping"
fi

# ── Step 3: FilterMutectCalls ────────────────────────────────────────────────
if [ ! -f "${FILTERED_VCF}" ]; then
    echo "[$(date '+%H:%M:%S')] Step 3: FilterMutectCalls"
    ${GATK} FilterMutectCalls \
        -R "${REF}" \
        -V "${RAW_VCF}" \
        --tmp-dir "${TMPDIR}" \
        -O "${FILTERED_VCF}"
    echo "[$(date '+%H:%M:%S')] Step 3 done"
else
    echo "[$(date '+%H:%M:%S')] Step 3: FilterMutectCalls already done, skipping"
fi

# ── Step 4: PASS 필터링 ──────────────────────────────────────────────────────
echo "[$(date '+%H:%M:%S')] Step 4: PASS filter"
${BCFTOOLS} view -f PASS "${FILTERED_VCF}" -O z -o "${PASS_VCF}"
${BCFTOOLS} index "${PASS_VCF}"

PASS_COUNT=$(${BCFTOOLS} stats "${PASS_VCF}" | grep "^SN" | grep "number of SNPs" | awk '{print $NF}')
echo "[$(date '+%H:%M:%S')] Step 4 done — PASS SNPs: ${PASS_COUNT}"

# ── Cleanup ──────────────────────────────────────────────────────────────────
rm -rf "${TMPDIR}"

echo "============================================"
echo "[$(date '+%H:%M:%S')] ${SAMPLE} COMPLETED"
echo "  dedup BAM : ${DEDUP_BAM}"
echo "  PASS VCF  : ${PASS_VCF}"
echo "============================================"
