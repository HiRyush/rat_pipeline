#!/bin/bash
set -euo pipefail

# Arm D: fastp 전처리 + Pooled SJ 2-pass STAR alignment
#
# 기존 Arm C 대비 추가 최적화:
#   1. fastp: adapter trimming + quality filter → 버려지던 reads 살림 (+1~3%p)
#   2. Pooled SJ 2-pass: 15샘플 SJ 통합 → junction 인근 coverage 향상 (+3~7%p)
#
# Arm C에서 Arm D로 개선 기대:
#   DP≥5 coverage: ~43-46% → ~48-55% (추정)
#
# 출력: results/arm_D/{sample}/{sample}.dedup.bam
#   → mutect2_scatter.sh에서 arm_D BAM 자동 감지 후 사용
#
# 리소스: 20코어 / 240GB RAM
#   Pass 1 STAR: 3병렬 × 6threads (RAM: 3×50GB=150GB)
#   Pass 2 STAR: 3병렬 × 6threads (pooled SJ 반영)
#   fastp:       20병렬 × 1thread

eval "$($HOME/RAT_project/miniforge3/bin/conda shell.bash hook 2>/dev/null)"
conda activate rnaseq

WORK="/home/yusanghyeon/RAT_project/PHMG_IT"
FASTQ_DIR="/media/yusanghyeon/30B4E366B4E32D52/01_Projects/korea_PHMG/RNA_fastq"
REF="${WORK}/reference/rn7.fa"
STAR_INDEX="${WORK}/reference/star_index"
LOG_DIR="${WORK}/logs/arm_D"
OUT="${WORK}/results/arm_D"

GATK="${HOME}/RAT_project/miniforge3/envs/rnaseq/bin/gatk"
SAMTOOLS="${HOME}/RAT_project/miniforge3/envs/rnaseq/bin/samtools"
FASTP="${HOME}/RAT_project/miniforge3/envs/rnaseq/bin/fastp"
STAR="${HOME}/RAT_project/miniforge3/envs/rnaseq/bin/STAR"

N_ALIGN=3     # STAR 병렬 수 (RAM: 3×50GB=150GB)
THREADS=6     # STAR threads per sample

CONTROLS=(C1 C2 C3 C4 C5)
TREATED=(P1 P2 P3 P4 P5 P6 P7 P8 P9 P10)
ALL_SAMPLES=("${CONTROLS[@]}" "${TREATED[@]}")

mkdir -p "${LOG_DIR}" "${OUT}/trimmed" "${OUT}/sj_pass1" "${OUT}/pooled_sj"

# ── 병렬 job 관리 ────────────────────────────────────────────────────────────
PIDS=(); FAILED=()

_wait_slot_n() {
    local N=$1
    while [ "${#PIDS[@]}" -ge "${N}" ]; do
        for i in "${!PIDS[@]}"; do
            pid="${PIDS[$i]%%:*}"; tag="${PIDS[$i]##*:}"
            if ! kill -0 "${pid}" 2>/dev/null; then
                wait "${pid}" || { echo "[FAIL] ${tag}"; FAILED+=("${tag}"); }
                unset 'PIDS[i]'
                PIDS=("${PIDS[@]+"${PIDS[@]}"}")
                return
            fi
        done
        sleep 3
    done
}

_wait_all() {
    for entry in "${PIDS[@]+"${PIDS[@]}"}"; do
        pid="${entry%%:*}"; tag="${entry##*:}"
        wait "${pid}" || { echo "[FAIL] ${tag}"; FAILED+=("${tag}"); }
    done
    PIDS=()
}

_check_failed() {
    [ ${#FAILED[@]} -gt 0 ] && { echo "ERROR: ${FAILED[*]}"; exit 1; }
    FAILED=()
}

echo "============================================"
echo "Arm D: fastp + Pooled SJ 2-pass STAR"
echo "Samples: ${#ALL_SAMPLES[@]} | STAR parallel: ${N_ALIGN}"
echo "Started: $(date)"
echo "============================================"

# 외장 HDD 마운트 확인 (fastp trimmed 파일이 없을 때만 필요)
NEED_HDD=0
for S in "${ALL_SAMPLES[@]}"; do
    [ ! -f "${OUT}/trimmed/${S}_1.fastq.gz" ] && NEED_HDD=1 && break
done
[ "${NEED_HDD}" -eq 1 ] && [ ! -d "${FASTQ_DIR}" ] && { echo "ERROR: 외장 HDD 미마운트: ${FASTQ_DIR}"; exit 1; }
[ ! -d "${FASTQ_DIR}" ] && echo "외장 HDD 미마운트 — trimmed 파일 이미 존재, 계속 진행"

# ════════════════════════════════════════════════════════════════════
# Step 0: fastp 전처리 (20병렬)
#   - adapter auto-detection + trimming
#   - quality filter: Q20 이하 base trim
#   - 너무 짧은 read 제거 (< 36bp)
# ════════════════════════════════════════════════════════════════════
echo ""
echo "[Step 0] fastp 전처리 — 20병렬"
for SAMPLE in "${ALL_SAMPLES[@]}"; do
    R1_TRIM="${OUT}/trimmed/${SAMPLE}_1.fastq.gz"
    [ -f "${R1_TRIM}" ] && { echo "  [SKIP] ${SAMPLE}"; continue; }

    R1="${FASTQ_DIR}/${SAMPLE}_1.fastq.gz"
    R2="${FASTQ_DIR}/${SAMPLE}_2.fastq.gz"
    [ ! -f "${R1}" ] && { echo "ERROR: FASTQ not found: ${R1}"; exit 1; }

    _wait_slot_n 20
    (
        ${FASTP} \
            -i "${R1}" -I "${R2}" \
            -o "${OUT}/trimmed/${SAMPLE}_1.fastq.gz" \
            -O "${OUT}/trimmed/${SAMPLE}_2.fastq.gz" \
            --detect_adapter_for_pe \
            --qualified_quality_phred 20 \
            --length_required 36 \
            --thread 2 \
            --json "${OUT}/trimmed/${SAMPLE}_fastp.json" \
            --html "${OUT}/trimmed/${SAMPLE}_fastp.html" \
            2>>"${LOG_DIR}/${SAMPLE}_fastp.log"
    ) &
    PIDS+=("$!:fastp_${SAMPLE}")
    echo "  [START] fastp ${SAMPLE}"
done
_wait_all; _check_failed
echo "[Step 0] 완료 — $(date)"

# ════════════════════════════════════════════════════════════════════
# Step 1: STAR 1-pass (${N_ALIGN}병렬)
#   → SJ.out.tab 수집 목적
#   Arm C 동일 파라미터 (aggressive)
# ════════════════════════════════════════════════════════════════════
echo ""
echo "[Step 1] STAR 1-pass (SJ 수집) — ${N_ALIGN}병렬"
for SAMPLE in "${ALL_SAMPLES[@]}"; do
    SJ_FILE="${OUT}/sj_pass1/${SAMPLE}_SJ.out.tab"
    [ -f "${SJ_FILE}" ] && { echo "  [SKIP] ${SAMPLE}"; continue; }

    mkdir -p "${OUT}/sj_pass1/${SAMPLE}_tmp"
    _wait_slot_n "${N_ALIGN}"
    (
        ${STAR} \
            --genomeDir "${STAR_INDEX}" \
            --readFilesIn "${OUT}/trimmed/${SAMPLE}_1.fastq.gz" \
                          "${OUT}/trimmed/${SAMPLE}_2.fastq.gz" \
            --readFilesCommand zcat \
            --outFileNamePrefix "${OUT}/sj_pass1/${SAMPLE}_tmp/" \
            --outSAMtype BAM SortedByCoordinate \
            --runThreadN "${THREADS}" \
            --outBAMsortingThreadN 2 \
            --limitBAMsortRAM 50000000000 \
            --outFilterMultimapNmax 50 \
            --outFilterScoreMinOverLread 0.33 \
            --outFilterMatchNminOverLread 0.33 \
            --alignEndsType Local \
            --winAnchorMultimapNmax 100 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --outSAMmode None \
            2>>"${LOG_DIR}/${SAMPLE}_star_pass1.log" && \
        cp "${OUT}/sj_pass1/${SAMPLE}_tmp/SJ.out.tab" "${SJ_FILE}" && \
        rm -rf "${OUT}/sj_pass1/${SAMPLE}_tmp"
    ) &
    PIDS+=("$!:STAR1_${SAMPLE}")
    echo "  [START] STAR 1-pass ${SAMPLE}"
done
_wait_all; _check_failed
echo "[Step 1] 완료 — $(date)"

# ════════════════════════════════════════════════════════════════════
# Step 2: Pooled SJ 생성
#   15개 샘플 SJ.out.tab 통합 → non-canonical junction 제거
#   → 모든 샘플의 junction 정보 공유
# ════════════════════════════════════════════════════════════════════
echo ""
echo "[Step 2] Pooled SJ 생성"
POOLED_SJ="${OUT}/pooled_sj/pooled_SJ.out.tab"

if [ ! -f "${POOLED_SJ}" ]; then
    cat "${OUT}/sj_pass1/"*"_SJ.out.tab" | \
        awk '($5 > 0 && $7 > 2 && $6 == 0)' | \
        cut -f1-6 | sort -u > "${POOLED_SJ}"
    SJ_COUNT=$(wc -l < "${POOLED_SJ}")
    echo "  Pooled junctions: ${SJ_COUNT}"
fi
echo "[Step 2] 완료 — $(date)"

# ════════════════════════════════════════════════════════════════════
# Step 3: STAR 2-pass with Pooled SJ (${N_ALIGN}병렬)
#   → arm_D final BAM 생성
# ════════════════════════════════════════════════════════════════════
echo ""
echo "[Step 3] STAR 2-pass (Pooled SJ) — ${N_ALIGN}병렬"
for SAMPLE in "${ALL_SAMPLES[@]}"; do
    FINAL_BAM="${OUT}/${SAMPLE}/${SAMPLE}.bam"
    [ -f "${FINAL_BAM}" ] && { echo "  [SKIP] ${SAMPLE}"; continue; }

    mkdir -p "${OUT}/${SAMPLE}"
    _wait_slot_n "${N_ALIGN}"
    (
        ${STAR} \
            --genomeDir "${STAR_INDEX}" \
            --readFilesIn "${OUT}/trimmed/${SAMPLE}_1.fastq.gz" \
                          "${OUT}/trimmed/${SAMPLE}_2.fastq.gz" \
            --readFilesCommand zcat \
            --outFileNamePrefix "${OUT}/${SAMPLE}/" \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattrRGline "ID:${SAMPLE}" "SM:${SAMPLE}" "PL:ILLUMINA" "LB:lib1" \
            --runThreadN "${THREADS}" \
            --outBAMsortingThreadN 2 \
            --limitBAMsortRAM 50000000000 \
            --sjdbFileChrStartEnd "${POOLED_SJ}" \
            --outFilterMultimapNmax 50 \
            --outFilterScoreMinOverLread 0.33 \
            --outFilterMatchNminOverLread 0.33 \
            --alignEndsType Local \
            --winAnchorMultimapNmax 100 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            2>>"${LOG_DIR}/${SAMPLE}_star_pass2.log" && \
        mv "${OUT}/${SAMPLE}/Aligned.sortedByCoord.out.bam" "${FINAL_BAM}" && \
        ${SAMTOOLS} index "${FINAL_BAM}" \
            2>>"${LOG_DIR}/${SAMPLE}_star_pass2.log"
    ) &
    PIDS+=("$!:STAR2_${SAMPLE}")
    echo "  [START] STAR 2-pass ${SAMPLE}"
done
_wait_all; _check_failed
echo "[Step 3] 완료 — $(date)"

# ════════════════════════════════════════════════════════════════════
# Step 4: MarkDuplicates (20병렬)
# ════════════════════════════════════════════════════════════════════
echo ""
echo "[Step 4] MarkDuplicates — 20병렬"
for SAMPLE in "${ALL_SAMPLES[@]}"; do
    DEDUP_BAM="${OUT}/${SAMPLE}/${SAMPLE}.dedup.bam"
    [ -f "${DEDUP_BAM}" ] && { echo "  [SKIP] ${SAMPLE}"; continue; }

    _wait_slot_n 20
    (
        ${GATK} MarkDuplicates \
            -I "${OUT}/${SAMPLE}/${SAMPLE}.bam" \
            -O "${DEDUP_BAM}" \
            -M "${OUT}/${SAMPLE}/${SAMPLE}.dup_metrics.txt" \
            --TMP_DIR "${OUT}/${SAMPLE}" \
            --VALIDATION_STRINGENCY SILENT \
            2>>"${LOG_DIR}/${SAMPLE}_markdup.log" && \
        ${SAMTOOLS} index "${DEDUP_BAM}" \
            2>>"${LOG_DIR}/${SAMPLE}_markdup.log" && \
        rm -f "${OUT}/${SAMPLE}/${SAMPLE}.bam" \
              "${OUT}/${SAMPLE}/${SAMPLE}.bam.bai"
    ) &
    PIDS+=("$!:MarkDup_${SAMPLE}")
    echo "  [START] MarkDup ${SAMPLE}"
done
_wait_all; _check_failed
echo "[Step 4] 완료 — $(date)"

# ════════════════════════════════════════════════════════════════════
# 완료 요약
# ════════════════════════════════════════════════════════════════════
echo ""
echo "============================================"
echo "Arm D alignment 완료: $(date)"
echo ""
echo "arm_D BAM 위치: ${OUT}/{sample}/{sample}.dedup.bam"
echo ""
echo "다음 단계:"
echo "  mutect2_scatter.sh 재실행 시 arm_D BAM 자동 감지"
echo "  Arm C vs Arm D coverage 비교:"
echo "  bash ${WORK}/scripts/evaluate_coverage.sh"
echo "Disk: $(df -h /home/yusanghyeon | tail -1 | awk '{print $4}') free"
echo "============================================"
