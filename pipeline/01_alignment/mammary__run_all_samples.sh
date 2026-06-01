#!/bin/bash
set -euo pipefail

# mammary_cancer: 23개 샘플 병렬 실행
#
# 리소스: 20코어 / 240GB RAM
#   N_PARALLEL=3  샘플 동시 실행 수
#   THREADS=6     STAR + GATK threads per sample
#   → 총 18 threads / ~102GB RAM 사용 (STAR 30GB × 3)
#
# PHMG_IT MuTect2와 동시 실행 시 합산 ~134GB RAM < 240GB ✓

eval "$($HOME/RAT_project/miniforge3/bin/conda shell.bash hook 2>/dev/null)"
conda activate rnaseq

WORK="/home/yusanghyeon/RAT_project/mammary_cancer"
SCRIPT="${WORK}/scripts/run_rnaseq_variant_calling.sh"
LOG_DIR="${WORK}/logs"
RESULT_DIR="${WORK}/results/rn6"

N_PARALLEL=3
THREADS=6

mkdir -p "${LOG_DIR}" "${RESULT_DIR}"

SAMPLES=(
    SRR33625148 SRR33625149 SRR33625150 SRR33625151  # M20 (4개)
    SRR33625152 SRR33625153 SRR33625154              # M12 (3개)
    SRR33625155                                       # M6  (1개)
    SRR33625156 SRR33625157 SRR33625158 SRR33625159  # M3  (8개)
    SRR33625160 SRR33625161 SRR33625162 SRR33625163
    SRR33625164 SRR33625165 SRR33625166 SRR33625167  # M1  (7개)
    SRR33625168 SRR33625169 SRR33625170
)

TOTAL=${#SAMPLES[@]}
echo "============================================"
echo "mammary_cancer RNA-seq Variant Calling (병렬 실행)"
echo "Samples: ${TOTAL} | Parallel: ${N_PARALLEL} | Threads/sample: ${THREADS}"
echo "Reference: rn6 (Rnor_6.0)"
echo "Started: $(date)"
echo "============================================"

# ── 병렬 실행 ────────────────────────────────────────────────────────────────
PIDS=()
RUNNING=0
FAILED_SAMPLES=()

for SRR in "${SAMPLES[@]}"; do
    # 완료 체크
    if [ -f "${RESULT_DIR}/${SRR}.filtered.vcf" ]; then
        echo "[SKIP] ${SRR} — 이미 완료"
        continue
    fi

    # FASTQ 존재 확인
    if [ ! -f "${WORK}/fastq/${SRR}_1.fastq.gz" ]; then
        echo "[SKIP] ${SRR} — FASTQ 없음"
        continue
    fi

    # 슬롯 대기
    while [ ${RUNNING} -ge ${N_PARALLEL} ]; do
        for pid_idx in "${!PIDS[@]}"; do
            pid_info="${PIDS[$pid_idx]}"
            pid="${pid_info%%:*}"
            sname="${pid_info##*:}"
            if ! kill -0 "${pid}" 2>/dev/null; then
                wait "${pid}" && echo "[DONE] ${sname}" || { echo "[FAIL] ${sname}"; FAILED_SAMPLES+=("${sname}"); }
                unset 'PIDS[pid_idx]'
                PIDS=("${PIDS[@]+"${PIDS[@]}"}")
                RUNNING=$((RUNNING - 1))
                break
            fi
        done
        sleep 5
    done

    echo "[START] ${SRR}"
    bash "${SCRIPT}" "${SRR}" "${THREADS}" \
        > "${LOG_DIR}/${SRR}_pipeline.log" 2>&1 &
    PIDS+=("$!:${SRR}")
    RUNNING=$((RUNNING + 1))
    echo "  실행 중: ${RUNNING}/${N_PARALLEL} | Disk: $(df -h /home/yusanghyeon | tail -1 | awk '{print $4}') free"
done

# 남은 작업 완료 대기
for pid_info in "${PIDS[@]+"${PIDS[@]}"}"; do
    pid="${pid_info%%:*}"
    sname="${pid_info##*:}"
    wait "${pid}" && echo "[DONE] ${sname}" || { echo "[FAIL] ${sname}"; FAILED_SAMPLES+=("${sname}"); }
done

echo ""
echo "============================================"
echo "완료: $(date)"
DONE=$(ls "${RESULT_DIR}"/*.filtered.vcf 2>/dev/null | wc -l)
echo "성공: ${DONE}/${TOTAL}"
if [ ${#FAILED_SAMPLES[@]} -gt 0 ]; then
    echo "실패: ${FAILED_SAMPLES[*]}"
    echo "로그 확인: ${LOG_DIR}/"
fi
echo "Disk: $(df -h /home/yusanghyeon | tail -1 | awk '{print $4}') free"
echo "============================================"
