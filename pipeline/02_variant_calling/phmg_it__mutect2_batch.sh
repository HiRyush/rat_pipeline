#!/bin/bash
set -euo pipefail

# MuTect2 batch: 15개 샘플 병렬 실행 + coworker differential analysis
#
# 리소스: 20코어 / 240GB RAM
#   N_PARALLEL=4  샘플 동시 실행 수
#   THREADS=4     MuTect2 --native-pair-hmm-threads per sample
#   → 총 16 threads / ~32GB RAM 사용
#
# Step 1: 15개 샘플 MuTect2 병렬
# Step 2: coworker run_pipeline.py (P1-P10 target / C1-C5 reference)

eval "$($HOME/RAT_project/miniforge3/bin/conda shell.bash hook 2>/dev/null)"
conda activate rnaseq

WORK="/home/yusanghyeon/RAT_project/PHMG_IT"
SCRIPT="${WORK}/scripts/mutect2_per_sample.sh"
LOG_DIR="${WORK}/logs/mutect2"
PASS_DIR="${WORK}/results/mutect2/pass"
COWORKER="${WORK}/coworker"

N_PARALLEL=4
THREADS=4

mkdir -p "${LOG_DIR}"

CONTROLS=(C1 C2 C3 C4 C5)
TREATED=(P1 P2 P3 P4 P5 P6 P7 P8 P9 P10)
ALL_SAMPLES=("${CONTROLS[@]}" "${TREATED[@]}")
TOTAL=${#ALL_SAMPLES[@]}

echo "============================================"
echo "PHMG_IT MuTect2 Batch (병렬 실행)"
echo "Samples: ${TOTAL} | Parallel: ${N_PARALLEL} | Threads/sample: ${THREADS}"
echo "Started: $(date)"
echo "============================================"

# ── Step 1: MuTect2 병렬 실행 ────────────────────────────────────────────────
PIDS=()
RUNNING=0
FAILED_SAMPLES=()

for SAMPLE in "${ALL_SAMPLES[@]}"; do
    if [ -f "${PASS_DIR}/${SAMPLE}.pass.vcf.gz" ]; then
        echo "[SKIP] ${SAMPLE} — 이미 완료"
        continue
    fi

    # 최대 N_PARALLEL 유지: 슬롯이 찰 때까지 대기
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

    echo "[START] ${SAMPLE}"
    bash "${SCRIPT}" "${SAMPLE}" "${THREADS}" \
        > "${LOG_DIR}/${SAMPLE}.log" 2>&1 &
    PIDS+=("$!:${SAMPLE}")
    RUNNING=$((RUNNING + 1))
done

# 남은 작업 완료 대기
for pid_info in "${PIDS[@]+"${PIDS[@]}"}"; do
    pid="${pid_info%%:*}"
    sname="${pid_info##*:}"
    wait "${pid}" && echo "[DONE] ${sname}" || { echo "[FAIL] ${sname}"; FAILED_SAMPLES+=("${sname}"); }
done

echo ""
if [ ${#FAILED_SAMPLES[@]} -gt 0 ]; then
    echo "ERROR: 실패한 샘플 — ${FAILED_SAMPLES[*]}"
    echo "로그 확인: ${LOG_DIR}/"
    exit 1
fi
echo "MuTect2 전체 완료: $(date)"
echo "Disk: $(df -h /home/yusanghyeon | tail -1 | awk '{print $4}') free"

# ── Step 2: coworker differential analysis ───────────────────────────────────
echo ""
echo "============================================"
echo "coworker differential analysis"
echo "  target    : P1-P10 (PHMG 처리군)"
echo "  reference : C1-C5  (Control)"
echo "  파라미터  : somatic (min_alt_ratio=0.1, max_ref_alt_ratio=0.05)"
echo "============================================"

TARGET_VCFS=()
for s in "${TREATED[@]}"; do TARGET_VCFS+=("${PASS_DIR}/${s}.pass.vcf.gz"); done

REF_VCFS=()
for s in "${CONTROLS[@]}"; do REF_VCFS+=("${PASS_DIR}/${s}.pass.vcf.gz"); done

OUTPUT_DIR="${WORK}/results/mutect2/differential"
mkdir -p "${OUTPUT_DIR}"

python "${COWORKER}/run_pipeline.py" \
    --mode independent \
    --target  "${TARGET_VCFS[@]}" \
    --reference "${REF_VCFS[@]}" \
    --output "${OUTPUT_DIR}" \
    --min-coverage 5 \
    --min-qual 10 \
    --min-alt-ratio 0.1 \
    --max-ref-alt-ratio 0.05 \
    --min-recurrence 2 \
    --min-delta-alt-ratio 0.1

echo ""
echo "============================================"
echo "완료: $(date)"
echo "결과 위치: ${OUTPUT_DIR}"
echo ""
# Phase 7 (bcftools) 결과와 비교
PHASE7="${WORK}/results/differential/phmg_vs_control/differential_snps.csv"
MUTECT2_SNP="${OUTPUT_DIR}/differential_snps.csv"
if [ -f "${PHASE7}" ] && [ -f "${MUTECT2_SNP}" ]; then
    P7=$(($(wc -l < "${PHASE7}") - 1))
    M2=$(($(wc -l < "${MUTECT2_SNP}") - 1))
    echo "=== Phase 7 (bcftools) vs MuTect2 ==="
    printf "  Phase 7 differential SNPs : %d\n" "${P7}"
    printf "  MuTect2 differential SNPs : %d\n" "${M2}"
    if [ "${M2}" -gt "${P7}" ]; then
        DIFF=$((M2 - P7))
        printf "  → MuTect2 +%d SNPs (%.1f%%↑)\n" "${DIFF}" "$(echo "scale=1; ${DIFF}*100/${P7}" | bc)"
    fi
fi
echo "============================================"
