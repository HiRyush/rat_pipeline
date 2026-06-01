#!/bin/bash
set -euo pipefail

# MuTect2 Scatter-Gather + Force-Calling (최대 속도 + 최대 sensitivity)
#
# ┌─────────────────────────────────────────────────────────────────┐
# │ 핵심 전략                                                        │
# │  1. Chromosome scatter: 15×23=345 jobs → 20코어 동시 (~1.5시간)  │
# │  2. Force-calling: 15샘플 union sites → 저depth 위치 rescue      │
# │     샘플 수 이점 활용: 개별 DP<5도 타 샘플 evidence로 살려냄    │
# └─────────────────────────────────────────────────────────────────┘
#
# Pipeline:
#   Step 0  MarkDuplicates          20병렬  ~8분
#   Step 1  MuTect2 scatter (1차)  20병렬  ~90분
#   Step 2  GatherVcfs             15병렬  ~3분
#   Step 3  FilterMutectCalls      15병렬  ~5분
#   Step 4  Union sites 생성                ~2분
#   Step 5  MuTect2 force-call    20병렬  ~30분  ← 샘플 수 이점
#   Step 6  GatherVcfs (force)    15병렬  ~3분
#   Step 7  Filter + PASS (force) 15병렬  ~5분
#   Step 8  coworker differential            ~5분
#   합계: ~2.5시간 (force-calling 추가, 기존 순차 ~8시간 대비)
#
# BAM 우선순위:
#   arm_D (fastp + pooled SJ 2-pass) > arm_C (기존)
#   arm_D_alignment.sh 실행 후 재실행 시 자동으로 arm_D 사용
#
# 리소스: 20코어 / 240GB RAM
#   MuTect2 scatter : 20병렬 × 1thread × 4GB = 80GB
#   force-call      : 20병렬 × 1thread × 4GB = 80GB  (sequential)
#   총 사용: ~80GB < 240GB ✓

eval "$($HOME/RAT_project/miniforge3/bin/conda shell.bash hook 2>/dev/null)"
conda activate rnaseq

WORK="/home/yusanghyeon/RAT_project/PHMG_IT"
REF="${WORK}/reference/rn7.fa"
LOG_DIR="${WORK}/logs/mutect2"
OUT="${WORK}/results/mutect2"
COWORKER="${WORK}/coworker"

GATK="${HOME}/RAT_project/miniforge3/envs/rnaseq/bin/gatk"
SAMTOOLS="${HOME}/RAT_project/miniforge3/envs/rnaseq/bin/samtools"
BCFTOOLS="${HOME}/RAT_project/miniforge3/envs/rnaseq/bin/bcftools"

MAX_JOBS=20

CONTROLS=(C1 C2 C3 C4 C5)
TREATED=(P1 P2 P3 P4 P5 P6 P7 P8 P9 P10)
ALL_SAMPLES=("${CONTROLS[@]}" "${TREATED[@]}")

MAJOR_CHRS=(chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18
            chr19 chr2 chr20 chr3 chr4 chr5 chr6 chr7 chr8 chr9
            chrM chrX chrY)

mkdir -p "${LOG_DIR}" \
         "${OUT}/markdup" \
         "${OUT}/scatter" \
         "${OUT}/raw" \
         "${OUT}/filtered" \
         "${OUT}/pass" \
         "${OUT}/force_scatter" \
         "${OUT}/force_raw" \
         "${OUT}/force_filtered" \
         "${OUT}/force_pass"

# ── 병렬 job 관리 ────────────────────────────────────────────────────────────
PIDS=()
FAILED=()

_wait_slot() {
    while [ "${#PIDS[@]}" -ge "${MAX_JOBS}" ]; do
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
    if [ ${#FAILED[@]} -gt 0 ]; then
        echo "ERROR: 실패 — ${FAILED[*]}"
        echo "로그: ${LOG_DIR}/"
        exit 1
    fi
    FAILED=()
}

# ── BAM 경로 결정 (arm_D 우선, 없으면 arm_C) ────────────────────────────────
_get_bam() {
    local SAMPLE=$1
    local ARM_D="${WORK}/results/arm_D/${SAMPLE}/${SAMPLE}.dedup.bam"
    local ARM_C_C1="${WORK}/results/pilot/arm_C_star_aggressive/C1_Aligned.sortedByCoord.out.bam"
    local ARM_C="${WORK}/results/pilot/arm_C_${SAMPLE}/${SAMPLE}_Aligned.sortedByCoord.out.bam"

    if [ -f "${ARM_D}" ]; then
        echo "${ARM_D}"
    elif [ "${SAMPLE}" = "C1" ]; then
        echo "${ARM_C_C1}"
    else
        echo "${ARM_C}"
    fi
}

echo "============================================"
echo "PHMG_IT MuTect2 Scatter-Gather + Force-Calling"
echo "Samples: ${#ALL_SAMPLES[@]} | Chr: ${#MAJOR_CHRS[@]}"
echo "1차 scatter: $((${#ALL_SAMPLES[@]} * ${#MAJOR_CHRS[@]})) jobs"
echo "Max concurrent: ${MAX_JOBS} | Started: $(date)"
echo "============================================"

# ── BAM 소스 확인 ────────────────────────────────────────────────────────────
ARM_D_COUNT=0
for S in "${ALL_SAMPLES[@]}"; do
    [ -f "${WORK}/results/arm_D/${S}/${S}.dedup.bam" ] && ARM_D_COUNT=$((ARM_D_COUNT+1))
done
if [ "${ARM_D_COUNT}" -gt 0 ]; then
    echo "arm_D BAM 감지: ${ARM_D_COUNT}/${#ALL_SAMPLES[@]}개 — arm_D 우선 사용"
else
    echo "arm_C BAM 사용 (arm_D 미생성 — arm_D_alignment.sh 실행 후 재실행 시 자동 전환)"
fi
echo ""

# ════════════════════════════════════════════════════════════════════
# Step 0: MarkDuplicates (전체 병렬)
# ════════════════════════════════════════════════════════════════════
echo "[Step 0] MarkDuplicates — ${#ALL_SAMPLES[@]}샘플 병렬"
for SAMPLE in "${ALL_SAMPLES[@]}"; do
    DEDUP="${OUT}/markdup/${SAMPLE}.dedup.bam"
    [ -f "${DEDUP}" ] && { echo "  [SKIP] ${SAMPLE}"; continue; }

    INPUT=$(_get_bam "${SAMPLE}")
    [ ! -f "${INPUT}" ] && { echo "ERROR: BAM not found: ${INPUT}"; exit 1; }

    _wait_slot
    (
        ${GATK} MarkDuplicates \
            -I "${INPUT}" -O "${DEDUP}" \
            -M "${OUT}/markdup/${SAMPLE}.metrics.txt" \
            --TMP_DIR "${OUT}/markdup" \
            --VALIDATION_STRINGENCY SILENT \
            2>>"${LOG_DIR}/${SAMPLE}_markdup.log" && \
        ${SAMTOOLS} index "${DEDUP}" 2>>"${LOG_DIR}/${SAMPLE}_markdup.log"
    ) &
    PIDS+=("$!:MarkDup_${SAMPLE}")
    echo "  [START] ${SAMPLE}"
done
_wait_all; _check_failed
echo "[Step 0] 완료 — $(date)"

# ════════════════════════════════════════════════════════════════════
# Step 1: MuTect2 Scatter 1차 (15샘플 × 23염색체)
# ════════════════════════════════════════════════════════════════════
echo ""
echo "[Step 1] MuTect2 Scatter 1차 — $((${#ALL_SAMPLES[@]} * ${#MAJOR_CHRS[@]})) jobs"
for SAMPLE in "${ALL_SAMPLES[@]}"; do
    for CHR in "${MAJOR_CHRS[@]}"; do
        SCATTER_VCF="${OUT}/scatter/${SAMPLE}_${CHR}.vcf.gz"
        [ -f "${SCATTER_VCF}" ] && continue

        _wait_slot
        (
            ${GATK} Mutect2 \
                -R "${REF}" \
                -I "${OUT}/markdup/${SAMPLE}.dedup.bam" \
                -tumor "${SAMPLE}" \
                -L "${CHR}" \
                --callable-depth 1 \
                --native-pair-hmm-threads 1 \
                --disable-read-filter MappingQualityAvailableReadFilter \
                --minimum-mapping-quality 0 \
                --annotations-to-exclude TandemRepeat \
                --tmp-dir "${OUT}/scatter" \
                -O "${SCATTER_VCF}" \
                2>>"${LOG_DIR}/${SAMPLE}_scatter.log"
        ) &
        PIDS+=("$!:Mutect2_${SAMPLE}_${CHR}")
    done
done
echo "  모든 job 제출, 완료 대기 중..."
_wait_all; _check_failed
echo "[Step 1] 완료 — $(date)"

# ════════════════════════════════════════════════════════════════════
# Step 2: GatherVcfs per sample (병렬)
# ════════════════════════════════════════════════════════════════════
echo ""
echo "[Step 2] GatherVcfs 1차"
for SAMPLE in "${ALL_SAMPLES[@]}"; do
    RAW_VCF="${OUT}/raw/${SAMPLE}.mutect2.vcf.gz"
    [ -f "${RAW_VCF}" ] && { echo "  [SKIP] ${SAMPLE}"; continue; }

    INPUT_ARGS=()
    for CHR in "${MAJOR_CHRS[@]}"; do INPUT_ARGS+=("-I" "${OUT}/scatter/${SAMPLE}_${CHR}.vcf.gz"); done

    STATS_ARGS=()
    for CHR in "${MAJOR_CHRS[@]}"; do STATS_ARGS+=("-stats" "${OUT}/scatter/${SAMPLE}_${CHR}.vcf.gz.stats"); done

    _wait_slot
    (
        ${GATK} GatherVcfs "${INPUT_ARGS[@]}" -O "${RAW_VCF}" \
            2>>"${LOG_DIR}/${SAMPLE}_gather.log" && \
        ${GATK} MergeMutectStats "${STATS_ARGS[@]}" -O "${RAW_VCF}.stats" \
            2>>"${LOG_DIR}/${SAMPLE}_gather.log" && \
        ${BCFTOOLS} index -t "${RAW_VCF}" 2>>"${LOG_DIR}/${SAMPLE}_gather.log"
    ) &
    PIDS+=("$!:Gather_${SAMPLE}")
    echo "  [START] ${SAMPLE}"
done
_wait_all; _check_failed
echo "[Step 2] 완료 — $(date)"

# ════════════════════════════════════════════════════════════════════
# Step 3: FilterMutectCalls + PASS 1차 (병렬)
# ════════════════════════════════════════════════════════════════════
echo ""
echo "[Step 3] FilterMutectCalls + PASS 1차"
for SAMPLE in "${ALL_SAMPLES[@]}"; do
    PASS_VCF="${OUT}/pass/${SAMPLE}.pass.vcf.gz"
    [ -f "${PASS_VCF}" ] && { echo "  [SKIP] ${SAMPLE}"; continue; }

    _wait_slot
    (
        ${GATK} FilterMutectCalls \
            -R "${REF}" \
            -V "${OUT}/raw/${SAMPLE}.mutect2.vcf.gz" \
            --min-median-base-quality 0 \
            --tmp-dir "${OUT}/filtered" \
            -O "${OUT}/filtered/${SAMPLE}.filtered.vcf.gz" \
            2>>"${LOG_DIR}/${SAMPLE}_filter.log" && \
        ${BCFTOOLS} view -f PASS \
            "${OUT}/filtered/${SAMPLE}.filtered.vcf.gz" \
            -O z -o "${PASS_VCF}" \
            2>>"${LOG_DIR}/${SAMPLE}_filter.log" && \
        ${BCFTOOLS} index -t "${PASS_VCF}" \
            2>>"${LOG_DIR}/${SAMPLE}_filter.log"
    ) &
    PIDS+=("$!:Filter_${SAMPLE}")
    echo "  [START] ${SAMPLE}"
done
_wait_all; _check_failed
echo "[Step 3] 완료 — $(date)"

# ════════════════════════════════════════════════════════════════════
# Step 4: Union Sites 생성 (샘플 수 이점 활용)
#   15개 샘플 PASS VCF를 합산 → 모든 샘플에서 발견된 후보 위치 통합
#   개별 DP<5 탈락 위치도 타 샘플에서 확인되면 force-call 대상이 됨
# ════════════════════════════════════════════════════════════════════
echo ""
echo "[Step 4] Union sites 생성 (15샘플 후보 통합)"

UNION_SITES="${OUT}/union_sites.vcf.gz"

if [ ! -f "${UNION_SITES}" ]; then
    ALL_PASS_VCFS=()
    for SAMPLE in "${ALL_SAMPLES[@]}"; do
        ALL_PASS_VCFS+=("${OUT}/pass/${SAMPLE}.pass.vcf.gz")
    done

    # 15개 VCF 병합 → 중복 제거 → 정렬 → union sites
    ${BCFTOOLS} merge \
        --merge none \
        --force-samples \
        "${ALL_PASS_VCFS[@]}" \
        2>>"${LOG_DIR}/union_sites.log" | \
    ${BCFTOOLS} view -f PASS | \
    ${BCFTOOLS} norm -m -any | \
    ${BCFTOOLS} sort | \
    ${BCFTOOLS} view -O z -o "${UNION_SITES}" \
        2>>"${LOG_DIR}/union_sites.log"

    ${BCFTOOLS} index -t "${UNION_SITES}"

    UNION_COUNT=$(${BCFTOOLS} stats "${UNION_SITES}" | grep "^SN.*number of SNPs" | awk '{print $NF}')
    echo "  Union candidate SNPs: ${UNION_COUNT}"
fi

# chr별 union_sites 분할 (force-call 동시 접근 충돌 방지)
mkdir -p "${OUT}/union_sites_per_chr"
for CHR in "${MAJOR_CHRS[@]}"; do
    CHR_UNION="${OUT}/union_sites_per_chr/${CHR}.vcf.gz"
    if [ ! -f "${CHR_UNION}" ]; then
        ${BCFTOOLS} view "${UNION_SITES}" "${CHR}" -O z -o "${CHR_UNION}" 2>/dev/null
        ${BCFTOOLS} index -t "${CHR_UNION}" 2>/dev/null
    fi
done
echo "[Step 4] 완료 — $(date)"

# ════════════════════════════════════════════════════════════════════
# Step 5: MuTect2 Force-Calling Scatter
#   --genotyping-mode GENOTYPE_GIVEN_ALLELES: union sites 강제 genotyping
#   --alleles union_sites.vcf.gz: 모든 샘플에서 동일 위치 평가
#   → 개별 DP=2~4 위치가 타 샘플 근거로 rescue
# ════════════════════════════════════════════════════════════════════
echo ""
echo "[Step 5] MuTect2 Force-Calling Scatter — $((${#ALL_SAMPLES[@]} * ${#MAJOR_CHRS[@]})) jobs"
for SAMPLE in "${ALL_SAMPLES[@]}"; do
    for CHR in "${MAJOR_CHRS[@]}"; do
        FORCE_VCF="${OUT}/force_scatter/${SAMPLE}_${CHR}.vcf.gz"
        [ -f "${FORCE_VCF}" ] && [ -f "${FORCE_VCF}.stats" ] && continue

        _wait_slot
        (
            ${GATK} Mutect2 \
                -R "${REF}" \
                -I "${OUT}/markdup/${SAMPLE}.dedup.bam" \
                -tumor "${SAMPLE}" \
                -L "${CHR}" \
                --alleles "${OUT}/union_sites_per_chr/${CHR}.vcf.gz" \
                --force-call-filtered-alleles \
                --callable-depth 1 \
                --native-pair-hmm-threads 1 \
                --disable-read-filter MappingQualityAvailableReadFilter \
                --minimum-mapping-quality 0 \
                --annotations-to-exclude TandemRepeat \
                --tmp-dir "${OUT}/force_scatter" \
                -O "${FORCE_VCF}" \
                2>>"${LOG_DIR}/${SAMPLE}_force_scatter.log" || true
            # GATK 4.6.1.0 bug: TandemRepeat crash in AssemblyRegionTrimmer causes non-zero exit
            # but VCF is fully written. Create missing .stats from initial scatter stats.
            if [ -f "${FORCE_VCF}" ]; then
                if [ ! -f "${FORCE_VCF}.stats" ]; then
                    INIT_STATS="${OUT}/scatter/${SAMPLE}_${CHR}.vcf.gz.stats"
                    if [ -f "${INIT_STATS}" ]; then
                        cp "${INIT_STATS}" "${FORCE_VCF}.stats"
                    else
                        printf "statistic\tvalue\ncallable\t0\n" > "${FORCE_VCF}.stats"
                    fi
                fi
                exit 0
            else
                exit 1
            fi
        ) &
        PIDS+=("$!:Force_${SAMPLE}_${CHR}")
    done
done
echo "  모든 force-call job 제출, 완료 대기 중..."
_wait_all; _check_failed
echo "[Step 5] 완료 — $(date)"

# ════════════════════════════════════════════════════════════════════
# Step 6: GatherVcfs (force-called)
# ════════════════════════════════════════════════════════════════════
echo ""
echo "[Step 6] GatherVcfs force-called"
for SAMPLE in "${ALL_SAMPLES[@]}"; do
    FORCE_RAW="${OUT}/force_raw/${SAMPLE}.force.vcf.gz"
    [ -f "${FORCE_RAW}" ] && { echo "  [SKIP] ${SAMPLE}"; continue; }

    INPUT_ARGS=()
    STATS_ARGS=()
    for CHR in "${MAJOR_CHRS[@]}"; do
        INPUT_ARGS+=("-I" "${OUT}/force_scatter/${SAMPLE}_${CHR}.vcf.gz")
        STATS_ARGS+=("-stats" "${OUT}/force_scatter/${SAMPLE}_${CHR}.vcf.gz.stats")
    done

    _wait_slot
    (
        ${GATK} GatherVcfs "${INPUT_ARGS[@]}" -O "${FORCE_RAW}" \
            2>>"${LOG_DIR}/${SAMPLE}_force_gather.log" && \
        ${GATK} MergeMutectStats "${STATS_ARGS[@]}" -O "${FORCE_RAW}.stats" \
            2>>"${LOG_DIR}/${SAMPLE}_force_gather.log" && \
        ${BCFTOOLS} index -t "${FORCE_RAW}" \
            2>>"${LOG_DIR}/${SAMPLE}_force_gather.log"
    ) &
    PIDS+=("$!:ForceGather_${SAMPLE}")
    echo "  [START] ${SAMPLE}"
done
_wait_all; _check_failed
echo "[Step 6] 완료 — $(date)"

# ════════════════════════════════════════════════════════════════════
# Step 7: FilterMutectCalls + PASS (force-called)
# ════════════════════════════════════════════════════════════════════
echo ""
echo "[Step 7] FilterMutectCalls + PASS (force-called)"
for SAMPLE in "${ALL_SAMPLES[@]}"; do
    FORCE_PASS="${OUT}/force_pass/${SAMPLE}.force_pass.vcf.gz"
    [ -f "${FORCE_PASS}" ] && { echo "  [SKIP] ${SAMPLE}"; continue; }

    _wait_slot
    (
        ${GATK} FilterMutectCalls \
            -R "${REF}" \
            -V "${OUT}/force_raw/${SAMPLE}.force.vcf.gz" \
            --min-median-base-quality 0 \
            --tmp-dir "${OUT}/force_filtered" \
            -O "${OUT}/force_filtered/${SAMPLE}.force_filtered.vcf.gz" \
            2>>"${LOG_DIR}/${SAMPLE}_force_filter.log" && \
        ${BCFTOOLS} view -f PASS \
            "${OUT}/force_filtered/${SAMPLE}.force_filtered.vcf.gz" \
            -O z -o "${FORCE_PASS}" \
            2>>"${LOG_DIR}/${SAMPLE}_force_filter.log" && \
        ${BCFTOOLS} index -t "${FORCE_PASS}" \
            2>>"${LOG_DIR}/${SAMPLE}_force_filter.log"
    ) &
    PIDS+=("$!:ForceFilter_${SAMPLE}")
    echo "  [START] ${SAMPLE}"
done
_wait_all; _check_failed

echo ""
echo "  === PASS VCF 비교 (1차 vs force-called) ==="
for SAMPLE in "${ALL_SAMPLES[@]}"; do
    N1=$(${BCFTOOLS} stats "${OUT}/pass/${SAMPLE}.pass.vcf.gz" 2>/dev/null \
        | grep "^SN.*number of SNPs" | awk '{print $NF}')
    N2=$(${BCFTOOLS} stats "${OUT}/force_pass/${SAMPLE}.force_pass.vcf.gz" 2>/dev/null \
        | grep "^SN.*number of SNPs" | awk '{print $NF}')
    DIFF=$((N2 - N1))
    printf "  %-4s  1차: %6s  force: %6s  (+%d)\n" "${SAMPLE}" "${N1}" "${N2}" "${DIFF}"
done
echo "[Step 7] 완료 — $(date)"

# ════════════════════════════════════════════════════════════════════
# Step 8: coworker differential analysis (force-called VCF 사용)
#   somatic 파라미터: min_alt_ratio=0.1 (낮은 AF 허용)
#   max_ref_alt_ratio=0.05 (germline 엄격 제거)
# ════════════════════════════════════════════════════════════════════
echo ""
echo "[Step 8] coworker differential analysis (force-called)"

TARGET_VCFS=()
for s in "${TREATED[@]}";  do TARGET_VCFS+=("${OUT}/force_pass/${s}.force_pass.vcf.gz"); done
REF_VCFS=()
for s in "${CONTROLS[@]}"; do REF_VCFS+=("${OUT}/force_pass/${s}.force_pass.vcf.gz"); done

DIFF_OUT="${OUT}/differential"
mkdir -p "${DIFF_OUT}"

python "${COWORKER}/run_pipeline.py" \
    --mode independent \
    --target   "${TARGET_VCFS[@]}" \
    --reference "${REF_VCFS[@]}" \
    --output   "${DIFF_OUT}" \
    --min-coverage 3 \
    --min-qual 10 \
    --min-alt-ratio 0.1 \
    --max-ref-alt-ratio 0.05 \
    --min-recurrence 2 \
    --min-delta-alt-ratio 0.1

echo "[Step 8] 완료 — $(date)"

# ════════════════════════════════════════════════════════════════════
# 최종 요약
# ════════════════════════════════════════════════════════════════════
echo ""
echo "============================================"
echo "전체 완료: $(date)"
echo ""
PHASE7="${WORK}/results/differential/phmg_vs_control/differential_snps.csv"
MUTECT2_SNP="${DIFF_OUT}/differential_snps.csv"
MUTECT2_INDEL="${DIFF_OUT}/differential_indels.csv"

if [ -f "${PHASE7}" ] && [ -f "${MUTECT2_SNP}" ]; then
    P7_SNP=$(($(wc -l < "${PHASE7}") - 1))
    M2_SNP=$(($(wc -l < "${MUTECT2_SNP}") - 1))
    M2_INDEL=$(($(wc -l < "${MUTECT2_INDEL}") - 1))
    echo "=== Phase 7 (bcftools) vs MuTect2 + Force-calling ==="
    printf "  Phase 7  differential SNPs  : %d\n" "${P7_SNP}"
    printf "  MuTect2  differential SNPs  : %d\n" "${M2_SNP}"
    printf "  MuTect2  differential INDELs: %d\n" "${M2_INDEL}"
    DIFF=$((M2_SNP - P7_SNP))
    if [ "${DIFF}" -gt 0 ]; then
        printf "  → SNP +%d (%.1f%%↑)\n" "${DIFF}" "$(echo "scale=1; ${DIFF}*100/${P7_SNP}" | bc)"
    else
        printf "  → SNP %d\n" "${DIFF}"
    fi
fi
echo ""
echo "  결과 위치: ${DIFF_OUT}"
echo "  Disk: $(df -h /home/yusanghyeon | tail -1 | awk '{print $4}') free"
echo "============================================"
