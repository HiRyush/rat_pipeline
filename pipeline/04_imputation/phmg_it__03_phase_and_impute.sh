#!/bin/bash
set -euo pipefail

export PATH="/home/yusanghyeon/RAT_project/miniforge3/envs/rnaseq/bin:$PATH"

WORK="/home/yusanghyeon/RAT_project/PHMG_IT/imputation"
REF="/home/yusanghyeon/RAT_project/PHMG_IT/reference/rn7.fa"
BEAGLE="${WORK}/tools/beagle.jar"
PANEL_DIR="${WORK}/ref_panel"
RNA_GL_DIR="${WORK}/rna_gl"
IMPUTED_DIR="${WORK}/imputed"
EVAL_DIR="${WORK}/evaluation"
GT_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/ground_truth"

SAMPLES=(C1 C2 C3 C4 C5 P1 P2 P3 P4 P5 P6 P7 P8 P9 P10)
CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10
        chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20
        chrX)

mkdir -p "${PANEL_DIR}/phased" "${IMPUTED_DIR}" "${EVAL_DIR}"

echo "============================================"
echo "Phase + Impute Pipeline"
echo "Started: $(date)"
echo "============================================"

# ─── Step A: Fill missing + Phase reference panel ────────────────────────────
echo ""
echo "[Step A] Reference panel: missing→0/0 + Beagle phasing"

for CHR in "${CHROMS[@]}"; do
    PHASED="${PANEL_DIR}/phased/hrdp_${CHR}_phased.vcf.gz"
    [ -f "${PHASED}" ] && continue

    MERGED="${PANEL_DIR}/merged/hrdp_${CHR}.vcf.gz"
    FILLED="${PANEL_DIR}/phased/hrdp_${CHR}_filled.vcf.gz"

    # missing → 0/0
    if [ ! -f "${FILLED}" ]; then
        bcftools view "${MERGED}" | \
            awk 'BEGIN{OFS="\t"} /^#/{print;next} {for(i=10;i<=NF;i++){if($i ~ /^\.\/.:/){sub(/^\.\/\./,"0/0",$i)}else if($i=="./.")$i="0/0"};print}' | \
            bcftools view -Oz -o "${FILLED}"
        bcftools index -t "${FILLED}"
    fi

    # Phase
    echo "  Phasing ${CHR}..."
    java -Xmx4g -jar "${BEAGLE}" \
        gt="${FILLED}" \
        chrom="${CHR}" \
        out="${PANEL_DIR}/phased/hrdp_${CHR}_phased" \
        nthreads=4 2>/dev/null

    # filled 삭제 (디스크 절약)
    rm -f "${FILLED}" "${FILLED}.tbi"
done
echo "[Step A] 완료 — $(date)"

# ─── Step B: Imputation (per-chr, 병렬) ──────────────────────────────────────
echo ""
echo "[Step B] Beagle imputation (per-chr parallel)"

N_PARALLEL=8
PIDS=()

for CHR in "${CHROMS[@]}"; do
    PHASED="${PANEL_DIR}/phased/hrdp_${CHR}_phased.vcf.gz"
    [ ! -f "${PHASED}" ] && continue

    for SAMPLE in "${SAMPLES[@]}"; do
        IMPUTED_OUT="${IMPUTED_DIR}/${SAMPLE}_${CHR}"
        [ -f "${IMPUTED_OUT}.vcf.gz" ] && continue

        GL_CHR="${RNA_GL_DIR}/${SAMPLE}_${CHR}_gl.vcf.gz"
        [ ! -f "${GL_CHR}" ] && continue

        # 슬롯 대기
        while [ ${#PIDS[@]} -ge ${N_PARALLEL} ]; do
            NEW_PIDS=()
            for pid in "${PIDS[@]}"; do
                kill -0 "${pid}" 2>/dev/null && NEW_PIDS+=("${pid}")
            done
            PIDS=("${NEW_PIDS[@]+"${NEW_PIDS[@]}"}")
            [ ${#PIDS[@]} -ge ${N_PARALLEL} ] && sleep 2
        done

        (
            java -Xmx2g -jar "${BEAGLE}" \
                ref="${PHASED}" \
                gt="${GL_CHR}" \
                chrom="${CHR}" \
                out="${IMPUTED_OUT}" \
                nthreads=2 2>/dev/null
        ) &
        PIDS+=("$!")
    done
done
echo "  모든 job 제출 완료, 대기 중..."
wait
echo "[Step B] 완료 — $(date)"

# ─── Step C: 샘플별 merge ───────────────────────────────────────────────────
echo ""
echo "[Step C] 샘플별 chromosome merge"

for SAMPLE in "${SAMPLES[@]}"; do
    FINAL="${IMPUTED_DIR}/${SAMPLE}_imputed.vcf.gz"
    [ -f "${FINAL}" ] && continue

    CHR_FILES=()
    for CHR in "${CHROMS[@]}"; do
        F="${IMPUTED_DIR}/${SAMPLE}_${CHR}.vcf.gz"
        [ -f "${F}" ] && CHR_FILES+=("${F}")
    done

    if [ ${#CHR_FILES[@]} -gt 0 ]; then
        bcftools concat "${CHR_FILES[@]}" -Oz -o "${FINAL}"
        bcftools index -t "${FINAL}"
        N=$(bcftools view -H "${FINAL}" | wc -l)
        N_ALT=$(bcftools view -H "${FINAL}" | awk '$10 !~ /^0\|0/' | wc -l)
        echo "  ${SAMPLE}: ${N} total, ${N_ALT} ALT (imputed)"
    fi
done
echo "[Step C] 완료 — $(date)"

# ─── Step D: Evaluate ───────────────────────────────────────────────────────
echo ""
echo "[Step D] DNA WGS 대비 정확도 평가"
echo ""

RESULT="${EVAL_DIR}/imputation_results.tsv"
echo -e "Sample\tDNA_SNPs\tBefore\tAfter\tCov_Before(%)\tCov_After(%)" > "${RESULT}"

printf "%-6s | %10s | %10s | %10s | %8s | %8s\n" "Sample" "DNA_SNPs" "Before" "After" "Before%" "After%"
printf "%s\n" "-------|------------|------------|------------|----------|----------"

for SAMPLE in "${SAMPLES[@]}"; do
    DNA_VCF="${GT_DIR}/${SAMPLE}_dna_snps.vcf.gz"
    RNA_BEFORE="/home/yusanghyeon/RAT_project/PHMG_IT/results/mutect2/force_pass/${SAMPLE}.force_pass.vcf.gz"
    RNA_AFTER="${IMPUTED_DIR}/${SAMPLE}_imputed.vcf.gz"

    [ ! -f "${RNA_AFTER}" ] && continue

    DNA_TOTAL=$(bcftools view -H "${DNA_VCF}" | wc -l)

    # Before: RNA PASS vs DNA
    TMPD=$(mktemp -d)
    bcftools isec -p "${TMPD}/b" "${DNA_VCF}" "${RNA_BEFORE}" 2>/dev/null
    BEFORE=$(grep -vc "^#" "${TMPD}/b/0002.vcf" 2>/dev/null || echo 0)
    rm -rf "${TMPD}/b"

    # After: imputed ALT vs DNA
    bcftools view -i 'GT!="0|0" && GT!="0/0"' "${RNA_AFTER}" -Oz -o "${TMPD}/alt.vcf.gz" 2>/dev/null
    bcftools index -t "${TMPD}/alt.vcf.gz" 2>/dev/null
    bcftools isec -p "${TMPD}/a" "${DNA_VCF}" "${TMPD}/alt.vcf.gz" 2>/dev/null
    AFTER=$(grep -vc "^#" "${TMPD}/a/0002.vcf" 2>/dev/null || echo 0)
    rm -rf "${TMPD}"

    PCT_B=$(awk "BEGIN {printf \"%.1f\", ${BEFORE}/${DNA_TOTAL}*100}")
    PCT_A=$(awk "BEGIN {printf \"%.1f\", ${AFTER}/${DNA_TOTAL}*100}")

    printf "%-6s | %10d | %10d | %10d | %7s%% | %7s%%\n" "${SAMPLE}" "${DNA_TOTAL}" "${BEFORE}" "${AFTER}" "${PCT_B}" "${PCT_A}"
    echo -e "${SAMPLE}\t${DNA_TOTAL}\t${BEFORE}\t${AFTER}\t${PCT_B}\t${PCT_A}" >> "${RESULT}"
done

echo ""
echo "============================================"
echo "전체 완료: $(date)"
echo "결과: ${RESULT}"
echo "Disk: $(df -h /home/yusanghyeon | tail -1 | awk '{print $4}') free"
echo "============================================"
