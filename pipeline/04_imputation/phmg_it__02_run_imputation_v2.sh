#!/bin/bash
set -euo pipefail

# ============================================================
# RNA-seq Genotype Imputation Pipeline v2 (optimized)
# ============================================================
# 변경점: per-chromosome 병렬 처리, -T (targets) 사용으로 속도 대폭 개선
# ============================================================

eval "$($HOME/RAT_project/miniforge3/bin/conda shell.bash hook 2>/dev/null)"
conda activate rnaseq

WORK="/home/yusanghyeon/RAT_project/PHMG_IT/imputation"
REF="/home/yusanghyeon/RAT_project/PHMG_IT/reference/rn7.fa"
BEAGLE="${WORK}/tools/beagle.jar"
BAM_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/mutect2/markdup"
GT_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/ground_truth"
PANEL_DIR="${WORK}/ref_panel/merged"
RNA_GL_DIR="${WORK}/rna_gl"
IMPUTED_DIR="${WORK}/imputed"
EVAL_DIR="${WORK}/evaluation"

mkdir -p "${RNA_GL_DIR}" "${IMPUTED_DIR}" "${EVAL_DIR}"

SAMPLES=(C1 C2 C3 C4 C5 P1 P2 P3 P4 P5 P6 P7 P8 P9 P10)
CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10
        chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20
        chrX)

N_PARALLEL=8

echo "============================================"
echo "RNA-seq Genotype Imputation Pipeline v2"
echo "Started: $(date)"
echo "============================================"

# ─── Step 2: GL 추출 (per-chr, 병렬) ────────────────────────────────────────
echo ""
echo "[Step 2] RNA-seq BAM에서 genotype likelihood 추출 (per-chr parallel)"

PIDS=()

for SAMPLE in "${SAMPLES[@]}"; do
    for CHR in "${CHROMS[@]}"; do
        GL_CHR="${RNA_GL_DIR}/${SAMPLE}_${CHR}_gl.vcf.gz"
        [ -f "${GL_CHR}" ] && continue

        PANEL="${PANEL_DIR}/hrdp_${CHR}.vcf.gz"
        [ ! -f "${PANEL}" ] && continue
        BAM="${BAM_DIR}/${SAMPLE}.dedup.bam"

        # 슬롯 대기
        while [ ${#PIDS[@]} -ge ${N_PARALLEL} ]; do
            NEW_PIDS=()
            for pid in "${PIDS[@]}"; do
                if kill -0 "${pid}" 2>/dev/null; then
                    NEW_PIDS+=("${pid}")
                fi
            done
            PIDS=("${NEW_PIDS[@]+"${NEW_PIDS[@]}"}")
            [ ${#PIDS[@]} -ge ${N_PARALLEL} ] && sleep 2
        done

        (
            bcftools mpileup \
                -f "${REF}" \
                -T "${PANEL}" \
                -I \
                -a "FORMAT/AD,FORMAT/DP" \
                --min-MQ 0 \
                "${BAM}" 2>/dev/null | \
            bcftools call -m -Oz -o "${GL_CHR}" 2>/dev/null && \
            bcftools index -t "${GL_CHR}" 2>/dev/null
        ) &
        PIDS+=("$!")
    done
    echo "  ${SAMPLE} — all chr jobs submitted"
done

wait
echo "[Step 2] 완료 — $(date)"

# ─── Step 2b: 샘플별 chr merge ──────────────────────────────────────────────
echo ""
echo "[Step 2b] 샘플별 chromosome GL merge"

for SAMPLE in "${SAMPLES[@]}"; do
    MERGED_GL="${RNA_GL_DIR}/${SAMPLE}_rna_gl.vcf.gz"
    [ -f "${MERGED_GL}" ] && continue

    CHR_FILES=()
    for CHR in "${CHROMS[@]}"; do
        CHR_FILE="${RNA_GL_DIR}/${SAMPLE}_${CHR}_gl.vcf.gz"
        [ -f "${CHR_FILE}" ] && CHR_FILES+=("${CHR_FILE}")
    done

    if [ ${#CHR_FILES[@]} -gt 0 ]; then
        bcftools concat "${CHR_FILES[@]}" -Oz -o "${MERGED_GL}" 2>/dev/null
        bcftools index -t "${MERGED_GL}" 2>/dev/null
        N=$(bcftools view -H "${MERGED_GL}" 2>/dev/null | wc -l)
        echo "  ${SAMPLE}: ${N} sites"
    fi
done
echo "[Step 2b] 완료 — $(date)"

# ─── Step 3: Beagle imputation (per-chr, 병렬) ──────────────────────────────
echo ""
echo "[Step 3] Beagle phasing + imputation"

PIDS=()

for CHR in "${CHROMS[@]}"; do
    PANEL="${PANEL_DIR}/hrdp_${CHR}.vcf.gz"
    [ ! -f "${PANEL}" ] && continue

    for SAMPLE in "${SAMPLES[@]}"; do
        IMPUTED_OUT="${IMPUTED_DIR}/${SAMPLE}_${CHR}"
        [ -f "${IMPUTED_OUT}.vcf.gz" ] && continue

        GL_CHR="${RNA_GL_DIR}/${SAMPLE}_${CHR}_gl.vcf.gz"
        [ ! -f "${GL_CHR}" ] && continue

        # 슬롯 대기
        while [ ${#PIDS[@]} -ge ${N_PARALLEL} ]; do
            NEW_PIDS=()
            for pid in "${PIDS[@]}"; do
                if kill -0 "${pid}" 2>/dev/null; then
                    NEW_PIDS+=("${pid}")
                fi
            done
            PIDS=("${NEW_PIDS[@]+"${NEW_PIDS[@]}"}")
            [ ${#PIDS[@]} -ge ${N_PARALLEL} ] && sleep 2
        done

        (
            java -Xmx2g -jar "${BEAGLE}" \
                ref="${PANEL}" \
                gt="${GL_CHR}" \
                chrom="${CHR}" \
                out="${IMPUTED_OUT}" \
                nthreads=2 \
                2>/dev/null
        ) &
        PIDS+=("$!")
    done
done
echo "  모든 imputation job 제출 완료, 대기 중..."
wait
echo "[Step 3] 완료 — $(date)"

# ─── Step 4: 샘플별 imputed chromosome merge ────────────────────────────────
echo ""
echo "[Step 4] 샘플별 imputed VCF merge"

for SAMPLE in "${SAMPLES[@]}"; do
    FINAL="${IMPUTED_DIR}/${SAMPLE}_imputed.vcf.gz"
    [ -f "${FINAL}" ] && continue

    CHR_FILES=()
    for CHR in "${CHROMS[@]}"; do
        CHR_VCF="${IMPUTED_DIR}/${SAMPLE}_${CHR}.vcf.gz"
        [ -f "${CHR_VCF}" ] && CHR_FILES+=("${CHR_VCF}")
    done

    if [ ${#CHR_FILES[@]} -gt 0 ]; then
        bcftools concat "${CHR_FILES[@]}" -Oz -o "${FINAL}" 2>/dev/null
        bcftools index -t "${FINAL}" 2>/dev/null
        N=$(bcftools view -H "${FINAL}" 2>/dev/null | wc -l)
        echo "  ${SAMPLE}: ${N} imputed variants"
    fi
done
echo "[Step 4] 완료 — $(date)"

# ─── Step 5: Evaluate ───────────────────────────────────────────────────────
echo ""
echo "[Step 5] DNA WGS 대비 정확도 평가"
echo ""
echo "Sample | DNA_SNPs | Before | After | Cov_Before(%) | Cov_After(%)"
echo "-------|----------|--------|-------|---------------|-------------"

RESULT="${EVAL_DIR}/imputation_results.tsv"
echo -e "Sample\tDNA_SNPs\tBefore\tAfter\tCov_Before\tCov_After" > "${RESULT}"

for SAMPLE in "${SAMPLES[@]}"; do
    DNA_VCF="${GT_DIR}/${SAMPLE}_dna_snps.vcf.gz"
    RNA_BEFORE="/home/yusanghyeon/RAT_project/PHMG_IT/results/mutect2/force_pass/${SAMPLE}.force_pass.vcf.gz"
    RNA_AFTER="${IMPUTED_DIR}/${SAMPLE}_imputed.vcf.gz"

    [ ! -f "${RNA_AFTER}" ] && continue

    DNA_TOTAL=$(bcftools view -H "${DNA_VCF}" 2>/dev/null | wc -l)

    # Before: RNA PASS vs DNA
    TMPDIR_B=$(mktemp -d)
    bcftools isec -p "${TMPDIR_B}" "${DNA_VCF}" "${RNA_BEFORE}" 2>/dev/null
    BEFORE=$(grep -v "^#" "${TMPDIR_B}/0002.vcf" 2>/dev/null | wc -l)
    rm -rf "${TMPDIR_B}"

    # After: imputed (ALT only) vs DNA
    TMPDIR_A=$(mktemp -d)
    bcftools view -i 'GT="alt"' "${RNA_AFTER}" -Oz -o "${TMPDIR_A}/imp_alt.vcf.gz" 2>/dev/null
    bcftools index -t "${TMPDIR_A}/imp_alt.vcf.gz" 2>/dev/null
    bcftools isec -p "${TMPDIR_A}/isec" "${DNA_VCF}" "${TMPDIR_A}/imp_alt.vcf.gz" 2>/dev/null
    AFTER=$(grep -v "^#" "${TMPDIR_A}/isec/0002.vcf" 2>/dev/null | wc -l)
    rm -rf "${TMPDIR_A}"

    PCT_B=$(awk "BEGIN {printf \"%.1f\", ${BEFORE}/${DNA_TOTAL}*100}")
    PCT_A=$(awk "BEGIN {printf \"%.1f\", ${AFTER}/${DNA_TOTAL}*100}")

    echo "${SAMPLE} | ${DNA_TOTAL} | ${BEFORE} | ${AFTER} | ${PCT_B}% | ${PCT_A}%"
    echo -e "${SAMPLE}\t${DNA_TOTAL}\t${BEFORE}\t${AFTER}\t${PCT_B}\t${PCT_A}" >> "${RESULT}"
done

echo ""
echo "============================================"
echo "전체 완료: $(date)"
echo "결과: ${RESULT}"
echo "Disk: $(df -h /home/yusanghyeon | tail -1 | awk '{print $4}') free"
echo "============================================"
