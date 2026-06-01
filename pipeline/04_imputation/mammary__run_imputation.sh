#!/bin/bash
set -euo pipefail

export PATH="/home/yusanghyeon/RAT_project/miniforge3/envs/rnaseq/bin:$PATH"

WORK="/home/yusanghyeon/RAT_project/mammary_cancer"
BEAGLE="/home/yusanghyeon/RAT_project/PHMG_IT/imputation/tools/beagle.jar"
REF="${WORK}/reference/rn6/Rnor_6.0.fa"
GT_DIR="${WORK}/ground_truth"
RNA_DIR="${WORK}/results/rn6"
IMP_DIR="${WORK}/imputation"
PANEL_DIR="${IMP_DIR}/ref_panel"
IMPUTED_DIR="${IMP_DIR}/imputed"
EVAL_DIR="${IMP_DIR}/evaluation"

mkdir -p "${PANEL_DIR}/phased" "${IMPUTED_DIR}" "${EVAL_DIR}"

CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 X)

# ─── Verified sample mapping (SRA API, 2026-04-24) ──────────────────────────
declare -A SAMPLE_MAP=(
    [SRR33625148]=GSM8994506_20m_4_tumor.pass.vcf
    [SRR33625149]=GSM8994505_20m_3_tumor.pass.vcf
    [SRR33625150]=GSM8994504_20m_2_tumor.pass.vcf
    [SRR33625151]=GSM8994503_20m_1_tumor.pass.vcf
    [SRR33625152]=GSM8994502_12m_3_tumor.pass.vcf
    [SRR33625153]=GSM8994501_12m_2_tumor.pass.vcf
    [SRR33625154]=GSM8994500_12m_1_tumor.pass.vcf
    [SRR33625155]=GSM8994499_6m_2_tumor.pass.vcf
    [SRR33625156]=GSM8994497_3m_8_tumor.pass.vcf
    [SRR33625157]=GSM8994496_3m_7_tumor.pass.vcf
    [SRR33625158]=GSM8994495_3m_6_tumor.pass.vcf
    [SRR33625159]=GSM8994494_3m_5_tumor.pass.vcf
    [SRR33625160]=GSM8994493_3m_4_tumor.pass.vcf
    [SRR33625161]=GSM8994492_3m_3_tumor.pass.vcf
    [SRR33625162]=GSM8994491_3m_2_tumor.pass.vcf
    [SRR33625163]=GSM8994490_3m_1_tumor.pass.vcf
    [SRR33625164]=GSM8994489_1m_8_tumor.pass.vcf
    [SRR33625165]=GSM8994488_1m_7_tumor.pass.vcf
    [SRR33625166]=GSM8994487_1m_6_tumor.pass.vcf
    [SRR33625167]=GSM8994486_1m_5_tumor.pass.vcf
    [SRR33625168]=GSM8994484_1m_3_tumor.pass.vcf
    [SRR33625169]=GSM8994483_1m_2_tumor.pass.vcf
    [SRR33625170]=GSM8994482_1m_1_tumor.pass.vcf
)
SAMPLES=(${!SAMPLE_MAP[@]})

echo "============================================"
echo "Mammary Cancer Imputation Pipeline"
echo "Samples: ${#SAMPLES[@]} | Reference: rn6"
echo "Started: $(date)"
echo "============================================"

# ─── Step 1: Download HRDP rn6 48-sample joint SNP VCF ──────────────────────
echo ""
echo "[Step 1] HRDP rn6 panel 다운로드"
PANEL_VCF="${PANEL_DIR}/HRDP_48smp_rn6_SNPs_PASS.vcf.gz"
if [ ! -f "${PANEL_VCF}" ]; then
    wget -O "${PANEL_VCF}" \
        "https://download.rgd.mcw.edu/strain_specific_variants/Dwinell_MCW_HybridRatDiversityProgram/Dec2021/Rnor6/HRDP_48smp_HPJoint_gatk4_rn6_SNPs_HF_PASS.vcf.gz"
    wget -q -O "${PANEL_VCF}.tbi" \
        "https://download.rgd.mcw.edu/strain_specific_variants/Dwinell_MCW_HybridRatDiversityProgram/Dec2021/Rnor6/HRDP_48smp_HPJoint_gatk4_rn6_SNPs_HF_PASS.vcf.gz.tbi"
    echo "  Downloaded: $(ls -lh ${PANEL_VCF} | awk '{print $5}')"
else
    echo "  [SKIP] 이미 다운로드됨"
fi
echo "[Step 1] 완료 — $(date)"

# ─── Step 2: Phase reference panel (per-chr) ─────────────────────────────────
echo ""
echo "[Step 2] Reference panel phasing (per-chr)"
for CHR in "${CHROMS[@]}"; do
    PHASED="${PANEL_DIR}/phased/hrdp_rn6_${CHR}_phased.vcf.gz"
    [ -f "${PHASED}" ] && continue

    echo "  Phasing chr${CHR}..."
    # HRDP panel uses "chr" prefix → strip to match rn6 naming (1,2,3...)
    # biallelic SNP만 추출 + missing → 0/0 + chr prefix 제거
    CHR_MAP="${PANEL_DIR}/phased/chr_strip.txt"
    if [ ! -f "${CHR_MAP}" ]; then
        for i in $(seq 1 20); do echo "chr${i} ${i}"; done > "${CHR_MAP}"
        echo "chrX X" >> "${CHR_MAP}"
        echo "chrY Y" >> "${CHR_MAP}"
        echo "chrM MT" >> "${CHR_MAP}"
    fi
    bcftools view -r "chr${CHR}" -v snps -m2 -M2 "${PANEL_VCF}" | \
        bcftools annotate --rename-chrs "${CHR_MAP}" | \
        awk 'BEGIN{OFS="\t"} /^#/{print;next} {for(i=10;i<=NF;i++){if($i ~ /^\.\/.:/){sub(/^\.\/\./,"0/0",$i)}else if($i=="./.")$i="0/0"};print}' | \
        bcftools view -Oz -o "${PANEL_DIR}/phased/hrdp_rn6_${CHR}_filled.vcf.gz"
    bcftools index -t "${PANEL_DIR}/phased/hrdp_rn6_${CHR}_filled.vcf.gz"

    java -Xmx4g -jar "${BEAGLE}" \
        gt="${PANEL_DIR}/phased/hrdp_rn6_${CHR}_filled.vcf.gz" \
        chrom="${CHR}" \
        out="${PANEL_DIR}/phased/hrdp_rn6_${CHR}_phased" \
        nthreads=4 2>/dev/null

    rm -f "${PANEL_DIR}/phased/hrdp_rn6_${CHR}_filled.vcf.gz" \
          "${PANEL_DIR}/phased/hrdp_rn6_${CHR}_filled.vcf.gz.tbi"
done
echo "[Step 2] 완료 — $(date)"

# ─── Step 3: Prepare RNA-seq VCFs (hard call, per-chr, bgzip) ────────────────
echo ""
echo "[Step 3] RNA-seq VCF 준비 (per-chr split + bgzip)"
RNA_GZ_DIR="${IMP_DIR}/rna_input"
mkdir -p "${RNA_GZ_DIR}"

for SRR in "${SAMPLES[@]}"; do
    FILTERED="${RNA_DIR}/${SRR}.filtered.vcf"
    [ ! -f "${FILTERED}" ] && continue

    for CHR in "${CHROMS[@]}"; do
        CHR_VCF="${RNA_GZ_DIR}/${SRR}_${CHR}.vcf.gz"
        [ -f "${CHR_VCF}" ] && continue

        # PASS SNP만, biallelic
        grep "^#\|^${CHR}[[:space:]]" "${FILTERED}" | \
            awk '$0 ~ /^#/ || $7=="PASS"' | \
            bcftools view -v snps -m2 -M2 -Oz -o "${CHR_VCF}" 2>/dev/null
        bcftools index -t "${CHR_VCF}" 2>/dev/null
    done
done
echo "[Step 3] 완료 — $(date)"

# ─── Step 4: Beagle imputation (per-chr, parallel) ──────────────────────────
echo ""
echo "[Step 4] Beagle imputation"

N_PARALLEL=8
PIDS=()

for CHR in "${CHROMS[@]}"; do
    PHASED="${PANEL_DIR}/phased/hrdp_rn6_${CHR}_phased.vcf.gz"
    [ ! -f "${PHASED}" ] && continue

    for SRR in "${SAMPLES[@]}"; do
        OUT="${IMPUTED_DIR}/${SRR}_${CHR}"
        [ -f "${OUT}.vcf.gz" ] && continue

        INPUT="${RNA_GZ_DIR}/${SRR}_${CHR}.vcf.gz"
        [ ! -f "${INPUT}" ] && continue

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
                gt="${INPUT}" \
                chrom="${CHR}" \
                out="${OUT}" \
                nthreads=2 2>/dev/null
        ) &
        PIDS+=("$!")
    done
done
echo "  모든 job 제출, 대기 중..."
wait
echo "[Step 4] 완료 — $(date)"

# ─── Step 5: Merge + Evaluate ───────────────────────────────────────────────
echo ""
echo "[Step 5] 평가 (DNA WES ground truth 대비)"

RESULT="${EVAL_DIR}/imputation_results.tsv"
echo -e "SRR\tTumor\tWES_SNPs\tBefore\tAfter\tCov_Before(%)\tCov_After(%)" > "${RESULT}"

printf "%-14s | %-6s | %8s | %8s | %8s | %7s | %7s\n" "SRR" "Tumor" "WES_SNPs" "Before" "After" "Before%" "After%"
printf "%s\n" "---------------|--------|----------|----------|----------|---------|--------"

for SRR in $(echo "${SAMPLES[@]}" | tr ' ' '\n' | sort); do
    WES_VCF_NAME="${SAMPLE_MAP[${SRR}]}"
    WES_VCF="${GT_DIR}/${WES_VCF_NAME}"
    TUMOR=$(echo "${WES_VCF_NAME}" | sed 's/GSM[0-9]*_//;s/_tumor.pass.vcf//')

    # WES VCF → bgzip + index (SNP only)
    TMPD=$(mktemp -d)
    grep "^#\|^[0-9XY]" "${WES_VCF}" | bcftools view -v snps -m2 -M2 -Oz -o "${TMPD}/wes.vcf.gz" 2>/dev/null
    bcftools index -t "${TMPD}/wes.vcf.gz" 2>/dev/null
    WES_TOTAL=$(bcftools view -H "${TMPD}/wes.vcf.gz" | wc -l)

    # Before: RNA-seq PASS SNP vs WES
    RNA_BEFORE="${RNA_DIR}/${SRR}.filtered.vcf"
    bcftools view -v snps -m2 -M2 "${RNA_BEFORE}" 2>/dev/null | \
        awk '$0 ~ /^#/ || $7=="PASS"' | \
        bcftools view -Oz -o "${TMPD}/rna.vcf.gz" 2>/dev/null
    bcftools index -t "${TMPD}/rna.vcf.gz" 2>/dev/null
    bcftools isec -p "${TMPD}/b" "${TMPD}/wes.vcf.gz" "${TMPD}/rna.vcf.gz" 2>/dev/null
    BEFORE=$(grep -vc "^#" "${TMPD}/b/0002.vcf" 2>/dev/null || echo 0)

    # After: imputed merge → ALT vs WES
    CHR_FILES=()
    for CHR in "${CHROMS[@]}"; do
        F="${IMPUTED_DIR}/${SRR}_${CHR}.vcf.gz"
        [ -f "${F}" ] && CHR_FILES+=("${F}")
    done

    if [ ${#CHR_FILES[@]} -gt 0 ]; then
        for F in "${CHR_FILES[@]}"; do
            [ ! -f "${F}.tbi" ] && bcftools index -t "${F}" 2>/dev/null
        done
        bcftools concat -a "${CHR_FILES[@]}" -Oz -o "${TMPD}/imputed.vcf.gz" 2>/dev/null
        bcftools index -t "${TMPD}/imputed.vcf.gz" 2>/dev/null
        bcftools view -i 'GT!="0|0" && GT!="0/0"' "${TMPD}/imputed.vcf.gz" -Oz -o "${TMPD}/imp_alt.vcf.gz" 2>/dev/null
        bcftools index -t "${TMPD}/imp_alt.vcf.gz" 2>/dev/null
        bcftools isec -p "${TMPD}/a" "${TMPD}/wes.vcf.gz" "${TMPD}/imp_alt.vcf.gz" 2>/dev/null
        AFTER=$(grep -vc "^#" "${TMPD}/a/0002.vcf" 2>/dev/null || echo 0)
    else
        AFTER=0
    fi

    rm -rf "${TMPD}"

    PCT_B=$(awk "BEGIN {printf \"%.1f\", ${BEFORE}/${WES_TOTAL}*100}")
    PCT_A=$(awk "BEGIN {printf \"%.1f\", ${AFTER}/${WES_TOTAL}*100}")

    printf "%-14s | %-6s | %8d | %8d | %8d | %6s%% | %6s%%\n" "${SRR}" "${TUMOR}" "${WES_TOTAL}" "${BEFORE}" "${AFTER}" "${PCT_B}" "${PCT_A}"
    echo -e "${SRR}\t${TUMOR}\t${WES_TOTAL}\t${BEFORE}\t${AFTER}\t${PCT_B}\t${PCT_A}" >> "${RESULT}"
done

echo ""
echo "============================================"
echo "전체 완료: $(date)"
echo "결과: ${RESULT}"
echo "Disk: $(df -h /home/yusanghyeon | tail -1 | awk '{print $4}') free"
echo "============================================"
