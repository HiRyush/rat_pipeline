#!/bin/bash
set -euo pipefail

export PATH="/home/yusanghyeon/RAT_project/miniforge3/envs/rnaseq/bin:$PATH"

WORK="/home/yusanghyeon/RAT_project/mammary_cancer"
BEAGLE="/home/yusanghyeon/RAT_project/PHMG_IT/imputation/tools/beagle.jar"
REF="${WORK}/reference/rn6/Rnor_6.0.fa"
STAR_INDEX="${WORK}/reference/star_rn6"
FASTQ_DIR="${WORK}/fastq"
GT_DIR="${WORK}/ground_truth"
IMP_DIR="${WORK}/imputation"
PANEL_DIR="${IMP_DIR}/ref_panel/phased"
GL_DIR="${IMP_DIR}/rna_gl"
IMPUTED_DIR="${IMP_DIR}/imputed_gl"
EVAL_DIR="${IMP_DIR}/evaluation"
BAM_DIR="${IMP_DIR}/bam_tmp"
DEDUP_DIR="${WORK}/results/dedup_bam"

mkdir -p "${GL_DIR}" "${IMPUTED_DIR}" "${EVAL_DIR}" "${BAM_DIR}" "${DEDUP_DIR}"

CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 X)

# Verified sample mapping (SRA API, 2026-04-24)
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

SAMPLES=(SRR33625148 SRR33625149 SRR33625150 SRR33625151
         SRR33625152 SRR33625153 SRR33625154 SRR33625155
         SRR33625156 SRR33625157 SRR33625158 SRR33625159
         SRR33625160 SRR33625161 SRR33625162 SRR33625163
         SRR33625164 SRR33625165 SRR33625166 SRR33625167
         SRR33625168 SRR33625169 SRR33625170)

THREADS=6
BATCH_SIZE=1

echo "============================================"
echo "Mammary Cancer GL-based Imputation Pipeline"
echo "Samples: ${#SAMPLES[@]} | Batch: ${BATCH_SIZE}"
echo "Started: $(date)"
echo "============================================"

# ─── Process in batches: STAR → MarkDup → GL extract → delete BAM ───────────
echo ""
echo "[Phase 1] STAR + GL extraction (batch=${BATCH_SIZE}, BAM 즉시 삭제)"

DONE_SAMPLES=0
TOTAL=${#SAMPLES[@]}

for ((i=0; i<TOTAL; i+=BATCH_SIZE)); do
    BATCH=("${SAMPLES[@]:i:BATCH_SIZE}")
    echo ""
    echo "  === Batch $((i/BATCH_SIZE+1)): ${BATCH[*]} ==="

    # ── STAR alignment (parallel within batch) ──
    PIDS=()
    for SRR in "${BATCH[@]}"; do
        GL_DONE="${GL_DIR}/${SRR}_1_gl.vcf.gz"
        [ -f "${GL_DONE}" ] && { echo "  [SKIP] ${SRR} — GL already exists"; continue; }

        TMPDIR="${BAM_DIR}/${SRR}"
        mkdir -p "${TMPDIR}"

        BAM="${TMPDIR}/${SRR}.sorted.bam"
        if [ ! -f "${BAM}" ]; then
            echo "  [STAR] ${SRR}"
            (
                STAR \
                    --genomeDir "${STAR_INDEX}" \
                    --readFilesIn "${FASTQ_DIR}/${SRR}_1.fastq.gz" "${FASTQ_DIR}/${SRR}_2.fastq.gz" \
                    --readFilesCommand zcat \
                    --outFileNamePrefix "${TMPDIR}/" \
                    --outSAMtype BAM SortedByCoordinate \
                    --twopassMode Basic \
                    --outSAMattrRGline "ID:${SRR}" "SM:${SRR}" "PL:ILLUMINA" \
                    --runThreadN "${THREADS}" \
                    --outBAMsortingThreadN 4 \
                    --limitBAMsortRAM 15000000000 \
                    --outFilterMultimapNmax 50 \
                    --outFilterScoreMinOverLread 0.33 \
                    --outFilterMatchNminOverLread 0.33 \
                    --alignEndsType Local \
                    --winAnchorMultimapNmax 100 \
                    --alignIntronMax 1000000 \
                    --alignMatesGapMax 1000000 \
                    > "${TMPDIR}/star.log" 2>&1 && \
                mv "${TMPDIR}/Aligned.sortedByCoord.out.bam" "${BAM}" && \
                samtools index "${BAM}"
            ) &
            PIDS+=("$!:${SRR}")
        fi
    done

    # STAR 완료 대기
    for pid_info in "${PIDS[@]+"${PIDS[@]}"}"; do
        pid="${pid_info%%:*}"
        srr="${pid_info##*:}"
        wait "${pid}" && echo "  [STAR done] ${srr}" || echo "  [STAR FAIL] ${srr}"
    done

    # ── MarkDuplicates + GL extraction + BAM 삭제 (sequential per sample) ──
    for SRR in "${BATCH[@]}"; do
        GL_DONE="${GL_DIR}/${SRR}_1_gl.vcf.gz"
        [ -f "${GL_DONE}" ] && continue

        TMPDIR="${BAM_DIR}/${SRR}"
        BAM="${TMPDIR}/${SRR}.sorted.bam"
        DEDUP="${TMPDIR}/${SRR}.dedup.bam"

        [ ! -f "${BAM}" ] && { echo "  [SKIP] ${SRR} — no BAM"; continue; }

        # MarkDuplicates
        DEDUP_KEEP="${DEDUP_DIR}/${SRR}.dedup.bam"
        if [ -f "${DEDUP_KEEP}" ]; then
            DEDUP="${DEDUP_KEEP}"
            echo "  [SKIP MarkDup] ${SRR} — dedup BAM already exists"
        else
            echo "  [MarkDup] ${SRR}"
            gatk MarkDuplicates \
                -I "${BAM}" \
                -O "${DEDUP}" \
                -M "${TMPDIR}/markdup_metrics.txt" \
                --TMP_DIR "${TMPDIR}" \
                2>"${TMPDIR}/markdup.log"
            samtools index "${DEDUP}"
            # dedup BAM을 보존 디렉토리로 복사
            cp "${DEDUP}" "${DEDUP_KEEP}"
            cp "${DEDUP}.bai" "${DEDUP_KEEP}.bai" 2>/dev/null || \
                samtools index "${DEDUP_KEEP}"
            echo "  [SAVED] ${DEDUP_KEEP}"
        fi
        rm -f "${BAM}" "${BAM}.bai"  # sorted BAM만 삭제

        # GL extraction (per-chr parallel)
        echo "  [GL] ${SRR}"
        GL_PIDS=()
        for CHR in "${CHROMS[@]}"; do
            GL_CHR="${GL_DIR}/${SRR}_${CHR}_gl.vcf.gz"
            [ -f "${GL_CHR}" ] && continue

            PANEL="${PANEL_DIR}/hrdp_rn6_${CHR}_phased.vcf.gz"
            [ ! -f "${PANEL}" ] && continue

            (
                bcftools mpileup \
                    -f "${REF}" \
                    -T "${PANEL}" \
                    -I \
                    -a "FORMAT/AD,FORMAT/DP" \
                    --min-MQ 0 \
                    "${DEDUP}" 2>/dev/null | \
                bcftools call -m -Oz -o "${GL_CHR}" 2>/dev/null && \
                bcftools index -t "${GL_CHR}" 2>/dev/null
            ) &
            GL_PIDS+=("$!")
            # 최대 8 병렬
            while [ ${#GL_PIDS[@]} -ge 8 ]; do
                NEW=()
                for p in "${GL_PIDS[@]}"; do kill -0 "$p" 2>/dev/null && NEW+=("$p"); done
                GL_PIDS=("${NEW[@]+"${NEW[@]}"}")
                [ ${#GL_PIDS[@]} -ge 8 ] && sleep 2
            done
        done
        wait  # 모든 GL chr 완료

        # tmp 정리 (dedup BAM은 이미 보존됨, sorted BAM/STAR tmp만 삭제)
        rm -rf "${TMPDIR}"
        DONE_SAMPLES=$((DONE_SAMPLES + 1))
        echo "  [DONE] ${SRR} — GL extracted, BAM deleted (${DONE_SAMPLES}/${TOTAL})"
        echo "    Disk: $(df -h /home/yusanghyeon | tail -1 | awk '{print $4}') free"
    done
done
echo ""
echo "[Phase 1] 완료 — $(date)"

# ─── Phase 2: Beagle imputation (per-chr parallel) ──────────────────────────
echo ""
echo "[Phase 2] Beagle imputation"

N_PARALLEL=8
PIDS=()

for CHR in "${CHROMS[@]}"; do
    PHASED="${PANEL_DIR}/hrdp_rn6_${CHR}_phased.vcf.gz"
    [ ! -f "${PHASED}" ] && continue

    for SRR in "${SAMPLES[@]}"; do
        OUT="${IMPUTED_DIR}/${SRR}_${CHR}"
        [ -f "${OUT}.vcf.gz" ] && continue

        GL="${GL_DIR}/${SRR}_${CHR}_gl.vcf.gz"
        [ ! -f "${GL}" ] && continue

        while [ ${#PIDS[@]} -ge ${N_PARALLEL} ]; do
            NEW=()
            for p in "${PIDS[@]}"; do kill -0 "$p" 2>/dev/null && NEW+=("$p"); done
            PIDS=("${NEW[@]+"${NEW[@]}"}")
            [ ${#PIDS[@]} -ge ${N_PARALLEL} ] && sleep 2
        done

        (
            java -Xmx2g -jar "${BEAGLE}" \
                ref="${PHASED}" \
                gt="${GL}" \
                chrom="${CHR}" \
                out="${OUT}" \
                nthreads=2 2>/dev/null
        ) &
        PIDS+=("$!")
    done
done
echo "  모든 job 제출, 대기 중..."
wait
echo "[Phase 2] 완료 — $(date)"

# ─── Phase 3: Evaluate ──────────────────────────────────────────────────────
echo ""
echo "[Phase 3] 평가 (WES ground truth 대비)"

RESULT="${EVAL_DIR}/imputation_gl_results.tsv"
echo -e "SRR\tTumor\tWES_SNPs\tBefore\tAfter_HardCall\tAfter_GL\tCov_Before(%)\tCov_HC(%)\tCov_GL(%)" > "${RESULT}"

printf "%-14s | %-6s | %8s | %8s | %8s | %8s | %7s | %7s | %7s\n" \
    "SRR" "Tumor" "WES_SNPs" "Before" "After_HC" "After_GL" "Bef%" "HC%" "GL%"
printf "%s\n" "---------------|--------|----------|----------|----------|----------|---------|---------|--------"

for SRR in $(echo "${SAMPLES[@]}" | tr ' ' '\n' | sort); do
    WES_VCF_NAME="${SAMPLE_MAP[${SRR}]}"
    WES_VCF="${GT_DIR}/${WES_VCF_NAME}"
    TUMOR=$(echo "${WES_VCF_NAME}" | sed 's/GSM[0-9]*_//;s/_tumor.pass.vcf//')

    TMPD=$(mktemp -d)

    # WES → bgzip SNP only
    grep "^#\|^[0-9XY]" "${WES_VCF}" | bcftools view -v snps -m2 -M2 -Oz -o "${TMPD}/wes.vcf.gz" 2>/dev/null
    bcftools index -t "${TMPD}/wes.vcf.gz" 2>/dev/null
    WES_TOTAL=$(bcftools view -H "${TMPD}/wes.vcf.gz" | wc -l)

    # Before: RNA-seq PASS
    RNA_BEFORE="${WORK}/results/rn6/${SRR}.filtered.vcf"
    bcftools view -v snps -m2 -M2 "${RNA_BEFORE}" 2>/dev/null | \
        awk '$0 ~ /^#/ || $7=="PASS"' | \
        bcftools view -Oz -o "${TMPD}/rna.vcf.gz" 2>/dev/null
    bcftools index -t "${TMPD}/rna.vcf.gz" 2>/dev/null
    bcftools isec -p "${TMPD}/b" "${TMPD}/wes.vcf.gz" "${TMPD}/rna.vcf.gz" 2>/dev/null
    BEFORE=$(grep -vc "^#" "${TMPD}/b/0002.vcf" 2>/dev/null || echo 0)

    # After hard-call (from previous run)
    HC_RESULT=$(grep "^${SRR}" "${EVAL_DIR}/imputation_results.tsv" 2>/dev/null | cut -f4)
    [ -z "${HC_RESULT}" ] && HC_RESULT=0

    # After GL-based
    CHR_FILES=()
    for CHR in "${CHROMS[@]}"; do
        F="${IMPUTED_DIR}/${SRR}_${CHR}.vcf.gz"
        [ -f "${F}" ] && CHR_FILES+=("${F}")
    done

    AFTER_GL=0
    if [ ${#CHR_FILES[@]} -gt 0 ]; then
        for F in "${CHR_FILES[@]}"; do
            [ ! -f "${F}.tbi" ] && bcftools index -t "${F}" 2>/dev/null
        done
        bcftools concat -a "${CHR_FILES[@]}" -Oz -o "${TMPD}/imputed.vcf.gz" 2>/dev/null
        bcftools index -t "${TMPD}/imputed.vcf.gz" 2>/dev/null
        bcftools view -i 'GT!="0|0" && GT!="0/0"' "${TMPD}/imputed.vcf.gz" -Oz -o "${TMPD}/alt.vcf.gz" 2>/dev/null
        bcftools index -t "${TMPD}/alt.vcf.gz" 2>/dev/null
        bcftools isec -p "${TMPD}/a" "${TMPD}/wes.vcf.gz" "${TMPD}/alt.vcf.gz" 2>/dev/null
        AFTER_GL=$(grep -vc "^#" "${TMPD}/a/0002.vcf" 2>/dev/null || echo 0)
    fi

    rm -rf "${TMPD}"

    PCT_B=$(awk "BEGIN {printf \"%.1f\", ${BEFORE}/${WES_TOTAL}*100}")
    PCT_HC=$(awk "BEGIN {printf \"%.1f\", ${HC_RESULT}/${WES_TOTAL}*100}")
    PCT_GL=$(awk "BEGIN {printf \"%.1f\", ${AFTER_GL}/${WES_TOTAL}*100}")

    printf "%-14s | %-6s | %8d | %8d | %8d | %8d | %6s%% | %6s%% | %6s%%\n" \
        "${SRR}" "${TUMOR}" "${WES_TOTAL}" "${BEFORE}" "${HC_RESULT}" "${AFTER_GL}" "${PCT_B}" "${PCT_HC}" "${PCT_GL}"
    echo -e "${SRR}\t${TUMOR}\t${WES_TOTAL}\t${BEFORE}\t${HC_RESULT}\t${AFTER_GL}\t${PCT_B}\t${PCT_HC}\t${PCT_GL}" >> "${RESULT}"
done

echo ""
echo "============================================"
echo "전체 완료: $(date)"
echo "결과: ${RESULT}"
echo "Disk: $(df -h /home/yusanghyeon | tail -1 | awk '{print $4}') free"
echo "============================================"
