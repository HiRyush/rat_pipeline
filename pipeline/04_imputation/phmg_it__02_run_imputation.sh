#!/bin/bash
set -euo pipefail

# ============================================================
# RNA-seq Genotype Imputation Pipeline
# ============================================================
# Reference panel: HRDP (75 strains, rn7)
# Input: RNA-seq BAMs (15 PHMG_IT samples)
# Output: Imputed genotypes evaluated against DNA WGS ground truth
# ============================================================

eval "$($HOME/RAT_project/miniforge3/bin/conda shell.bash hook 2>/dev/null)"
conda activate rnaseq

WORK="/home/yusanghyeon/RAT_project/PHMG_IT/imputation"
REF="/home/yusanghyeon/RAT_project/PHMG_IT/reference/rn7.fa"
BEAGLE="${WORK}/tools/beagle.jar"
HRDP_DIR="${WORK}/ref_panel/hrdp_vcf"
BAM_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/mutect2/markdup"
GT_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/ground_truth"
RNA_GL_DIR="${WORK}/rna_gl"
PANEL_DIR="${WORK}/ref_panel"
IMPUTED_DIR="${WORK}/imputed"
EVAL_DIR="${WORK}/evaluation"

mkdir -p "${RNA_GL_DIR}" "${PANEL_DIR}/merged" "${PANEL_DIR}/phased" \
         "${IMPUTED_DIR}" "${EVAL_DIR}"

SAMPLES=(C1 C2 C3 C4 C5 P1 P2 P3 P4 P5 P6 P7 P8 P9 P10)

# rn7 major chromosomes
CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10
        chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20
        chrX)

# chromosome rename map (HRDP: 1→chr1)
CHR_MAP="${WORK}/ref_panel/chr_rename.txt"
if [ ! -f "${CHR_MAP}" ]; then
    for i in $(seq 1 20); do echo "${i} chr${i}"; done > "${CHR_MAP}"
    echo "X chrX" >> "${CHR_MAP}"
    echo "Y chrY" >> "${CHR_MAP}"
    echo "MT chrM" >> "${CHR_MAP}"
fi

echo "============================================"
echo "RNA-seq Genotype Imputation Pipeline"
echo "Started: $(date)"
echo "============================================"

# ─── Step 1: Merge HRDP VCFs into multi-sample reference panel (per-chr) ─────
echo ""
echo "[Step 1] HRDP reference panel 구축 (per-chromosome merge + chr rename)"

for CHR in "${CHROMS[@]}"; do
    MERGED="${PANEL_DIR}/merged/hrdp_${CHR}.vcf.gz"
    [ -f "${MERGED}" ] && [ -f "${MERGED}.tbi" ] && continue

    # HRDP chromosome name (without "chr" prefix)
    HRDP_CHR="${CHR#chr}"

    echo "  Merging ${CHR}..."

    # 각 HRDP VCF에서 해당 chromosome 추출 → 임시 파일들
    TMPDIR=$(mktemp -d "${WORK}/ref_panel/tmp_${CHR}_XXXX")
    VCF_LIST="${TMPDIR}/vcf_list.txt"
    > "${VCF_LIST}"

    for VCF in "${HRDP_DIR}"/*_rnBN7_SNPs_PASS.vcf.gz; do
        STRAIN=$(basename "${VCF}" _rnBN7_SNPs_PASS.vcf.gz)
        OUT_VCF="${TMPDIR}/${STRAIN}.vcf.gz"

        # 해당 chromosome 추출 + biallelic SNP만 + chr prefix 추가
        bcftools view -r "${HRDP_CHR}" "${VCF}" 2>/dev/null | \
            bcftools annotate --rename-chrs "${CHR_MAP}" 2>/dev/null | \
            bcftools view -v snps -m2 -M2 -O z -o "${OUT_VCF}" 2>/dev/null && \
        bcftools index -t "${OUT_VCF}" 2>/dev/null && \
        echo "${OUT_VCF}" >> "${VCF_LIST}" || true
    done

    N_VCF=$(wc -l < "${VCF_LIST}")
    if [ "${N_VCF}" -gt 0 ]; then
        bcftools merge -l "${VCF_LIST}" -O z -o "${MERGED}" 2>/dev/null
        bcftools index -t "${MERGED}" 2>/dev/null
        echo "    → ${CHR}: ${N_VCF} strains merged"
    fi
    rm -rf "${TMPDIR}"
done
echo "[Step 1] 완료 — $(date)"

# ─── Step 2: Extract genotype likelihoods from RNA-seq BAMs ──────────────────
echo ""
echo "[Step 2] RNA-seq BAM에서 genotype likelihood 추출"

for SAMPLE in "${SAMPLES[@]}"; do
    GL_VCF="${RNA_GL_DIR}/${SAMPLE}_rna_gl.vcf.gz"
    [ -f "${GL_VCF}" ] && continue

    BAM="${BAM_DIR}/${SAMPLE}.dedup.bam"
    echo "  [START] ${SAMPLE}"

    # 모든 HRDP variant 위치에서 allele depth 추출
    # -I: indel 무시 (SNP만), -a: AD,DP,INFO/AD 포함
    # -T: target sites (reference panel의 variant 위치)
    bcftools mpileup \
        -f "${REF}" \
        -I \
        -a "FORMAT/AD,FORMAT/DP" \
        -R <(for CHR in "${CHROMS[@]}"; do
                 M="${PANEL_DIR}/merged/hrdp_${CHR}.vcf.gz"
                 [ -f "${M}" ] && bcftools query -f '%CHROM\t%POS\n' "${M}" 2>/dev/null
             done | sort -k1,1V -k2,2n | awk '{print $1"\t"$2-1"\t"$2}') \
        "${BAM}" 2>/dev/null | \
    bcftools call -m -Oz -o "${GL_VCF}" 2>/dev/null

    bcftools index -t "${GL_VCF}" 2>/dev/null
    N_VAR=$(bcftools view -H "${GL_VCF}" 2>/dev/null | wc -l)
    echo "  [DONE] ${SAMPLE}: ${N_VAR} sites with GL"
done
echo "[Step 2] 완료 — $(date)"

# ─── Step 3: Phase reference panel + Impute ──────────────────────────────────
echo ""
echo "[Step 3] Beagle phasing + imputation (per-chromosome)"

for CHR in "${CHROMS[@]}"; do
    PANEL="${PANEL_DIR}/merged/hrdp_${CHR}.vcf.gz"
    [ ! -f "${PANEL}" ] && continue

    for SAMPLE in "${SAMPLES[@]}"; do
        IMPUTED_VCF="${IMPUTED_DIR}/${SAMPLE}_${CHR}"
        [ -f "${IMPUTED_VCF}.vcf.gz" ] && continue

        GL_VCF="${RNA_GL_DIR}/${SAMPLE}_rna_gl.vcf.gz"

        # Beagle: ref=reference panel, gt=target with GL
        # Beagle does phasing of ref internally if unphased
        java -Xmx4g -jar "${BEAGLE}" \
            ref="${PANEL}" \
            gt=<(bcftools view -r "${CHR}" "${GL_VCF}" 2>/dev/null) \
            chrom="${CHR}" \
            out="${IMPUTED_VCF}" \
            2>/dev/null || true
    done
    echo "  ${CHR} done"
done
echo "[Step 3] 완료 — $(date)"

# ─── Step 4: Merge imputed chromosomes per sample ────────────────────────────
echo ""
echo "[Step 4] 샘플별 chromosome merge"

for SAMPLE in "${SAMPLES[@]}"; do
    FINAL="${IMPUTED_DIR}/${SAMPLE}_imputed.vcf.gz"
    [ -f "${FINAL}" ] && continue

    INPUT_ARGS=()
    for CHR in "${CHROMS[@]}"; do
        CHR_VCF="${IMPUTED_DIR}/${SAMPLE}_${CHR}.vcf.gz"
        [ -f "${CHR_VCF}" ] && INPUT_ARGS+=("${CHR_VCF}")
    done

    if [ ${#INPUT_ARGS[@]} -gt 0 ]; then
        bcftools concat "${INPUT_ARGS[@]}" -O z -o "${FINAL}" 2>/dev/null
        bcftools index -t "${FINAL}" 2>/dev/null
        echo "  ${SAMPLE}: $(bcftools view -H "${FINAL}" 2>/dev/null | wc -l) imputed variants"
    fi
done
echo "[Step 4] 완료 — $(date)"

# ─── Step 5: Evaluate against DNA WGS ground truth ──────────────────────────
echo ""
echo "[Step 5] DNA WGS 대비 정확도 평가"

# DNA sample name mapping
declare -A DNA_MAP=(
    [C1]=P5S40W-2 [C2]=P5S54W-2 [C3]=P5S54W-3 [C4]=P5S54W-6 [C5]=P5S54W-8
    [P1]=P5H35W-11 [P2]=P5H40W-2 [P3]=P5H40W-16 [P4]=P5H48W-12 [P5]=P5H49W-9
    [P6]=P5H49W-18 [P7]=P5H52W-8 [P8]=P5H54W-2 [P9]=P5H54W-19 [P10]=P5M44W-11
)

echo "Sample | DNA_SNPs | Before_Imputation | After_Imputation | Coverage_Before(%) | Coverage_After(%)" > "${EVAL_DIR}/imputation_results.tsv"

for SAMPLE in "${SAMPLES[@]}"; do
    DNA_VCF="${GT_DIR}/${SAMPLE}_dna_snps.vcf.gz"
    RNA_BEFORE="/home/yusanghyeon/RAT_project/PHMG_IT/results/mutect2/force_pass/${SAMPLE}.force_pass.vcf.gz"
    RNA_AFTER="${IMPUTED_DIR}/${SAMPLE}_imputed.vcf.gz"

    [ ! -f "${RNA_AFTER}" ] && continue

    DNA_TOTAL=$(bcftools view -H "${DNA_VCF}" 2>/dev/null | wc -l)

    # Before imputation: RNA PASS vs DNA
    TMPDIR_B=$(mktemp -d)
    bcftools isec -p "${TMPDIR_B}" "${DNA_VCF}" "${RNA_BEFORE}" 2>/dev/null
    BEFORE=$(grep -v "^#" "${TMPDIR_B}/0002.vcf" 2>/dev/null | wc -l)
    rm -rf "${TMPDIR_B}"

    # After imputation: imputed vs DNA
    # 먼저 imputed VCF에서 ALT allele가 있는 것만 추출
    TMPDIR_A=$(mktemp -d)
    bcftools view -i 'GT="alt"' "${RNA_AFTER}" -O z -o "${TMPDIR_A}/imputed_alt.vcf.gz" 2>/dev/null
    bcftools index -t "${TMPDIR_A}/imputed_alt.vcf.gz" 2>/dev/null
    bcftools isec -p "${TMPDIR_A}/isec" "${DNA_VCF}" "${TMPDIR_A}/imputed_alt.vcf.gz" 2>/dev/null
    AFTER=$(grep -v "^#" "${TMPDIR_A}/isec/0002.vcf" 2>/dev/null | wc -l)
    rm -rf "${TMPDIR_A}"

    PCT_B=$(awk "BEGIN {printf \"%.1f\", ${BEFORE}/${DNA_TOTAL}*100}")
    PCT_A=$(awk "BEGIN {printf \"%.1f\", ${AFTER}/${DNA_TOTAL}*100}")

    echo "${SAMPLE} | ${DNA_TOTAL} | ${BEFORE} | ${AFTER} | ${PCT_B} | ${PCT_A}"
    echo "${SAMPLE} | ${DNA_TOTAL} | ${BEFORE} | ${AFTER} | ${PCT_B} | ${PCT_A}" >> "${EVAL_DIR}/imputation_results.tsv"
done

echo ""
echo "============================================"
echo "전체 완료: $(date)"
echo "결과: ${EVAL_DIR}/imputation_results.tsv"
echo "Disk: $(df -h /home/yusanghyeon | tail -1 | awk '{print $4}') free"
echo "============================================"
