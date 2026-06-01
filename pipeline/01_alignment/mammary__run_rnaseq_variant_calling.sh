#!/bin/bash
set -euo pipefail

# RNA-seq Variant Calling Pipeline (Mammary Cancer / NMU Rat Model)
# Reference: rn6 (Rnor_6.0) — to match Strelka ground truth VCF coordinates
#
# Key differences from GATK default best practices:
#   - SplitNCigarReads REMOVED: causes -54% coverage loss (PHMG Phase 1 finding)
#   - STAR Arm C params: aggressive multi-mapping for maximum coverage
#   - HaplotypeCaller: soft-clipped bases ALLOWED (consistent with Arm C)
#
# Usage: ./run_rnaseq_variant_calling.sh <SRR_ID> [THREADS]

SRR=$1
THREADS=${2:-16}

WORK="/home/yusanghyeon/RAT_project/mammary_cancer"
REF="${WORK}/reference/rn6/Rnor_6.0.fa"
GTF="${WORK}/reference/rn6/Rnor_6.0.104.gtf"
STAR_INDEX="${WORK}/reference/star_rn6"
KNOWN_SITES="${WORK}/known_sites/rattus_norvegicus.vcf.gz"
FASTQ_DIR="${WORK}/fastq"
RESULT_DIR="${WORK}/results/rn6"
TMPDIR="${WORK}/tmp/${SRR}"

mkdir -p "${TMPDIR}" "${RESULT_DIR}"

FASTQ_R1="${FASTQ_DIR}/${SRR}_1.fastq.gz"
FASTQ_R2="${FASTQ_DIR}/${SRR}_2.fastq.gz"

if [ ! -f "${FASTQ_R1}" ] || [ ! -f "${FASTQ_R2}" ]; then
    echo "ERROR: FASTQ files not found for ${SRR}"
    exit 1
fi

echo "============================================"
echo "[$(date '+%H:%M:%S')] Starting ${SRR}"
echo "============================================"

# ─── Step 1: STAR 2-pass alignment (Arm C: aggressive multi-mapping) ───────────
echo "[$(date '+%H:%M:%S')] Step 1: STAR alignment"
STAR --genomeDir "${STAR_INDEX}" \
     --readFilesIn "${FASTQ_R1}" "${FASTQ_R2}" \
     --readFilesCommand zcat \
     --outFileNamePrefix "${TMPDIR}/" \
     --outSAMtype BAM SortedByCoordinate \
     --twopassMode Basic \
     --outSAMattrRGline "ID:${SRR}" "SM:${SRR}" "PL:ILLUMINA" \
     --runThreadN "${THREADS}" \
     --outBAMsortingThreadN 4 \
     --limitBAMsortRAM 30000000000 \
     --outFilterMultimapNmax 50 \
     --outFilterScoreMinOverLread 0.33 \
     --outFilterMatchNminOverLread 0.33 \
     --alignEndsType Local \
     --winAnchorMultimapNmax 100 \
     --alignIntronMax 1000000 \
     --alignMatesGapMax 1000000
echo "[$(date '+%H:%M:%S')] Step 1 done"

# ─── Step 2: Mark duplicates ───────────────────────────────────────────────────
echo "[$(date '+%H:%M:%S')] Step 2: MarkDuplicates"
gatk MarkDuplicates \
     -I "${TMPDIR}/Aligned.sortedByCoord.out.bam" \
     -O "${TMPDIR}/${SRR}.dedup.bam" \
     -M "${TMPDIR}/${SRR}.metrics.txt" \
     --TMP_DIR "${TMPDIR}"

rm -f "${TMPDIR}/Aligned.sortedByCoord.out.bam"
samtools index "${TMPDIR}/${SRR}.dedup.bam"
echo "[$(date '+%H:%M:%S')] Step 2 done"

# NOTE: SplitNCigarReads intentionally SKIPPED
#       Reason: causes -54% DP≥5 coverage loss (hard-clips junction overhangs)
#       Validated in PHMG project Phase 1 (Arm E-pre vs E-post comparison)

# ─── Step 3: BQSR ─────────────────────────────────────────────────────────────
echo "[$(date '+%H:%M:%S')] Step 3: BQSR"
gatk BaseRecalibrator \
     -R "${REF}" \
     -I "${TMPDIR}/${SRR}.dedup.bam" \
     --known-sites "${KNOWN_SITES}" \
     -O "${TMPDIR}/${SRR}.recal.table" \
     --tmp-dir "${TMPDIR}"

gatk ApplyBQSR \
     -R "${REF}" \
     -I "${TMPDIR}/${SRR}.dedup.bam" \
     --bqsr-recal-file "${TMPDIR}/${SRR}.recal.table" \
     -O "${TMPDIR}/${SRR}.bqsr.bam" \
     --tmp-dir "${TMPDIR}"

rm -f "${TMPDIR}/${SRR}.dedup.bam" "${TMPDIR}/${SRR}.dedup.bam.bai"
echo "[$(date '+%H:%M:%S')] Step 3 done"

# ─── Step 4: HaplotypeCaller ───────────────────────────────────────────────────
echo "[$(date '+%H:%M:%S')] Step 4: HaplotypeCaller"
gatk HaplotypeCaller \
     -R "${REF}" \
     -I "${TMPDIR}/${SRR}.bqsr.bam" \
     -O "${TMPDIR}/${SRR}.gatk.vcf" \
     --standard-min-confidence-threshold-for-calling 10 \
     --min-base-quality-score 10 \
     --native-pair-hmm-threads "${THREADS}" \
     --disable-read-filter MappingQualityAvailableReadFilter \
     --minimum-mapping-quality 0 \
     --tmp-dir "${TMPDIR}"
echo "[$(date '+%H:%M:%S')] Step 4 done"

# ─── Step 5: Variant Filtration ───────────────────────────────────────────────
echo "[$(date '+%H:%M:%S')] Step 5: VariantFiltration"
gatk VariantFiltration \
     -R "${REF}" \
     -V "${TMPDIR}/${SRR}.gatk.vcf" \
     -O "${RESULT_DIR}/${SRR}.filtered.vcf" \
     --filter-name "QD" --filter-expression "QD < 2.0" \
     --filter-name "FS" --filter-expression "FS > 30.0" \
     --filter-name "DP" --filter-expression "DP < 5" \
     --tmp-dir "${TMPDIR}"
echo "[$(date '+%H:%M:%S')] Step 5 done"

# Save raw VCF as well
cp "${TMPDIR}/${SRR}.gatk.vcf" "${RESULT_DIR}/${SRR}.raw.vcf" 2>/dev/null || true

# Cleanup
rm -rf "${TMPDIR}"

echo "============================================"
echo "[$(date '+%H:%M:%S')] ${SRR} COMPLETED"
echo "============================================"
