#!/bin/bash
set -euo pipefail

# Step 2: MuTect2 per-sample variant calling (exonic only, tumor-only mode)
# Usage: ./02_mutect2_per_sample.sh <SRR_ID>
#
# 리소스 cap:
#   - Java heap: 8 GB
#   - Native threads: 2
#   - Total per process: ~8-10 GB RAM, 2-4 CPU threads
#
# PHMG_IT과의 차이:
#   - Exonic intervals only (WES validation 영역 맞춤)
#   - Tumor-only mode (matched normal 없음)
#   - 동일: GATK 4.6.1.0, --annotations-to-exclude TandemRepeat,
#         --disable-read-filter MappingQualityAvailableReadFilter

SRR=$1
WORK="/home/yusanghyeon/RAT_project/mammary_cancer"
REF="${WORK}/reference/rn6/Rnor_6.0.fa"
BAM="${WORK}/results/dedup_bam/${SRR}.dedup.bam"
INTERVALS="${WORK}/results/mutect2/intervals/exonic_rn6.interval_list"
OUT_RAW="${WORK}/results/mutect2/raw/${SRR}.mutect2.vcf.gz"
OUT_F1R2="${WORK}/results/mutect2/raw/${SRR}.f1r2.tar.gz"
LOG="${WORK}/results/mutect2/logs/${SRR}.mutect2.log"
TMPDIR="/mnt/data/mammary_tmp/${SRR}"

mkdir -p "$TMPDIR"

source /home/yusanghyeon/RAT_project/miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

echo "[$(date +%H:%M:%S)] [${SRR}] Starting MuTect2 (exonic, tumor-only)..." > "$LOG"

gatk --java-options "-Xmx8g -XX:ParallelGCThreads=2 -Djava.io.tmpdir=${TMPDIR}" Mutect2 \
  -R "$REF" \
  -I "$BAM" \
  -L "$INTERVALS" \
  --native-pair-hmm-threads 2 \
  --annotations-to-exclude TandemRepeat \
  --disable-read-filter MappingQualityAvailableReadFilter \
  --f1r2-tar-gz "$OUT_F1R2" \
  -O "$OUT_RAW" \
  --tmp-dir "$TMPDIR" \
  >> "$LOG" 2>&1

echo "[$(date +%H:%M:%S)] [${SRR}] MuTect2 done. Variants in raw VCF:" >> "$LOG"
zcat "$OUT_RAW" | grep -v "^#" | wc -l >> "$LOG"

# 정리
rm -rf "$TMPDIR"

echo "[$(date +%H:%M:%S)] [${SRR}] DONE"
