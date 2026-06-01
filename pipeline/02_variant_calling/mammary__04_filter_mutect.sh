#!/bin/bash
set -euo pipefail

# Step 4: FilterMutectCalls + PASS variant 추출
# Per-sample 처리

WORK="/home/yusanghyeon/RAT_project/mammary_cancer"
REF="${WORK}/reference/rn6/Rnor_6.0.fa"
RAW_DIR="${WORK}/results/mutect2/raw"
FILT_DIR="${WORK}/results/mutect2/filtered"
LOG_DIR="${WORK}/results/mutect2/logs"

source /home/yusanghyeon/RAT_project/miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

SAMPLES=(
  SRR33625170 SRR33625169 SRR33625168 SRR33625167 SRR33625166 SRR33625165
  SRR33625154 SRR33625153 SRR33625152
  SRR33625151 SRR33625150 SRR33625149 SRR33625148
)

filter_one() {
  local SRR=$1
  local RAW="${RAW_DIR}/${SRR}.mutect2.vcf.gz"
  local F1R2="${RAW_DIR}/${SRR}.f1r2.tar.gz"
  local FILT="${FILT_DIR}/${SRR}.filtered.vcf.gz"
  local PASS="${FILT_DIR}/${SRR}.pass.vcf.gz"
  local OB_MODEL="${FILT_DIR}/${SRR}.read-orientation-model.tar.gz"
  local LOG="${LOG_DIR}/${SRR}.filter.log"
  local TMPDIR="/mnt/data/mammary_tmp/${SRR}_filter"
  mkdir -p "$TMPDIR"

  echo "[$(date +%H:%M:%S)] [${SRR}] LearnReadOrientationModel..." > "$LOG"
  gatk --java-options "-Xmx4g -XX:ParallelGCThreads=1 -Djava.io.tmpdir=${TMPDIR}" LearnReadOrientationModel \
    -I "$F1R2" -O "$OB_MODEL" >> "$LOG" 2>&1

  echo "[$(date +%H:%M:%S)] [${SRR}] FilterMutectCalls..." >> "$LOG"
  # PHMG_IT과 동일하게 --min-median-base-quality 0 사용 (GATK 4.6.1.0에서 MBQ가 항상 0,0으로 보고되는 quirk 우회)
  gatk --java-options "-Xmx4g -XX:ParallelGCThreads=1 -Djava.io.tmpdir=${TMPDIR}" FilterMutectCalls \
    -R "$REF" \
    -V "$RAW" \
    --ob-priors "$OB_MODEL" \
    --min-median-base-quality 0 \
    -O "$FILT" \
    --tmp-dir "$TMPDIR" >> "$LOG" 2>&1

  echo "[$(date +%H:%M:%S)] [${SRR}] PASS variant 추출..." >> "$LOG"
  bcftools view -f PASS "$FILT" -Oz -o "$PASS"
  bcftools index -t "$PASS"

  local PASS_COUNT=$(bcftools view -H "$PASS" | wc -l)
  echo "[$(date +%H:%M:%S)] [${SRR}] PASS variants: ${PASS_COUNT}" >> "$LOG"

  rm -rf "$TMPDIR"
  echo "[$(date +%H:%M:%S)] [${SRR}] DONE: ${PASS_COUNT} PASS"
}

export -f filter_one
export RAW_DIR FILT_DIR LOG_DIR REF

# 4 parallel으로 처리 (FilterMutectCalls는 가벼우므로 OK)
printf '%s\n' "${SAMPLES[@]}" | xargs -I {} -P 4 bash -c 'filter_one "$@"' _ {}

echo ""
echo "=== Filter Summary ==="
for SRR in "${SAMPLES[@]}"; do
  PASS="${FILT_DIR}/${SRR}.pass.vcf.gz"
  if [ -f "$PASS" ]; then
    COUNT=$(bcftools view -H "$PASS" | wc -l)
    echo "  ${SRR}: ${COUNT} PASS variants"
  fi
done
