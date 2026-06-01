#!/bin/bash
set -euo pipefail

# Step 3: 13 mammary samples (M1 early + M12+M20 late)에 MuTect2 병렬 실행
# 리소스 안전장치:
#   - 동시 4 samples (각 8GB heap + 2 native threads = 32 GB + 8 thread)
#   - 시스템 잔여: VM (8 thread, 17GB) + 일반 (~10GB) → 충분
#   - 13 samples → 4 batches (4+4+4+1)

WORK="/home/yusanghyeon/RAT_project/mammary_cancer"
SCRIPT="${WORK}/scripts/02_mutect2_per_sample.sh"
PARALLEL=5

# Sample list: M1 (early, 6) + M12+M20 (late, 7) = 13
# Based on sample_mapping_verified.tsv
SAMPLES=(
  # M1 (1 month) — 6 samples
  SRR33625170  # M1_1
  SRR33625169  # M1_2
  SRR33625168  # M1_3
  SRR33625167  # M1_5
  SRR33625166  # M1_6
  SRR33625165  # M1_7
  # M12 (12 month) — 3 samples
  SRR33625154  # M12_1
  SRR33625153  # M12_2
  SRR33625152  # M12_3
  # M20 (20 month) — 4 samples
  SRR33625151  # M20_1
  SRR33625150  # M20_2
  SRR33625149  # M20_3
  SRR33625148  # M20_4
)

echo "[$(date +%H:%M:%S)] Starting parallel MuTect2 on ${#SAMPLES[@]} samples (${PARALLEL} parallel)"
echo "Sample list:"
printf '  %s\n' "${SAMPLES[@]}"

# parallel execution with backpressure
printf '%s\n' "${SAMPLES[@]}" | xargs -I {} -P "${PARALLEL}" bash "$SCRIPT" {}

echo "[$(date +%H:%M:%S)] All samples completed."

# Summary
echo ""
echo "=== Summary ==="
for SRR in "${SAMPLES[@]}"; do
  RAW="${WORK}/results/mutect2/raw/${SRR}.mutect2.vcf.gz"
  if [ -f "$RAW" ]; then
    COUNT=$(zcat "$RAW" | grep -v "^#" | wc -l)
    SIZE=$(du -h "$RAW" | cut -f1)
    echo "  ${SRR}: ${COUNT} variants (${SIZE})"
  else
    echo "  ${SRR}: MISSING!"
  fi
done
