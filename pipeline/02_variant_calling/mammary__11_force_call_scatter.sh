#!/bin/bash
# Mammary force-call scatter — PHMG_IT mutect2_scatter.sh Step 5 동등
# 13 samples × 21 chrs = 273 jobs, MAX_PAR 병렬.
set -uo pipefail  # NOTE: no -e, GATK 4.6.1.0 TandemRepeat bug → 비정상 종료지만 VCF는 완성

M=/home/yusanghyeon/RAT_project/mammary_cancer
REF=$M/reference/rn6/Rnor_6.0.fa
OUT=$M/results/mutect2/force_call
GATK=/home/yusanghyeon/RAT_project/miniforge3/envs/rnaseq/bin/gatk
SCRATCH=$OUT/scratch

source /home/yusanghyeon/RAT_project/miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

SAMPLES=(SRR33625170 SRR33625169 SRR33625168 SRR33625167 SRR33625166 SRR33625165 \
         SRR33625154 SRR33625153 SRR33625152 SRR33625151 SRR33625150 SRR33625149 SRR33625148)
CHRS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 X)

mkdir -p $OUT/force_scatter $OUT/logs $SCRATCH

MAX_PAR=20
_wait_slot() {
  while [ $(jobs -r | wc -l) -ge $MAX_PAR ]; do sleep 1; done
}

TOTAL=$((${#SAMPLES[@]} * ${#CHRS[@]}))
echo "[$(date +%H:%M:%S)] Force-call scatter starting: $TOTAL jobs, $MAX_PAR parallel"

for S in "${SAMPLES[@]}"; do
  BAM=$M/results/dedup_bam/${S}.dedup.bam
  for CHR in "${CHRS[@]}"; do
    OUT_VCF=$OUT/force_scatter/${S}_${CHR}.vcf.gz
    UNION=$OUT/union_sites_per_chr/${CHR}.vcf.gz
    if [ -f $OUT_VCF ] && [ -f ${OUT_VCF}.stats ]; then continue; fi
    _wait_slot
    (
      $GATK Mutect2 \
        -R $REF \
        -I $BAM \
        -tumor $S \
        -L $CHR \
        --alleles $UNION \
        --force-call-filtered-alleles \
        --callable-depth 1 \
        --native-pair-hmm-threads 1 \
        --disable-read-filter MappingQualityAvailableReadFilter \
        --minimum-mapping-quality 0 \
        --annotations-to-exclude TandemRepeat \
        --tmp-dir $SCRATCH \
        -O $OUT_VCF \
        > $OUT/logs/${S}_${CHR}.log 2>&1 || true
      # GATK 4.6.1.0 TandemRepeat bug: VCF 완성됐어도 .stats 누락 가능 → 빈 stats 생성
      if [ -f $OUT_VCF ] && [ ! -f ${OUT_VCF}.stats ]; then
        echo "stats placeholder due to TandemRepeat bug" > ${OUT_VCF}.stats
      fi
    ) &
  done
done

wait
echo "[$(date +%H:%M:%S)] All force-call jobs done."

# Check failures
FAIL=0
for S in "${SAMPLES[@]}"; do
  for CHR in "${CHRS[@]}"; do
    if [ ! -f $OUT/force_scatter/${S}_${CHR}.vcf.gz ]; then
      echo "  MISSING: ${S}_${CHR}"
      FAIL=$((FAIL+1))
    fi
  done
done
echo "  Missing: $FAIL"

# Gather per-sample
echo "[$(date +%H:%M:%S)] Gathering per-sample force-called VCFs..."
mkdir -p $OUT/force_called
for S in "${SAMPLES[@]}"; do
  OUT_GATHERED=$OUT/force_called/${S}.force.vcf.gz
  if [ -f $OUT_GATHERED ]; then echo "  $S already gathered"; continue; fi
  CHR_VCFS=()
  for CHR in "${CHRS[@]}"; do
    CHR_VCFS+=($OUT/force_scatter/${S}_${CHR}.vcf.gz)
  done
  bcftools concat -a "${CHR_VCFS[@]}" -Oz -o $OUT_GATHERED 2>$OUT/logs/gather_${S}.log
  bcftools index -t $OUT_GATHERED 2>/dev/null
  N=$(bcftools view -H $OUT_GATHERED 2>/dev/null | wc -l)
  echo "  $S: $N variants"
done
echo "[$(date +%H:%M:%S)] Done."
