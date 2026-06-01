#!/bin/bash
# Mammary force-called VCF filter + differential (PHMG_IT mutect2_scatter Step 7-8 동등)
set -uo pipefail

M=/home/yusanghyeon/RAT_project/mammary_cancer
REF=$M/reference/rn6/Rnor_6.0.fa
OUT=$M/results/mutect2/force_call
COWORKER=/home/yusanghyeon/RAT_project/PHMG_IT/coworker
GATK=/home/yusanghyeon/RAT_project/miniforge3/envs/rnaseq/bin/gatk
BCFTOOLS=/home/yusanghyeon/RAT_project/miniforge3/envs/rnaseq/bin/bcftools

source /home/yusanghyeon/RAT_project/miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

SAMPLES=(SRR33625170 SRR33625169 SRR33625168 SRR33625167 SRR33625166 SRR33625165 \
         SRR33625154 SRR33625153 SRR33625152 SRR33625151 SRR33625150 SRR33625149 SRR33625148)
EARLY=(SRR33625170 SRR33625169 SRR33625168 SRR33625167 SRR33625166 SRR33625165)
LATE=(SRR33625154 SRR33625153 SRR33625152 SRR33625151 SRR33625150 SRR33625149 SRR33625148)
CHRS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 X)

mkdir -p $OUT/force_filtered $OUT/force_pass $OUT/force_raw $OUT/logs $OUT/diff_target $OUT/diff_ref

# === Step 1: Per-sample MergeMutectStats + Re-gather VCF with GATK GatherVcfs ===
echo "[$(date +%H:%M:%S)] Step 1: Gather VCFs + Merge stats per sample"
for S in "${SAMPLES[@]}"; do
  OUT_RAW=$OUT/force_raw/${S}.force.vcf.gz
  if [ -f $OUT_RAW.stats ] && [ -f $OUT_RAW ]; then continue; fi
  STATS_ARGS=""
  VCF_ARGS=""
  for CHR in "${CHRS[@]}"; do
    STATS_ARGS="$STATS_ARGS -stats $OUT/force_scatter/${S}_${CHR}.vcf.gz.stats"
    VCF_ARGS="$VCF_ARGS -I $OUT/force_scatter/${S}_${CHR}.vcf.gz"
  done
  $GATK GatherVcfs $VCF_ARGS -O $OUT_RAW > $OUT/logs/gather_${S}.log 2>&1
  $GATK MergeMutectStats $STATS_ARGS -O $OUT_RAW.stats >> $OUT/logs/gather_${S}.log 2>&1
  $BCFTOOLS index -t $OUT_RAW 2>/dev/null
  N=$($BCFTOOLS view -H $OUT_RAW 2>/dev/null | wc -l)
  echo "  $S: $N variants"
done

# === Step 2: FilterMutectCalls + PASS ===
echo "[$(date +%H:%M:%S)] Step 2: FilterMutectCalls + PASS (parallel 8)"
MAX_PAR=8
_wait_slot() {
  while [ $(jobs -r | wc -l) -ge $MAX_PAR ]; do sleep 1; done
}
for S in "${SAMPLES[@]}"; do
  FORCE_PASS=$OUT/force_pass/${S}.force_pass.vcf.gz
  if [ -f $FORCE_PASS ]; then continue; fi
  _wait_slot
  (
    $GATK FilterMutectCalls \
      -R $REF \
      -V $OUT/force_raw/${S}.force.vcf.gz \
      --min-median-base-quality 0 \
      --tmp-dir $OUT/force_filtered \
      -O $OUT/force_filtered/${S}.force_filtered.vcf.gz \
      > $OUT/logs/filter_${S}.log 2>&1 && \
    $BCFTOOLS view -f PASS $OUT/force_filtered/${S}.force_filtered.vcf.gz -Oz -o $FORCE_PASS 2>>$OUT/logs/filter_${S}.log && \
    $BCFTOOLS index -t $FORCE_PASS 2>>$OUT/logs/filter_${S}.log
  ) &
done
wait
echo "  PASS counts:"
for S in "${SAMPLES[@]}"; do
  N=$($BCFTOOLS view -H $OUT/force_pass/${S}.force_pass.vcf.gz 2>/dev/null | wc -l)
  printf "    %s  %s\n" "$S" "$N"
done

# === Step 3: coworker differential (late=target, early=reference) ===
echo "[$(date +%H:%M:%S)] Step 3: coworker differential analysis"
TARGET_VCFS=()
for S in "${LATE[@]}";  do TARGET_VCFS+=($OUT/force_pass/${S}.force_pass.vcf.gz); done
REF_VCFS=()
for S in "${EARLY[@]}"; do REF_VCFS+=($OUT/force_pass/${S}.force_pass.vcf.gz); done

DIFF_OUT=$OUT/differential_force
mkdir -p $DIFF_OUT

cd $COWORKER && python run_pipeline.py \
  --mode independent \
  --target   "${TARGET_VCFS[@]}" \
  --reference "${REF_VCFS[@]}" \
  --output   $DIFF_OUT \
  --min-coverage 3 \
  --min-qual 10 \
  --min-alt-ratio 0.1 \
  --max-ref-alt-ratio 0.05 \
  --min-recurrence 2 \
  --min-delta-alt-ratio 0.1 \
  > $OUT/logs/differential.log 2>&1

echo "[$(date +%H:%M:%S)] Done."
echo "  Differential SNPs: $(($(wc -l < $DIFF_OUT/differential_snps.csv) - 1))"
echo "  Differential INDELs: $(($(wc -l < $DIFF_OUT/differential_indels.csv) - 1))"
