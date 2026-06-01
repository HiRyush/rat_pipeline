#!/bin/bash
# 01_build_inputs.sh
# Build input sets for recall decomposition:
#   - DNA truth BED (treat>=2 & ctrl=0, Option A definition)
#   - RNA captured BED (947 TRUE_SOMATIC positions)
#   - Missed BED = DNA truth - RNA captured
#
# Output: data/{dna_truth.bed, rna_captured.bed, missed.bed, summary.txt}

set -euo pipefail

ROOT=/home/yusanghyeon/RAT_project/PHMG_IT
BASE=$ROOT/results/recall_decomposition
DATA=$BASE/data

source $ROOT/../miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

DNA_COUNTS=$ROOT/results/dna_truth_coverage/dna_truth_counts_OptA.tsv.gz
RNA_TSV=$ROOT/results/dna_truth_coverage/rna_captured.clean.tsv   # CRLF-stripped version

# Threshold for treatment recurrence in DNA truth
TREAT_MIN="${1:-2}"

echo "[1/4] Building DNA truth BED (treat>=${TREAT_MIN} & ctrl=0, Option A)..."
zcat $DNA_COUNTS \
  | awk -F'\t' -v t=$TREAT_MIN 'NR>1 && $5==0 && $6>=t {print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' \
  | sort -k1,1 -k2,2n -u > $DATA/dna_truth.bed
echo "  DNA truth positions: $(wc -l < $DATA/dna_truth.bed)"

echo "[2/4] Building RNA captured BED (TRUE_SOMATIC only)..."
awk -F'\t' 'NR>1 && $11=="TRUE_SOMATIC" {print $1"\t"$2-1"\t"$2"\t"$3"\t"$4}' $RNA_TSV \
  | sort -k1,1 -k2,2n -u > $DATA/rna_captured.bed
echo "  RNA captured positions: $(wc -l < $DATA/rna_captured.bed)"

echo "[3/4] Building missed set (DNA truth minus captured, position-level)..."
# Use chr_pos keys for set difference (allele-agnostic at position level for bucket analysis,
# since RNA may have called a different ALT than DNA — those count as "captured" for recall purposes
# but we want to be strict on EXACT match here for fairness).
# Strict allele-match version:
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' $DATA/dna_truth.bed | sort > /tmp/dna_keys.tsv
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' $DATA/rna_captured.bed | sort > /tmp/rna_keys.tsv
comm -23 /tmp/dna_keys.tsv /tmp/rna_keys.tsv > $DATA/missed.bed
echo "  Missed positions: $(wc -l < $DATA/missed.bed)"

# Captured (DNA truth ∩ RNA captured)
comm -12 /tmp/dna_keys.tsv /tmp/rna_keys.tsv > $DATA/captured.bed
echo "  Captured (DNA truth ∩ RNA captured): $(wc -l < $DATA/captured.bed)"

echo "[4/4] Building auxiliary universes for downstream reporting..."
# MuTect2 evaluated BED
cp $ROOT/results/dna_truth_coverage/mutect2_evaluated.bed $DATA/mutect2_evaluated.bed
# Transcript BED
cp $ROOT/results/dna_truth_coverage/transcript_regions.bed $DATA/transcript_regions.bed
# Gene BED with names
cp $ROOT/results/dna_truth_coverage/genes.bed $DATA/genes.bed

# Summary
{
    echo "Recall decomposition input summary"
    echo "Generated: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "Treatment recurrence threshold: >= $TREAT_MIN"
    echo
    echo "DNA truth positions:  $(wc -l < $DATA/dna_truth.bed)"
    echo "RNA captured (TRUE_SOMATIC): $(wc -l < $DATA/rna_captured.bed)"
    echo "Captured (intersection):  $(wc -l < $DATA/captured.bed)"
    echo "Missed (DNA - RNA):  $(wc -l < $DATA/missed.bed)"
    echo
    echo "Auxiliary BEDs:"
    echo "  mutect2_evaluated.bed: $(wc -l < $DATA/mutect2_evaluated.bed)"
    echo "  transcript_regions.bed: $(wc -l < $DATA/transcript_regions.bed)"
    echo "  genes.bed: $(wc -l < $DATA/genes.bed) lines"
} | tee $BASE/results/01_input_summary.txt

rm /tmp/dna_keys.tsv /tmp/rna_keys.tsv

echo "Done. Outputs in $DATA"
