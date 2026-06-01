#!/bin/bash
# build_dna_truth_treatment_specific.sh
# Build DNA treatment-specific SNP set as denominator for method recall.
# Input: 15 per-sample DNA SNP VCFs (Isaac variant caller output)
# Output:
#   dna_sites_long.tsv.gz  — long form (chr, pos, ref, alt, sample, group)
#   dna_truth_counts.tsv.gz — per-variant (chr, pos, ref, alt, ctrl_n, treat_n)
#   dna_truth_summary.txt   — set sizes under different (ctrl, treat) thresholds
#
# Usage:
#   bash build_dna_truth_treatment_specific.sh [CHROM_FILTER]
# If CHROM_FILTER given (e.g. "chr1"), restricts to that chromosome (for testing).

set -euo pipefail

GT_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/ground_truth"
OUT_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/dna_truth_coverage"
CHROM_FILTER="${1:-}"

mkdir -p "$OUT_DIR"

source /home/yusanghyeon/RAT_project/miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

LONG_TSV="$OUT_DIR/dna_sites_long.tsv"
COUNTS_TSV="$OUT_DIR/dna_truth_counts.tsv"
SUMMARY="$OUT_DIR/dna_truth_summary.txt"

if [[ -n "$CHROM_FILTER" ]]; then
    LONG_TSV="${LONG_TSV%.tsv}.${CHROM_FILTER}.tsv"
    COUNTS_TSV="${COUNTS_TSV%.tsv}.${CHROM_FILTER}.tsv"
    SUMMARY="${SUMMARY%.txt}.${CHROM_FILTER}.txt"
fi

CTRL_SAMPLES=(C1 C2 C3 C4 C5)
TREAT_SAMPLES=(P1 P2 P3 P4 P5 P6 P7 P8 P9 P10)

echo "[1/3] Extracting ALT-bearing sites from 15 DNA VCFs..."
> "$LONG_TSV"

extract_one() {
    local sample="$1" group="$2"
    local vcf="$GT_DIR/${sample}_dna_snps.vcf.gz"
    local region_arg=()
    if [[ -n "$CHROM_FILTER" ]]; then
        region_arg=(-r "$CHROM_FILTER")
    fi
    # Only PASS or LowGQX? Use all variants (PASS-only would discard recoverable),
    # but require non-ref GT. Split multi-allelic via norm.
    bcftools norm -m- -Ou "${region_arg[@]}" "$vcf" 2>/dev/null \
      | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' \
      | awk -v s="$sample" -v g="$group" 'BEGIN{OFS="\t"}{
            gt=$5
            # any allele == "1" → non-ref
            if(gt ~ /1/) print $1,$2,$3,$4,s,g
        }'
}

for s in "${CTRL_SAMPLES[@]}"; do
    echo "  $s (Control)"
    extract_one "$s" "C" >> "$LONG_TSV"
done
for s in "${TREAT_SAMPLES[@]}"; do
    echo "  $s (Treatment)"
    extract_one "$s" "T" >> "$LONG_TSV"
done

echo "  Long TSV rows: $(wc -l < "$LONG_TSV")"

echo "[2/3] Aggregating to per-variant counts..."
# Sort by chr,pos,ref,alt then aggregate
sort -k1,1 -k2,2n -k3,3 -k4,4 -T "$OUT_DIR" "$LONG_TSV" \
  | awk -F'\t' 'BEGIN{OFS="\t"; print "chrom","pos","ref","alt","ctrl_n","treat_n"}
    {
        key=$1"\t"$2"\t"$3"\t"$4
        if(key != prev && prev != ""){
            print pchrom, ppos, pref, palt, c, t
            c=0; t=0
        }
        if($6=="C") c++
        else if($6=="T") t++
        prev=key; pchrom=$1; ppos=$2; pref=$3; palt=$4
    }
    END{ if(prev!="") print pchrom, ppos, pref, palt, c, t }' \
  > "$COUNTS_TSV"

echo "  Counts TSV rows: $(($(wc -l < "$COUNTS_TSV") - 1))"

echo "[3/3] Summary of treatment-specific set sizes..."
{
    echo "# DNA treatment-specific SNP set sizes (PHMG_IT, ${CHROM_FILTER:-all chromosomes})"
    echo "# Generated: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "# Source: per-sample DNA SNP VCFs from $GT_DIR (Isaac caller, multi-allele split by bcftools norm)"
    echo "# Denominator definitions tested:"
    echo
    printf "%-40s %12s\n" "Definition" "Count"
    printf "%-40s %12s\n" "----------------------------------------" "------------"

    total=$(awk 'NR>1' "$COUNTS_TSV" | wc -l)
    printf "%-40s %12d\n" "Total unique (chr,pos,ref,alt)" "$total"

    for thresh in 1 2 3 5; do
        n=$(awk -F'\t' -v t="$thresh" 'NR>1 && $5==0 && $6>=t' "$COUNTS_TSV" | wc -l)
        printf "%-40s %12d\n" "ctrl_n=0 & treat_n>=${thresh}" "$n"
    done

    # Also: shared (would-be germline)
    shared=$(awk -F'\t' 'NR>1 && $5>0 && $6>0' "$COUNTS_TSV" | wc -l)
    printf "%-40s %12d\n" "ctrl_n>0 & treat_n>0 (shared/germline)" "$shared"
    ctrl_only=$(awk -F'\t' 'NR>1 && $5>0 && $6==0' "$COUNTS_TSV" | wc -l)
    printf "%-40s %12d\n" "ctrl_n>0 & treat_n=0 (ctrl-only)" "$ctrl_only"
} | tee "$SUMMARY"

# Compress outputs
gzip -f "$LONG_TSV"
gzip -f "$COUNTS_TSV"

echo "Done. Outputs in $OUT_DIR"
