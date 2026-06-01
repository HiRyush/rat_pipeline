#!/bin/bash
# 02b_reparse.sh
# Fix the parsing bug in Step 2: rebuild bucket_b1b2b3.tsv from raw pileup with correct column order.

set -euo pipefail

ROOT=/home/yusanghyeon/RAT_project/PHMG_IT
BASE=$ROOT/results/recall_decomposition
DATA=$BASE/data
RES=$BASE/results

source $ROOT/../miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

# Rebuild pileup_aggregated.tsv (from compressed raw_mpileup.tsv.gz)
echo "[1/3] Re-aggregating pileup (10 samples → total_dp + concatenated bases)..."
zcat $DATA/raw_mpileup.tsv.gz | awk 'BEGIN{OFS="\t"}{
    total_dp = 0
    bases_concat = ""
    for (s=0; s<10; s++) {
        col_dp   = 4 + 3*s
        col_base = 5 + 3*s
        dp = $(col_dp)
        if (dp == "") dp = 0
        total_dp += dp
        if (dp > 0) bases_concat = bases_concat $(col_base)
    }
    # Output: chrom, pos, ref_from_pileup, total_dp, bases
    print $1, $2, $3, total_dp, bases_concat
}' > $DATA/pileup_aggregated.tsv

echo "  Aggregated rows: $(wc -l < $DATA/pileup_aggregated.tsv)"

echo "[2/3] Joining with missed positions and reclassifying buckets (correctly this time)..."

# missed_3col.bed columns: chrom, start, end (=pos), ref, alt
# Key on chrom_pos (which is chrom_end)
awk 'BEGIN{OFS="\t"} {print $1"_"$3"\t"$4"\t"$5}' $DATA/missed_3col.bed | sort -k1,1 > /tmp/missed_alt.tsv
awk 'BEGIN{OFS="\t"} {print $1"_"$2"\t"$3"\t"$4"\t"$5}' $DATA/pileup_aggregated.tsv | sort -k1,1 > /tmp/pileup_kv.tsv

# After join: key | ref_missed | alt | ref_pileup | dp | bases  (6 fields)
join -t$'\t' -1 1 -2 1 /tmp/missed_alt.tsv /tmp/pileup_kv.tsv > /tmp/joined.tsv

echo "  Joined: $(wc -l < /tmp/joined.tsv)"

# Parse base column to count alt occurrences (handle pileup syntax)
awk 'BEGIN{OFS="\t"; print "chrom","pos","ref","alt","total_dp","total_alt_dp","bucket"}
{
    # Split key chrom_pos. Rat chroms: chr1..chr20, chrX, chrY, chrM — no underscores.
    split($1, k, "_")
    chrom = k[1]; pos = k[2]
    ref   = toupper($2)
    alt   = toupper($3)
    # $4 = ref_from_pileup, $5 = total_dp, $6 = bases
    dp    = $5 + 0
    bases = toupper($6)

    # Strip pileup syntax: ^X (start marker, skip 2 chars), $ (end marker, skip 1), +N/-N (indels)
    out_bases = ""
    i = 1
    L = length(bases)
    while (i <= L) {
        c = substr(bases, i, 1)
        if (c == "^") { i += 2; continue }
        if (c == "$") { i += 1; continue }
        if (c == "+" || c == "-") {
            j = i + 1
            num = ""
            while (j <= L && substr(bases, j, 1) ~ /[0-9]/) { num = num substr(bases, j, 1); j++ }
            i = j + (num+0)
            continue
        }
        out_bases = out_bases c
        i++
    }

    # Count alt allele in out_bases (single char comparison)
    alt_count = 0
    if (length(alt) == 1) {
        for (i = 1; i <= length(out_bases); i++) {
            if (substr(out_bases, i, 1) == alt) alt_count++
        }
    }

    bucket = (dp == 0) ? "B1" : (alt_count == 0 ? "B2" : "B3")
    print chrom, pos, ref, alt, dp, alt_count, bucket
}' /tmp/joined.tsv > $DATA/bucket_b1b2b3.tsv

n_joined=$(wc -l < /tmp/joined.tsv)
n_missed=$(wc -l < $DATA/missed.bed)
echo "  Output rows: $(($(wc -l < $DATA/bucket_b1b2b3.tsv) - 1))"

echo "[3/3] Bucket summary (corrected)..."
{
    echo "Bucket B1/B2/B3 classification summary (REPARSED)"
    echo "Generated: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo
    echo "Total missed: $n_missed"
    echo "Pileup-joined: $n_joined"
    n_unjoined=$((n_missed - n_joined))
    echo "Unjoined (no mpileup record) → B1: $n_unjoined"
    echo
    awk -F'\t' 'NR>1{c[$7]++; total++} END{
        for(b in c) printf "  %s (in joined): %d (%.1f%%)\n", b, c[b], c[b]*100/total
        printf "  TOTAL joined: %d\n", total
    }' $DATA/bucket_b1b2b3.tsv

    # Total B1 = unjoined + (joined with dp=0)
    n_b1_joined=$(awk -F'\t' 'NR>1 && $7=="B1"' $DATA/bucket_b1b2b3.tsv | wc -l)
    n_b1_total=$((n_unjoined + n_b1_joined))
    n_b2=$(awk -F'\t' 'NR>1 && $7=="B2"' $DATA/bucket_b1b2b3.tsv | wc -l)
    n_b3=$(awk -F'\t' 'NR>1 && $7=="B3"' $DATA/bucket_b1b2b3.tsv | wc -l)
    echo
    echo "=== Final bucket sizes (full missed set) ==="
    printf "  B1 (DP=0):  %d (%.1f%%)\n" $n_b1_total $(echo "scale=4; $n_b1_total*100/$n_missed" | bc)
    printf "  B2 (DP>0, alt=0):  %d (%.1f%%)\n" $n_b2 $(echo "scale=4; $n_b2*100/$n_missed" | bc)
    printf "  B3 (DP>0, alt>0):  %d (%.1f%%)\n" $n_b3 $(echo "scale=4; $n_b3*100/$n_missed" | bc)
    echo
    echo "=== Reachable recall (captured / (captured + B3)) ==="
    captured=$(wc -l < $DATA/captured.bed)
    reachable=$(echo "scale=4; $captured*100/($captured+$n_b3)" | bc)
    echo "  Captured: $captured"
    echo "  Captured + B3: $((captured + n_b3))"
    echo "  Reachable recall: ${reachable}%"
} | tee $RES/02_bucket_summary.txt

rm /tmp/missed_alt.tsv /tmp/pileup_kv.tsv /tmp/joined.tsv
gzip -f $DATA/pileup_aggregated.tsv

echo "Done."
