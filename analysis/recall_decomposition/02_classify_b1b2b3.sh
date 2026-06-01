#!/bin/bash
# 02_classify_b1b2b3.sh
# Classify missed DNA-truth positions into B1/B2/B3 via mpileup on treatment BAMs.
#
# B1: total_DP = 0 across all 10 treatment samples (not expressed)
# B2: total_DP > 0 but total_alt_DP = 0 (alt reads dropped/misaligned)
# B3: total_DP > 0 and total_alt_DP > 0 (true method-level miss)
#
# Output:
#   data/missed_pileup.tsv : per-position chrom pos ref alt total_dp total_alt_dp
#   data/bucket_b1b2b3.tsv : per-position with bucket label
#   results/02_bucket_summary.txt : counts

set -euo pipefail

ROOT=/home/yusanghyeon/RAT_project/PHMG_IT
BASE=$ROOT/results/recall_decomposition
DATA=$BASE/data
RES=$BASE/results
LOGS=$BASE/logs

source $ROOT/../miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

REF=$ROOT/reference/rn7.fa
TREAT_BAMS=(P1 P2 P3 P4 P5 P6 P7 P8 P9 P10)

# Build proper 3-col BED for mpileup. missed.bed already has BED-format coords (start=pos-1, end=pos)
# from Step 1's comm output, so we just pass through.
echo "[1/4] Preparing position list..."
cp $DATA/missed.bed $DATA/missed_3col.bed
echo "  Missed positions: $(wc -l < $DATA/missed_3col.bed)"

# Build BAM list arg
BAM_ARGS=""
for s in "${TREAT_BAMS[@]}"; do
    BAM_ARGS="$BAM_ARGS $ROOT/results/mutect2/markdup/${s}.dedup.bam"
done

echo "[2/4] Running multi-sample mpileup (10 treatment BAMs)..."
echo "  Output: $DATA/raw_mpileup.tsv"
# -A: count anomalous; -B: no BAQ; -q 0 -Q 0: include all reads
# -d 100000: deep pileup
samtools mpileup -f $REF -A -B -q 0 -Q 0 -d 100000 \
    -l $DATA/missed_3col.bed \
    $BAM_ARGS \
    2> $LOGS/mpileup.stderr \
    > $DATA/raw_mpileup.tsv
echo "  mpileup output lines: $(wc -l < $DATA/raw_mpileup.tsv)"

echo "[3/4] Parsing pileup to per-position DP/alt_DP across 10 samples..."
# Each row: chrom pos ref [DP base qual] x 10
# Sum total_DP over all 10. For total_alt_DP we need the alt allele per position.
# Join with missed_3col to get alt allele.

# First aggregate per-position from pileup
awk 'BEGIN{OFS="\t"}{
    total_dp = 0
    bases_concat = ""
    # 10 samples: columns 4,5,6 = sample1 dp,base,qual; 7,8,9 = sample2; ...
    for (s=0; s<10; s++) {
        col_dp   = 4 + 3*s
        col_base = 5 + 3*s
        dp = $(col_dp)
        if (dp == "") dp = 0
        total_dp += dp
        if (dp > 0) bases_concat = bases_concat $(col_base)
    }
    print $1, $2, $3, total_dp, bases_concat
}' $DATA/raw_mpileup.tsv > $DATA/pileup_aggregated.tsv

echo "  Aggregated rows: $(wc -l < $DATA/pileup_aggregated.tsv)"

# Join with missed positions to add alt allele
# Key on chrom_pos
awk 'BEGIN{OFS="\t"} {print $1"_"$3, $4, $5}' $DATA/missed_3col.bed | sort > /tmp/missed_alt.tsv
awk 'BEGIN{OFS="\t"} {print $1"_"$2, $3, $4, $5}' $DATA/pileup_aggregated.tsv | sort > /tmp/pileup_kv.tsv

# join: key, alt, total_dp, ref_from_pileup, bases
join -t$'\t' -1 1 -2 1 /tmp/missed_alt.tsv /tmp/pileup_kv.tsv > /tmp/joined.tsv

# Now parse bases string to count occurrences of alt allele.
# After join: key | ref_missed | alt | ref_pileup | total_dp | bases
# (missed_alt: $1=key, $2=ref, $3=alt;  pileup_kv: $1=key, $2=ref_pileup, $3=dp, $4=bases)
awk 'BEGIN{OFS="\t"; print "chrom","pos","ref","alt","total_dp","total_alt_dp","bucket"}
{
    split($1, k, "_")
    chrom=k[1]; pos=k[2]
    ref = toupper($2)        # ref from missed_alt
    alt = toupper($3)        # alt from missed_alt
    ref_pileup = toupper($4) # ref from pileup (sanity)
    dp  = $5 + 0
    bases = toupper($6)

    # Strip indel markers: +n[ACGTNacgtn]{n}
    # Simplification: remove ^X (start) and $ end, then count alt base occurrences
    # outside of indel insertion strings.
    # Pragmatic regex-free approach:
    out_bases = ""
    i = 1
    L = length(bases)
    while (i <= L) {
        c = substr(bases, i, 1)
        if (c == "^") { i += 2; continue }  # skip ^X
        if (c == "$") { i += 1; continue }  # skip $
        if (c == "+" || c == "-") {
            # parse digit run
            j = i + 1
            num = ""
            while (j <= L && substr(bases, j, 1) ~ /[0-9]/) {
                num = num substr(bases, j, 1); j++
            }
            i = j + (num+0)  # skip indel bases
            continue
        }
        out_bases = out_bases c
        i++
    }
    # count alt occurrences (single-character SNV alt)
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
echo "  Joined positions: $n_joined / $n_missed missed"

echo "[4/4] Bucket summary..."
{
    echo "Bucket B1/B2/B3 classification summary"
    echo "Generated: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "Source: mpileup of 10 treatment BAMs at missed DNA truth positions"
    echo
    echo "Total missed: $n_missed"
    echo "Pileup-joined: $n_joined  (unjoined = no record from mpileup, treated as B1 separately)"
    echo
    awk -F'\t' 'NR>1{c[$7]++; total++} END{
        for(b in c) printf "  %s: %d (%.1f%%)\n", b, c[b], c[b]*100/total
        printf "  TOTAL: %d\n", total
    }' $DATA/bucket_b1b2b3.tsv

    # Unjoined (no record at all in mpileup) — treat as B1
    n_unjoined=$((n_missed - n_joined))
    echo
    echo "Unjoined (no mpileup record at all) → reclassify as B1: $n_unjoined"
} | tee $RES/02_bucket_summary.txt

rm /tmp/missed_alt.tsv /tmp/pileup_kv.tsv /tmp/joined.tsv
gzip -f $DATA/raw_mpileup.tsv
gzip -f $DATA/pileup_aggregated.tsv

echo "Done. Bucket table: $DATA/bucket_b1b2b3.tsv"
