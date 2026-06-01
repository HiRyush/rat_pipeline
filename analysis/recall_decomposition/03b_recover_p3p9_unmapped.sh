#!/bin/bash
# 03b_recover_p3p9_unmapped.sh
# Recover B4 hits that were missed because P3 and P9 markdup BAMs had unmapped reads stripped.
# Source raw STAR BAMs (pilot/arm_C_P{3,9}) still contain the unmapped reads.
#
# Strategy:
#   1. Reuse existing alt k-mer set (b4_kmer_meta.tsv from Step 3)
#   2. Stream unmapped reads from P3 and P9 raw STAR BAMs
#   3. Add to existing b4_hits.tsv (cumulative)
#   4. Re-classify bucket_full.tsv with updated hits
#   5. Update summary

set -euo pipefail

ROOT=/home/yusanghyeon/RAT_project/PHMG_IT
BASE=$ROOT/results/recall_decomposition
DATA=$BASE/data
RES=$BASE/results
LOGS=$BASE/logs

source $ROOT/../miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

FLANK=15
KMER_LEN=$((2*FLANK + 1))
MIN_HITS=2

# Raw STAR BAMs (have unmapped reads that markdup BAM strips)
declare -A RAW_BAMS
RAW_BAMS[P3]=$ROOT/results/pilot/arm_C_P3/P3_Aligned.sortedByCoord.out.bam
RAW_BAMS[P9]=$ROOT/results/pilot/arm_C_P9/P9_Aligned.sortedByCoord.out.bam

for s in P3 P9; do
    if [[ ! -f "${RAW_BAMS[$s]}" ]]; then
        echo "ERROR: missing raw STAR BAM for $s" >&2; exit 1
    fi
done

if [[ ! -f $DATA/b1b2_flanks.tsv ]]; then
    echo "ERROR: b1b2_flanks.tsv missing — run Step 3 first." >&2; exit 1
fi

echo "[1/3] Rebuilding k-mer set from flanks + scanning P3+P9 raw STAR unmapped..."

python3 - <<PYEOF | tee $LOGS/p3p9_recovery.log
import sys, os, subprocess
from collections import defaultdict

FLANK = $FLANK
KMER = $KMER_LEN
DATA = "$DATA"
BAMS = {"P3": "${RAW_BAMS[P3]}", "P9": "${RAW_BAMS[P9]}"}

def revcomp(s):
    return s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]

# 1. Rebuild k-mer set from b1b2_flanks.tsv (same logic as Step 3)
kmer_to_keys = defaultdict(set)
position_info = {}
skipped = 0
with open(f"{DATA}/b1b2_flanks.tsv") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        name = parts[0].split("::")[0]
        seq = parts[1].upper()
        if len(seq) != KMER:
            skipped += 1; continue
        toks = name.split("_")
        if len(toks) < 5:
            skipped += 1; continue
        chrom = "_".join(toks[:-4])
        pos = toks[-4]; refb = toks[-3]; altb = toks[-2]; bucket = toks[-1]
        if seq[FLANK] != refb:
            skipped += 1; continue
        if altb not in "ACGT" or len(altb) != 1:
            skipped += 1; continue
        kmer_fwd = seq[:FLANK] + altb + seq[FLANK+1:]
        if "N" in kmer_fwd:
            skipped += 1; continue
        kmer_rc = revcomp(kmer_fwd)
        key = f"{chrom}:{pos}:{refb}:{altb}"
        position_info[key] = (chrom, pos, refb, altb, bucket)
        kmer_to_keys[kmer_fwd].add(key)
        kmer_to_keys[kmer_rc].add(key)
print(f"  Rebuilt {len(position_info)} positions, {len(kmer_to_keys)} unique k-mers (skipped {skipped})", file=sys.stderr)

# 2. Scan P3, P9 raw STAR BAMs
new_hits = defaultdict(int)
total_reads = 0
for sample, bam in BAMS.items():
    print(f"  Scanning {sample} raw STAR BAM unmapped...", file=sys.stderr)
    proc = subprocess.Popen(
        ["samtools", "view", "-f", "4", "-F", "256", bam],
        stdout=subprocess.PIPE, text=True
    )
    n = 0
    for line in proc.stdout:
        cols = line.split("\t", 11)
        if len(cols) < 10: continue
        seq = cols[9]
        if len(seq) < KMER: continue
        n += 1
        for i in range(len(seq) - KMER + 1):
            kmer = seq[i:i+KMER]
            if kmer in kmer_to_keys:
                for key in kmer_to_keys[kmer]:
                    new_hits[key] += 1
    proc.wait()
    total_reads += n
    print(f"    {sample}: {n} unmapped reads scanned, cumulative hits added: {sum(new_hits.values())}", file=sys.stderr)

print(f"  Total new reads scanned: {total_reads}", file=sys.stderr)
print(f"  Positions with NEW hits (from P3/P9): {len(new_hits)}", file=sys.stderr)

# 3. Merge with existing b4_hits.tsv: existing hits + new hits
existing = {}  # key -> (chrom, pos, ref, alt, orig_bucket, hits)
with open(f"{DATA}/b4_hits.tsv") as f:
    header = next(f)
    for line in f:
        chrom, pos, ref, alt, bucket, hits, is_b4 = line.rstrip("\n").split("\t")
        key = f"{chrom}:{pos}:{ref}:{alt}"
        existing[key] = (chrom, pos, ref, alt, bucket, int(hits))

# Backup original
import shutil
shutil.copy(f"{DATA}/b4_hits.tsv", f"{DATA}/b4_hits.preP3P9.tsv")

# Write updated table with combined hits
n_newly_b4 = 0
n_already_b4 = 0
with open(f"{DATA}/b4_hits.tsv", "w") as f:
    f.write(header)
    MIN_HITS = $MIN_HITS
    for key, (chrom, pos, ref, alt, bucket, hits_old) in existing.items():
        hits_new = hits_old + new_hits.get(key, 0)
        was_b4 = hits_old >= MIN_HITS
        is_b4 = hits_new >= MIN_HITS
        if (not was_b4) and is_b4:
            n_newly_b4 += 1
        elif was_b4:
            n_already_b4 += 1
        f.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{bucket}\t{hits_new}\t{int(is_b4)}\n")
    # New positions that weren't in existing (shouldn't happen, but check)
    for key, n in new_hits.items():
        if key not in existing and n >= MIN_HITS:
            chrom, pos, ref, alt, bucket = position_info[key]
            f.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{bucket}\t{n}\t1\n")
            n_newly_b4 += 1

print(f"  Newly classified B4 (was B1/B2, now B4): {n_newly_b4}", file=sys.stderr)
print(f"  Already B4 (kept): {n_already_b4}", file=sys.stderr)
PYEOF

echo "[2/3] Rebuilding bucket_full.tsv with updated B4..."
python3 - <<PYEOF
import csv
b4 = {}
with open("$DATA/b4_hits.tsv") as f:
    next(f)
    for row in csv.reader(f, delimiter="\t"):
        chrom, pos, ref, alt, bucket, hits, is_b4 = row
        b4[(chrom, pos, ref, alt)] = (int(hits), bucket)

import shutil
shutil.copy("$DATA/bucket_full.tsv", "$DATA/bucket_full.preP3P9.tsv")

with open("$DATA/bucket_full.tsv", "w") as out, open("$DATA/bucket_b1b2b3.tsv") as f:
    out.write("chrom\tpos\tref\talt\ttotal_dp\ttotal_alt_dp\tbucket_b1b2b3\tunmapped_kmer_hits\tbucket_final\n")
    next(f)
    for line in f:
        parts = line.rstrip("\n").split("\t")
        chrom, pos, ref, alt, dp, alt_dp, b = parts
        key = (chrom, pos, ref, alt)
        if key in b4:
            hits, _ = b4[key]
        else:
            hits = 0
        if b in ("B1","B2") and hits >= $MIN_HITS:
            final = "B4"
        else:
            final = b
        out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{dp}\t{alt_dp}\t{b}\t{hits}\t{final}\n")

from collections import Counter
cnt = Counter()
with open("$DATA/bucket_full.tsv") as f:
    next(f)
    for line in f:
        cnt[line.split("\t")[8].strip()] += 1
print("New bucket counts:")
for b in ("B1","B2","B3","B4"):
    print(f"  {b}: {cnt.get(b, 0)}")
PYEOF

echo "[3/3] Updated summary..."
{
    echo "B4 recovery from P3/P9 raw STAR BAMs"
    echo "Generated: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo
    n_missed=$(wc -l < $DATA/missed.bed)
    n_captured=$(wc -l < $DATA/captured.bed)
    echo "Captured: $n_captured  Missed: $n_missed"
    echo
    n_b1=$(awk -F'\t' 'NR>1 && $9=="B1"' $DATA/bucket_full.tsv | wc -l)
    n_b2=$(awk -F'\t' 'NR>1 && $9=="B2"' $DATA/bucket_full.tsv | wc -l)
    n_b3=$(awk -F'\t' 'NR>1 && $9=="B3"' $DATA/bucket_full.tsv | wc -l)
    n_b4=$(awk -F'\t' 'NR>1 && $9=="B4"' $DATA/bucket_full.tsv | wc -l)
    n_total=$((n_b1+n_b2+n_b3+n_b4))
    n_unjoined=$((n_missed - n_total))
    [[ $n_unjoined -gt 0 ]] && n_b1=$((n_b1 + n_unjoined))

    echo "Bucket sizes (updated with P3/P9 unmapped):"
    printf "  B1 (DP=0): %d (%.1f%%)\n" $n_b1 $(echo "scale=4; $n_b1*100/$n_missed" | bc)
    printf "  B2 (DP>0, alt=0): %d (%.1f%%)\n" $n_b2 $(echo "scale=4; $n_b2*100/$n_missed" | bc)
    printf "  B3 (DP>0, alt>0): %d (%.1f%%)\n" $n_b3 $(echo "scale=4; $n_b3*100/$n_missed" | bc)
    printf "  B4 (unmapped k-mer hit): %d (%.1f%%)\n" $n_b4 $(echo "scale=4; $n_b4*100/$n_missed" | bc)
} | tee $RES/03b_b4_recovery_summary.txt

echo "Done. Run Step 4 to regenerate recall report with updated B4."
