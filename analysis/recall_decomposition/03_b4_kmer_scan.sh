#!/bin/bash
# 03_b4_kmer_scan.sh
# Detect B4 positions: alt-allele reads that failed STAR alignment (in unmapped pool).
#
# Strategy (no pysam/kmc dependency — pure samtools+bedtools+python set ops):
#   1. Take B1+B2 missed positions
#   2. Build flank BED (pos-15 .. pos+15) — 31bp window centered on variant
#   3. bedtools getfasta → reference flanking sequences
#   4. python: swap center base ref→alt, generate fwd + revcomp 31-mer per position
#   5. samtools view -f 4 unmapped reads → stream sequences
#   6. For each read, slide 31-mer window, check hash-set membership
#   7. Aggregate per-position hit counts → call B4 if hits >= MIN_HITS

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

KMER_LEN=31
FLANK=15           # KMER_LEN/2 floor
MIN_HITS=2

if [[ ! -f $DATA/bucket_b1b2b3.tsv ]]; then
    echo "ERROR: bucket_b1b2b3.tsv missing." >&2; exit 1
fi

echo "[1/5] Selecting B1+B2 positions..."
awk -F'\t' 'NR>1 && ($7=="B1" || $7=="B2"){print $1"\t"$2"\t"$3"\t"$4"\t"$7}' $DATA/bucket_b1b2b3.tsv \
    > $DATA/b1b2_positions.tsv
n_b12=$(wc -l < $DATA/b1b2_positions.tsv)
echo "  B1+B2 positions: $n_b12"

# Also include the unjoined (no mpileup record) — these are also missed positions that need B4 check
# Get DNA truth positions that did NOT appear in bucket_b1b2b3.tsv
awk -F'\t' 'NR>1{print $1"\t"$2}' $DATA/bucket_b1b2b3.tsv | sort -u > /tmp/joined_keys.tsv
awk '{print $1"\t"$3}' $DATA/missed.bed | sort -u > /tmp/missed_keys.tsv
comm -23 /tmp/missed_keys.tsv /tmp/joined_keys.tsv > /tmp/unjoined_keys.tsv
n_unjoin=$(wc -l < /tmp/unjoined_keys.tsv)
echo "  Unjoined (treated as B1) added: $n_unjoin"

# Get ref/alt for unjoined positions from missed.bed
awk 'NR==FNR{key=$1"\t"$2; keys[key]=1; next} {key=$1"\t"$3; if(key in keys) print $1"\t"$3"\t"$4"\t"$5"\tB1"}' \
    /tmp/unjoined_keys.tsv $DATA/missed.bed >> $DATA/b1b2_positions.tsv

n_total=$(wc -l < $DATA/b1b2_positions.tsv)
echo "  Total to k-mer scan: $n_total"

echo "[2/5] Building flank BED (pos-${FLANK} .. pos+${FLANK})..."
awk -v f=$FLANK 'BEGIN{OFS="\t"}{
    start = $2 - 1 - f
    end   = $2 + f
    if (start < 0) start = 0
    print $1, start, end, $1"_"$2"_"$3"_"$4"_"$5
}' $DATA/b1b2_positions.tsv > $DATA/b1b2_flank.bed
echo "  Flank BED: $(wc -l < $DATA/b1b2_flank.bed) entries"

echo "[3/5] Extracting flanking reference sequences..."
bedtools getfasta -fi $REF -bed $DATA/b1b2_flank.bed -name -tab \
    > $DATA/b1b2_flanks.tsv 2>$LOGS/getfasta.stderr
echo "  Extracted: $(wc -l < $DATA/b1b2_flanks.tsv)"

echo "[4/5] Building alt 31-mer set + scanning unmapped reads..."
KMER_PY=$DATA/_kmer_scan.py
cat > $KMER_PY <<'PYEOF'
import sys, os, subprocess
from collections import defaultdict

FLANK = int(os.environ['FLANK'])
KMER = 2*FLANK + 1
ROOT = os.environ['ROOT']
DATA = os.environ['DATA']
LOGS = os.environ['LOGS']
TREAT = os.environ['TREAT'].split()

# 1. Build alt 31-mer set + mapping kmer → position_key
def revcomp(s):
    return s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]

kmer_to_keys = defaultdict(set)  # kmer -> set of position keys
position_info = {}                # key -> (chrom, pos, ref, alt, orig_bucket)

skipped = 0; built = 0
with open(f"{DATA}/b1b2_flanks.tsv") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        # name col format: chrom_pos_ref_alt_bucket  followed by ::chrom:start-end (bedtools -name appends)
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
        built += 1
print(f"  Built {built} positions, {len(kmer_to_keys)} unique kmers, skipped {skipped}", file=sys.stderr)

# 2. Stream unmapped reads from each treat BAM, scan for k-mer hits
hits = defaultdict(int)  # position_key -> total read-hit count
total_reads = 0

samtools = "samtools"
for sample in TREAT:
    bam = f"{ROOT}/results/mutect2/markdup/{sample}.dedup.bam"
    print(f"  Scanning unmapped reads from {sample}...", file=sys.stderr)
    proc = subprocess.Popen(
        [samtools, "view", "-f", "4", "-F", "256", bam],
        stdout=subprocess.PIPE, text=True
    )
    n_reads_sample = 0
    for line in proc.stdout:
        cols = line.split("\t", 11)
        if len(cols) < 10: continue
        seq = cols[9]
        if len(seq) < KMER: continue
        n_reads_sample += 1
        # Slide window
        for i in range(len(seq) - KMER + 1):
            kmer = seq[i:i+KMER]
            if kmer in kmer_to_keys:
                for key in kmer_to_keys[kmer]:
                    hits[key] += 1
    proc.wait()
    total_reads += n_reads_sample
    print(f"    {sample}: {n_reads_sample} reads scanned, cumulative hits to date: {sum(hits.values())}", file=sys.stderr)

print(f"  Total unmapped reads scanned: {total_reads}", file=sys.stderr)
print(f"  Total positions with hits: {len(hits)}", file=sys.stderr)

# 3. Emit per-position hit table + final bucket reclassification
MIN_HITS = int(os.environ.get('MIN_HITS', '2'))
out_hits = open(f"{DATA}/b4_hits.tsv", "w")
out_hits.write("chrom\tpos\tref\talt\torig_bucket\tunmapped_kmer_hits\tcalled_b4\n")
for key, info in position_info.items():
    chrom, pos, refb, altb, bucket = info
    n = hits.get(key, 0)
    is_b4 = (n >= MIN_HITS)
    out_hits.write(f"{chrom}\t{pos}\t{refb}\t{altb}\t{bucket}\t{n}\t{int(is_b4)}\n")
out_hits.close()
print("  Wrote b4_hits.tsv", file=sys.stderr)
PYEOF

KMER_LEN=$KMER_LEN FLANK=$FLANK ROOT=$ROOT DATA=$DATA LOGS=$LOGS \
    TREAT="${TREAT_BAMS[*]}" MIN_HITS=$MIN_HITS \
    python3 $KMER_PY 2>&1 | tee $LOGS/kmer_scan.log

echo "[5/5] Building final bucket table + summary..."
# Merge: bucket_b1b2b3.tsv + b4_hits.tsv → bucket_full.tsv
python3 - <<PYEOF
import csv
b4 = {}  # (chrom,pos,ref,alt) -> hits
with open("$DATA/b4_hits.tsv") as f:
    next(f)
    for row in csv.reader(f, delimiter="\t"):
        chrom, pos, ref, alt, bucket, hits, is_b4 = row
        b4[(chrom, pos, ref, alt)] = (int(hits), bucket)

# Read bucket_b1b2b3 and add b4 info
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

# Counts
from collections import Counter
cnt = Counter()
with open("$DATA/bucket_full.tsv") as f:
    next(f)
    for line in f:
        cnt[line.split("\t")[8].strip()] += 1
print("Bucket final counts:")
for b, n in sorted(cnt.items()):
    print(f"  {b}: {n}")
PYEOF

{
    echo "B4 detection + final bucket summary"
    echo "Generated: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "K-mer length: $KMER_LEN, MIN_HITS: $MIN_HITS"
    echo
    n_missed=$(wc -l < $DATA/missed.bed)
    n_captured=$(wc -l < $DATA/captured.bed)
    echo "Total missed: $n_missed  Captured: $n_captured"
    echo

    n_b1=$(awk -F'\t' 'NR>1 && $9=="B1"' $DATA/bucket_full.tsv | wc -l)
    n_b2=$(awk -F'\t' 'NR>1 && $9=="B2"' $DATA/bucket_full.tsv | wc -l)
    n_b3=$(awk -F'\t' 'NR>1 && $9=="B3"' $DATA/bucket_full.tsv | wc -l)
    n_b4=$(awk -F'\t' 'NR>1 && $9=="B4"' $DATA/bucket_full.tsv | wc -l)
    n_total_joined=$((n_b1 + n_b2 + n_b3 + n_b4))
    n_unjoined=$((n_missed - n_total_joined))

    echo "From joined positions:"
    printf "  B1: %d (%.1f%%)\n" $n_b1 $(echo "scale=4; $n_b1*100/$n_missed" | bc)
    printf "  B2: %d (%.1f%%)\n" $n_b2 $(echo "scale=4; $n_b2*100/$n_missed" | bc)
    printf "  B3: %d (%.1f%%)\n" $n_b3 $(echo "scale=4; $n_b3*100/$n_missed" | bc)
    printf "  B4: %d (%.1f%%)\n" $n_b4 $(echo "scale=4; $n_b4*100/$n_missed" | bc)
    if [[ $n_unjoined -gt 0 ]]; then
        echo "  Unjoined (no mpileup record): $n_unjoined (counted as B1)"
        n_b1=$((n_b1 + n_unjoined))
    fi
    echo
    echo "=== Recall under different denominator definitions ==="
    naive=$(echo "scale=4; $n_captured*100/$n_missed" | bc)
    naive=$(echo "scale=4; $n_captured*100/($n_captured+$n_missed)" | bc)
    reach=$(echo "scale=4; $n_captured*100/($n_captured+$n_b3)" | bc)
    reach_b2b3=$(echo "scale=4; $n_captured*100/($n_captured+$n_b2+$n_b3)" | bc)
    no_b4=$(echo "scale=4; $n_captured*100/($n_captured+$n_b1+$n_b2+$n_b3)" | bc)
    echo "  Naive (cap / total truth): $naive%"
    echo "  Without B1+B4 (cap / (cap+B2+B3)): $reach_b2b3%"
    echo "  Reachable (cap / (cap+B3)): $reach%"
} | tee $RES/03_b4_summary.txt

rm -f /tmp/joined_keys.tsv /tmp/missed_keys.tsv /tmp/unjoined_keys.tsv
echo "Done. Final bucket: $DATA/bucket_full.tsv"
