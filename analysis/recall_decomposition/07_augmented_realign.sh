#!/bin/bash
# 07_augmented_realign.sh
# Augmented mini-reference realignment for B4 verification + expansion.
#
# Idea: Build a small reference containing alt-allele windows at each B1/B2 position.
# Re-align unmapped reads (from all 10 treatment BAMs incl P3/P9 raw) to this mini-ref.
# Reads aligning uniquely to alt windows → evidence that the read carries the alt allele
# and was rejected by main genome alignment due to reference bias.
#
# Outputs:
#   data/mini_ref.fa         alt-allele windows
#   data/unmapped_all.fq.gz  pooled unmapped reads
#   data/mini_ref.sam.gz     BWA alignments
#   data/augmented_hits.tsv  per-position read counts from realignment
#   results/07_augmented_summary.txt

set -euo pipefail

ROOT=/home/yusanghyeon/RAT_project/PHMG_IT
BASE=$ROOT/results/recall_decomposition
DATA=$BASE/data
RES=$BASE/results
LOGS=$BASE/logs

source $ROOT/../miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

REF=$ROOT/reference/rn7.fa
HALF_WIN=100  # window = 201bp centered on alt
CENTER=$((HALF_WIN + 1))  # 1-based position of alt in window

# Treatment BAMs: prefer raw STAR BAM for P3/P9 (has unmapped), markdup for others
declare -A BAMS
for s in P1 P2 P4 P5 P6 P7 P8 P10; do
    BAMS[$s]=$ROOT/results/mutect2/markdup/${s}.dedup.bam
done
BAMS[P3]=$ROOT/results/pilot/arm_C_P3/P3_Aligned.sortedByCoord.out.bam
BAMS[P9]=$ROOT/results/pilot/arm_C_P9/P9_Aligned.sortedByCoord.out.bam

# ---- 1. Build mini-ref FASTA ----
echo "[1/6] Building augmented mini-reference (B1+B2 positions, ±${HALF_WIN}bp window with alt swap)..."

# Build BED of B1+B2 positions with their ref/alt
awk -F'\t' 'NR>1 && ($9=="B1" || $9=="B2" || $9=="B4"){
    start = $2 - 1 - '$HALF_WIN'
    end   = $2 + '$HALF_WIN'
    if (start < 0) start = 0
    print $1"\t"start"\t"end"\t"$1"_"$2"_"$3"_"$4
}' $DATA/bucket_full.tsv > $DATA/_mini_ref.bed
echo "  Windows to build: $(wc -l < $DATA/_mini_ref.bed)"

# Extract reference sequence for each window
bedtools getfasta -fi $REF -bed $DATA/_mini_ref.bed -name -tab > $DATA/_mini_ref_seq.tsv 2>$LOGS/getfasta_07.stderr

# Swap center to alt
python3 - <<PYEOF
import os
HALF = $HALF_WIN
out = open("$DATA/mini_ref.fa", "w")
n_built = 0; n_skip = 0
with open("$DATA/_mini_ref_seq.tsv") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        name = parts[0].split("::")[0]
        seq = parts[1].upper()
        if len(seq) != 2*HALF + 1:
            n_skip += 1; continue
        toks = name.split("_")
        if len(toks) < 4:
            n_skip += 1; continue
        chrom = "_".join(toks[:-3]); pos = toks[-3]; refb = toks[-2]; altb = toks[-1]
        if seq[HALF] != refb:
            n_skip += 1; continue
        if altb not in "ACGT" or len(altb) != 1:
            n_skip += 1; continue
        if "N" in seq:
            n_skip += 1; continue
        alt_seq = seq[:HALF] + altb + seq[HALF+1:]
        # Header: chrom_pos_ref_alt (unique window name)
        out.write(f">{chrom}_{pos}_{refb}_{altb}\n{alt_seq}\n")
        n_built += 1
out.close()
print(f"  Mini-ref built: {n_built} windows (skipped {n_skip})")
PYEOF
rm $DATA/_mini_ref.bed $DATA/_mini_ref_seq.tsv

echo "  Mini-ref FASTA size: $(du -h $DATA/mini_ref.fa | cut -f1)"

# ---- 2. BWA index ----
echo "[2/6] Building BWA index..."
bwa index $DATA/mini_ref.fa 2>$LOGS/bwa_index.stderr
samtools faidx $DATA/mini_ref.fa
echo "  Index files: $(ls $DATA/mini_ref.fa.* | wc -l)"

# ---- 3. Extract unmapped reads (FASTQ) from all 10 BAMs ----
echo "[3/6] Extracting unmapped reads from 10 treatment BAMs..."
UFQ=$DATA/unmapped_all.fq
> $UFQ
for s in P1 P2 P3 P4 P5 P6 P7 P8 P9 P10; do
    bam=${BAMS[$s]}
    echo "  $s ($(basename $bam))"
    # -f 4 unmapped; -F 256 primary only
    # Convert to FASTQ. Tag read name with sample prefix to keep traceable.
    samtools view -f 4 -F 256 "$bam" 2>/dev/null \
      | awk -v s=$s 'BEGIN{OFS="\n"}
        {
            # SAM cols: QNAME=1, SEQ=10, QUAL=11
            qname = $1
            seq   = $10
            qual  = $11
            if (length(seq) >= 50) print "@"s"_"qname"_"NR, seq, "+", qual
        }' >> $UFQ
done
n_reads=$(awk 'NR%4==1' $UFQ | wc -l)
echo "  Total unmapped reads pooled: $n_reads"

# ---- 4. BWA mem alignment ----
echo "[4/6] BWA mem unmapped vs mini-ref (multi-threaded)..."
THREADS=$(nproc 2>/dev/null || echo 4)
bwa mem -t $THREADS -k 19 -T 30 $DATA/mini_ref.fa $UFQ 2>$LOGS/bwa_mem.stderr \
  | samtools view -bS - > $DATA/mini_ref.bam
echo "  Alignment done. BAM: $(du -h $DATA/mini_ref.bam | cut -f1)"

# ---- 5. Parse alignments: count reads per window that cover center with alt ----
echo "[5/6] Parsing alignments..."
# Filter: MAPQ >= 20 (unique-ish), mapped (not flag 4), primary (not flag 256/2048)
# For each aligned read, check if alignment covers the center position of the window.
# Window is 201bp (HALF_WIN=100, center=position 101 1-based).
# Read RNAME = window name; POS = 1-based leftmost mapping. CIGAR gives length.
samtools view -F 260 -q 20 $DATA/mini_ref.bam 2>/dev/null \
  | awk -v half=$HALF_WIN -v ctr=$CENTER 'BEGIN{OFS="\t"}
    {
        rname = $3
        pos   = $4
        cigar = $6
        # Compute read length on reference
        n = 0
        ref_len = 0
        cig_copy = cigar
        # Parse CIGAR
        while (match(cig_copy, /^[0-9]+[MIDNSHPX=]/)) {
            tok = substr(cig_copy, RSTART, RLENGTH)
            l = substr(tok, 1, RLENGTH-1) + 0
            op = substr(tok, RLENGTH, 1)
            if (op ~ /[MDN=X]/) ref_len += l
            cig_copy = substr(cig_copy, RLENGTH+1)
        }
        end = pos + ref_len - 1
        # Does the alignment cover the center?
        if (pos <= ctr && end >= ctr) {
            print rname
        }
    }' \
  | sort | uniq -c | awk 'BEGIN{OFS="\t"; print "window","reads_covering_alt"} {print $2, $1}' \
  > $DATA/augmented_hits.tsv

n_windows_hit=$(($(wc -l < $DATA/augmented_hits.tsv) - 1))
echo "  Windows with >=1 alt-covering read: $n_windows_hit"

# ---- 6. Join with B4 + bucket info and summarize ----
echo "[6/6] Cross-check with k-mer salvage + buckets..."
python3 - <<PYEOF
from collections import Counter, defaultdict

# Load augmented hits
aug = {}  # (chrom, pos, ref, alt) -> read count
with open("$DATA/augmented_hits.tsv") as f:
    next(f)
    for line in f:
        win, n = line.rstrip("\n").split("\t")
        toks = win.split("_")
        chrom = "_".join(toks[:-3]); pos = toks[-3]; ref = toks[-2]; alt = toks[-1]
        aug[(chrom, pos, ref, alt)] = int(n)

# Load current bucket
bucket = {}
with open("$DATA/bucket_full.tsv") as f:
    next(f)
    for line in f:
        parts = line.rstrip("\n").split("\t")
        bucket[(parts[0], parts[1], parts[2], parts[3])] = parts[8]

# Load robust B4 (k-mer salvage)
robust = set()
with open("$DATA/b4_robust.tsv") as f:
    next(f)
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if parts[8] == "KEPT":
            robust.add((parts[0], parts[1], parts[2], parts[3]))

# Tier B (k-mer + imputation filter)
tier_B = set()
with open("$DATA/salvage_tier_B_moderate.tsv") as f:
    next(f)
    for line in f:
        parts = line.rstrip("\n").split("\t")
        tier_B.add((parts[0], parts[1], parts[2], parts[3]))

# Stats
THRESHOLDS = [1, 2, 3, 5, 10]
print(f"Augmented realignment positive (>=1 read): {len(aug)}")
print()
print("=== Augmented evidence stratified by read count ===")
for t in THRESHOLDS:
    n = sum(1 for v in aug.values() if v >= t)
    print(f"  >= {t} reads: {n}")

# Cross-tab: bucket vs augmented evidence
print()
print("=== Bucket × augmented evidence cross-tab ===")
cross = defaultdict(lambda: Counter())
for key, b in bucket.items():
    n_aug = aug.get(key, 0)
    if n_aug == 0: bin_ = "0"
    elif n_aug < 2: bin_ = "1"
    elif n_aug < 5: bin_ = "2-4"
    elif n_aug < 10: bin_ = "5-9"
    else: bin_ = ">=10"
    cross[b][bin_] += 1

print(f"{'Bucket':<10} {'aug=0':>10} {'aug=1':>10} {'aug=2-4':>10} {'aug=5-9':>10} {'aug>=10':>10}")
for b in ["B1","B2","B3","B4"]:
    c = cross[b]
    print(f"{b:<10} {c.get('0',0):>10} {c.get('1',0):>10} {c.get('2-4',0):>10} {c.get('5-9',0):>10} {c.get('>=10',0):>10}")

# Cross with k-mer salvage tiers
print()
print("=== K-mer salvage tier × augmented confirmation ===")
for t in THRESHOLDS:
    aug_pos = {k for k, v in aug.items() if v >= t}
    n_robust_conf = len(robust & aug_pos)
    n_tier_b_conf = len(tier_B & aug_pos)
    print(f"  augmented >= {t}: robust_B4 confirmed = {n_robust_conf}/{len(robust)}; tier_B confirmed = {n_tier_b_conf}/{len(tier_B)}")

# Positions newly identified by augmented (not in k-mer salvage)
new_positions = {k for k, v in aug.items() if v >= 2} - robust
b1_b2_new = {k for k in new_positions if bucket.get(k) in ("B1", "B2")}
b3_new    = {k for k in new_positions if bucket.get(k) == "B3"}
print()
print(f"=== Positions confirmed by augmented (>=2 reads) NOT in k-mer robust salvage ===")
print(f"  Total: {len(new_positions)}")
print(f"  From B1/B2 (true new salvage candidates): {len(b1_b2_new)}")
print(f"  From B3 (also covered in mapped pool): {len(b3_new)}")

# Save augmented-confirmed positions for downstream
with open("$DATA/augmented_confirmed.tsv", "w") as f:
    f.write("chrom\tpos\tref\talt\taug_read_count\torig_bucket\tin_kmer_robust\tin_tier_B\n")
    for k, v in sorted(aug.items(), key=lambda x: -x[1]):
        if v >= 2:
            b = bucket.get(k, "?")
            in_r = "Y" if k in robust else "N"
            in_b = "Y" if k in tier_B else "N"
            f.write(f"{k[0]}\t{k[1]}\t{k[2]}\t{k[3]}\t{v}\t{b}\t{in_r}\t{in_b}\n")

print()
print(f"Saved: $DATA/augmented_confirmed.tsv")
PYEOF

# Cleanup large intermediate
gzip -f $DATA/unmapped_all.fq
gzip -f $DATA/mini_ref.fa 2>/dev/null || true
# Keep BAM compressed already; index it for inspection
samtools sort -@ 4 $DATA/mini_ref.bam -o $DATA/mini_ref.sorted.bam 2>/dev/null
mv $DATA/mini_ref.sorted.bam $DATA/mini_ref.bam
samtools index $DATA/mini_ref.bam 2>/dev/null

echo "Done. See $DATA/augmented_confirmed.tsv and $DATA/augmented_hits.tsv"
