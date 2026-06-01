#!/bin/bash
# 09_full_imputation_filter.sh
# Query imputed VCF for ALL aug>=2 positions (not just k-mer subset) and recompute
# a properly imputation-filtered Aug-strict salvage set.

set -euo pipefail

ROOT=/home/yusanghyeon/RAT_project/PHMG_IT
BASE=$ROOT/results/recall_decomposition
DATA=$BASE/data
RES=$BASE/results

source $ROOT/../miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

IMP=$ROOT/results/imputed_differential/all15_imputed_merged.vcf.gz

# 1. Build BED of all aug>=2 positions (universe = augmented_confirmed.tsv, which is aug>=2)
echo "[1/3] Building BED of all aug>=2 positions..."
awk 'BEGIN{OFS="\t"} NR>1 {print $1, $2-1, $2, $3, $4}' $DATA/augmented_confirmed.tsv \
    | sort -k1,1 -k2,2n -u > $DATA/aug_all_positions.bed
echo "  Positions: $(wc -l < $DATA/aug_all_positions.bed)"

# 2. Query imputed VCF
echo "[2/3] Querying imputed VCF (15 samples: C1-C5, P1-P10)..."
bcftools view -R $DATA/aug_all_positions.bed -Ou $IMP 2>/dev/null \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT,]\n' \
  > $DATA/aug_all_imputed_gt.tsv
echo "  Imputed records: $(wc -l < $DATA/aug_all_imputed_gt.tsv)"

# 3. Classify by ctrl/treat ALT counts and apply filter
echo "[3/3] Classifying and recomputing Aug-strict..."
python3 - <<'PYEOF'
from collections import Counter

# Load aug>=2 positions (key -> aug_reads)
aug_pos = {}
with open("/home/yusanghyeon/RAT_project/PHMG_IT/results/recall_decomposition/data/augmented_confirmed.tsv") as f:
    next(f)
    for line in f:
        parts = line.rstrip("\n").split("\t")
        aug_pos[(parts[0], parts[1], parts[2], parts[3])] = int(parts[4])

def has_alt(gt):
    return "1" in gt and gt not in (".", "./.", ".|.")

# Sample order in merged VCF: C1-C5, P1-P10
imputed = {}  # key -> (n_ctrl_alt, n_treat_alt) or None if no record
with open("/home/yusanghyeon/RAT_project/PHMG_IT/results/recall_decomposition/data/aug_all_imputed_gt.tsv") as f:
    for line in f:
        chrom, pos, ref, alt, gts = line.rstrip("\n").split("\t")
        for actual_alt in alt.split(","):
            key = (chrom, pos, ref, actual_alt)
            if key in aug_pos:
                gt_list = [g for g in gts.split(",") if g]
                if len(gt_list) >= 15:
                    n_c = sum(1 for g in gt_list[:5] if has_alt(g))
                    n_t = sum(1 for g in gt_list[5:15] if has_alt(g))
                    imputed[key] = (n_c, n_t)

# Tier breakdown
in_panel = set(imputed.keys())
panel_novel = set(aug_pos.keys()) - in_panel
clean_in_panel = {k for k, v in imputed.items() if v[0] == 0}
germline_like = {k for k, v in imputed.items() if v[0] > 0}

print(f"Total aug>=2 positions: {len(aug_pos)}")
print(f"  Found in imputation panel: {len(in_panel)}")
print(f"    Clean (0 ctrl ALT): {len(clean_in_panel)}")
print(f"    Germline-like (≥1 ctrl ALT): {len(germline_like)}")
print(f"  Panel-novel (no imputation record): {len(panel_novel)}")
print()

# New tier definitions (proper imputation filter on ALL aug positions)
aug_strict_v2 = clean_in_panel | panel_novel   # everything NOT imputation-flagged
aug_strict_strict = clean_in_panel              # only positions explicitly in panel and clean

print(f"=== Aug-strict v2 (proper full imputation filter) ===")
print(f"  aug_strict_v2 (clean + panel-novel): {len(aug_strict_v2)}")
print(f"  aug_strict_strict (clean in panel only): {len(aug_strict_strict)}")
print()

# Save
DATA = "/home/yusanghyeon/RAT_project/PHMG_IT/results/recall_decomposition/data"
with open(f"{DATA}/aug_strict_v2.tsv", "w") as f:
    f.write("chrom\tpos\tref\talt\taug_reads\timputation_status\n")
    for k in sorted(aug_strict_v2):
        c, p, r, a = k
        if k in clean_in_panel:
            status = "in_panel_clean"
        else:
            status = "panel_novel"
        f.write(f"{c}\t{p}\t{r}\t{a}\t{aug_pos[k]}\t{status}\n")

with open(f"{DATA}/aug_germline_flagged.tsv", "w") as f:
    f.write("chrom\tpos\tref\talt\taug_reads\tn_ctrl_alt\tn_treat_alt\n")
    for k in sorted(germline_like):
        c, p, r, a = k
        nc, nt = imputed[k]
        f.write(f"{c}\t{p}\t{r}\t{a}\t{aug_pos[k]}\t{nc}\t{nt}\n")

# Recompute recall under new aug_strict
captured = set()
with open(f"{DATA}/captured.bed") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        captured.add((parts[0], parts[2], parts[3], parts[4]))
n_cap = len(captured)

b3_rows = []
with open(f"{DATA}/bucket_full.tsv") as f:
    next(f)
    for line in f:
        parts = line.rstrip("\n").split("\t")
        try:
            if parts[8] == "B3":
                b3_rows.append((int(parts[5]), int(parts[4])))
        except ValueError:
            continue

def b3_af(af): return sum(1 for a, d in b3_rows if d>0 and a/d >= af)
def b3_dp(dp): return sum(1 for a, d in b3_rows if a >= dp)

def recall(num, b3): return num / (num + b3) * 100

print(f"=== Recall comparison (with FULL imputation filter on aug) ===")
print(f"{'Operating point':<18} {'baseline':>10} {'old aug-strict':>14} {'NEW aug-strict v2':>18} {'NEW aug-strict (panel only)':>28}")
n_old = 1962  # prior incomplete
n_v2 = len(aug_strict_v2)
n_strict = len(aug_strict_strict)
print(f"{'(numerator)':<18} {n_cap:>10} {n_cap+n_old:>14} {n_cap+n_v2:>18} {n_cap+n_strict:>28}")
print()
for af in [0.05, 0.10, 0.20, 0.30]:
    b3 = b3_af(af)
    r_b = recall(n_cap, b3)
    r_old = recall(n_cap+n_old, b3)
    r_v2 = recall(n_cap+n_v2, b3)
    r_strict = recall(n_cap+n_strict, b3)
    print(f"AF>={af:<4}           {r_b:>9.2f}% {r_old:>13.2f}% {r_v2:>17.2f}% {r_strict:>27.2f}%")
print()
for dp in [1, 2, 3, 5, 10]:
    b3 = b3_dp(dp)
    r_b = recall(n_cap, b3)
    r_old = recall(n_cap+n_old, b3)
    r_v2 = recall(n_cap+n_v2, b3)
    r_strict = recall(n_cap+n_strict, b3)
    print(f"alt_DP>={dp:<5}         {r_b:>9.2f}% {r_old:>13.2f}% {r_v2:>17.2f}% {r_strict:>27.2f}%")

# Save matrix
with open(f"/home/yusanghyeon/RAT_project/PHMG_IT/results/recall_decomposition/results/09_recall_full_imputation.tsv", "w") as f:
    f.write("threshold_type\tthreshold\tb3_size\tbaseline\told_aug_strict\tnew_aug_strict_v2\tnew_aug_strict_panel_only\n")
    for tt, vals in [("af", [0.05, 0.10, 0.20, 0.30]), ("dp", [1, 2, 3, 5, 10])]:
        for t in vals:
            b3 = b3_af(t) if tt == "af" else b3_dp(t)
            f.write(f"{tt}\t{t}\t{b3}\t{recall(n_cap,b3):.2f}\t{recall(n_cap+n_old,b3):.2f}\t{recall(n_cap+n_v2,b3):.2f}\t{recall(n_cap+n_strict,b3):.2f}\n")

print(f"\nSaved: aug_strict_v2.tsv ({n_v2}), aug_germline_flagged.tsv ({len(germline_like)}), 09_recall_full_imputation.tsv")
PYEOF
