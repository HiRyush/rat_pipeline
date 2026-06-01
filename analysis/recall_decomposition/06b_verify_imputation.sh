#!/bin/bash
# 06b_verify_imputation.sh — cross-check salvaged B4 positions against imputed VCF.
# If imputation called the position as 0|1 or 1|1 in any CONTROL sample → germline-like,
# should NOT be claimed as somatic.

set -euo pipefail

ROOT=/home/yusanghyeon/RAT_project/PHMG_IT
BASE=$ROOT/results/recall_decomposition
DATA=$BASE/data
RES=$BASE/results

source $ROOT/../miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

IMP=$ROOT/results/imputed_differential/all15_imputed_merged.vcf.gz

# Build BED of salvaged (KEPT) positions
awk -F'\t' 'NR>1 && $9=="KEPT"{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4}' $DATA/b4_robust.tsv \
    > $DATA/salvaged_positions.bed
echo "Salvaged positions: $(wc -l < $DATA/salvaged_positions.bed)"

echo "[1/2] Querying imputed VCF at salvaged positions..."
# Use bcftools view -R to extract records at these positions
bcftools view -R $DATA/salvaged_positions.bed -Ou $IMP 2>/dev/null \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT,]\n' \
  > $DATA/salvaged_imputed_gt.tsv

echo "  Imputed records found at salvaged positions: $(wc -l < $DATA/salvaged_imputed_gt.tsv)"

echo "[2/2] Cross-check: count salvaged that show germline-like imputation in any CONTROL..."
# Sample order in merged VCF: C1 C2 C3 C4 C5 P1 P2 ... P10
# That's the order shown in the VCF header. Check first 5 GT (control) for non-ref.
python3 - <<'PYEOF'
from collections import Counter

# Load salvaged positions
salvaged = {}  # (chrom, pos, ref, alt) -> hits
with open("/home/yusanghyeon/RAT_project/PHMG_IT/results/recall_decomposition/data/b4_robust.tsv") as f:
    next(f)
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if parts[8] == "KEPT":
            salvaged[(parts[0], parts[1], parts[2], parts[3])] = int(parts[4])

# Parse imputed GT — order: C1..C5, P1..P10
def has_alt(gt):
    return "1" in gt and gt not in (".", "./.", ".|.", "")

ctrl_alt_cnt = Counter()   # number of salvaged positions with N control samples having ALT
treat_alt_cnt = Counter()
matched = 0
unmatched = set(salvaged.keys())

with open("/home/yusanghyeon/RAT_project/PHMG_IT/results/recall_decomposition/data/salvaged_imputed_gt.tsv") as f:
    for line in f:
        chrom, pos, ref, alt, gts = line.rstrip("\n").split("\t")
        # Handle multi-allelic ALT
        for actual_alt in alt.split(","):
            key = (chrom, pos, ref, actual_alt)
            if key in salvaged:
                matched += 1
                unmatched.discard(key)
                gt_list = [g for g in gts.split(",") if g]
                if len(gt_list) >= 15:
                    ctrl_gts = gt_list[:5]
                    treat_gts = gt_list[5:15]
                    n_ctrl_alt = sum(1 for g in ctrl_gts if has_alt(g))
                    n_treat_alt = sum(1 for g in treat_gts if has_alt(g))
                    ctrl_alt_cnt[n_ctrl_alt] += 1
                    treat_alt_cnt[n_treat_alt] += 1

print(f"Salvaged positions: {len(salvaged)}")
print(f"Matched in imputed VCF: {matched}")
print(f"Unmatched (no imputation record): {len(unmatched)}")
print()
print("=== Control samples with imputed ALT at salvaged positions ===")
for n in sorted(ctrl_alt_cnt):
    print(f"  {n} controls with ALT: {ctrl_alt_cnt[n]} positions")
print()
print("=== Treatment samples with imputed ALT at salvaged positions ===")
for n in sorted(treat_alt_cnt):
    print(f"  {n} treats with ALT: {treat_alt_cnt[n]} positions")

# Concerning: positions where imputation called ANY control as ALT (germline-like)
n_concerning = sum(v for k, v in ctrl_alt_cnt.items() if k > 0)
print()
print(f"⚠️  Concerning (>=1 control imputed as ALT → germline-like): {n_concerning} / {matched} ({n_concerning*100/max(matched,1):.1f}%)")
print(f"✅  Clean (no control imputed as ALT): {ctrl_alt_cnt.get(0,0)} / {matched} ({ctrl_alt_cnt.get(0,0)*100/max(matched,1):.1f}%)")

# Save flagged positions for follow-up
with open("/home/yusanghyeon/RAT_project/PHMG_IT/results/recall_decomposition/data/salvaged_imputed_flagged.tsv", "w") as out:
    out.write("chrom\tpos\tref\talt\tn_ctrl_alt\tn_treat_alt\timputed_GT\n")
    with open("/home/yusanghyeon/RAT_project/PHMG_IT/results/recall_decomposition/data/salvaged_imputed_gt.tsv") as f:
        for line in f:
            chrom, pos, ref, alt, gts = line.rstrip("\n").split("\t")
            for actual_alt in alt.split(","):
                key = (chrom, pos, ref, actual_alt)
                if key in salvaged:
                    gt_list = [g for g in gts.split(",") if g]
                    if len(gt_list) >= 15:
                        n_c = sum(1 for g in gt_list[:5] if has_alt(g))
                        n_t = sum(1 for g in gt_list[5:15] if has_alt(g))
                        if n_c > 0 or n_t > 0:
                            out.write(f"{chrom}\t{pos}\t{ref}\t{actual_alt}\t{n_c}\t{n_t}\t{gts}\n")
PYEOF

echo "Done. Flagged positions saved to $DATA/salvaged_imputed_flagged.tsv"
