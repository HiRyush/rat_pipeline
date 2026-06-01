#!/usr/bin/env python3
"""
06c_salvage_with_imputation_filter.py — Apply imputation cross-check to salvage.

Salvage tiers:
  A (strict):   imputation says NO control sample has ALT  (n_ctrl_alt = 0 in imputed VCF)
  B (moderate): A ∪ {positions with no imputation record} (panel-novel, assumed rare)
  C (raw):      all 1,476 robust B4 (no imputation filter)

Computes recall curves for each tier and writes summary.
"""
import os, sys, datetime, gzip
from collections import Counter

BASE = "/home/yusanghyeon/RAT_project/PHMG_IT/results/recall_decomposition"
DATA = f"{BASE}/data"
RES  = f"{BASE}/results"

# Load robust B4 (KEPT)
robust = set()
b4_meta = {}
with open(f"{DATA}/b4_robust.tsv") as f:
    next(f)
    for line in f:
        parts = line.rstrip("\n").split("\t")
        chrom, pos, ref, alt, hits, kmer, ent, hom, status, reason = parts
        if status == "KEPT":
            key = (chrom, pos, ref, alt)
            robust.add(key)
            b4_meta[key] = (int(hits), float(ent))

# Parse imputed GT records
def has_alt(g):
    return "1" in g and g not in (".", "./.", ".|.")

imputed_ctrl_alt = {}   # key -> n_ctrl_alt
imputed_records = set()
with open(f"{DATA}/salvaged_imputed_gt.tsv") as f:
    for line in f:
        chrom, pos, ref, alt, gts = line.rstrip("\n").split("\t")
        for actual_alt in alt.split(","):
            key = (chrom, pos, ref, actual_alt)
            if key in robust:
                imputed_records.add(key)
                gt_list = [g for g in gts.split(",") if g]
                if len(gt_list) >= 15:
                    imputed_ctrl_alt[key] = sum(1 for g in gt_list[:5] if has_alt(g))

# Build tiers
tier_A = {k for k in robust if k in imputed_records and imputed_ctrl_alt.get(k, 0) == 0}
tier_B_only_novel = {k for k in robust if k not in imputed_records}
tier_B = tier_A | tier_B_only_novel
tier_C = robust

print(f"Robust B4 (KEPT): {len(robust)}")
print(f"  In imputation panel: {len(imputed_records)}")
print(f"  Panel-novel: {len(tier_B_only_novel)}")
print(f"  Tier A (strict, clean ctrl): {len(tier_A)}")
print(f"  Tier B (moderate, A + novel): {len(tier_B)}")
print(f"  Tier C (aggressive, all): {len(tier_C)}")

# Load captured baseline
captured = set()
with open(f"{DATA}/captured.bed") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        captured.add((parts[0], parts[2], parts[3], parts[4]))
n_cap = len(captured)
print(f"\nBaseline captured: {n_cap}")

# Load B3 with AF
b3_rows = []  # (alt_dp, total_dp)
with open(f"{DATA}/bucket_full.tsv") as f:
    next(f)
    for line in f:
        parts = line.rstrip("\n").split("\t")
        try:
            if parts[8] == "B3":
                dp = int(parts[4]); alt_dp = int(parts[5])
                b3_rows.append((alt_dp, dp))
        except ValueError:
            continue

af_thresholds = [0.05, 0.10, 0.20, 0.30]
dp_thresholds = [1, 2, 3, 5, 10]

def b3_count(thr_kind, thr):
    if thr_kind == 'dp':
        return sum(1 for a, d in b3_rows if a >= thr)
    elif thr_kind == 'af':
        return sum(1 for a, d in b3_rows if d > 0 and (a/d) >= thr)

# Recall under each tier
print("\n=== Recall comparison ===")
print(f"{'Operating point':<20} {'baseline':>12} {'tier A':>12} {'tier B':>12} {'tier C':>12}")
print(f"{'(numerator)':<20} {n_cap:>12} {n_cap+len(tier_A):>12} {n_cap+len(tier_B):>12} {n_cap+len(tier_C):>12}")
print("-"*72)

out_rows = []
for af in af_thresholds:
    n_b3 = b3_count('af', af)
    r_b = n_cap / (n_cap + n_b3) * 100
    r_A = (n_cap + len(tier_A)) / (n_cap + len(tier_A) + n_b3) * 100
    r_B = (n_cap + len(tier_B)) / (n_cap + len(tier_B) + n_b3) * 100
    r_C = (n_cap + len(tier_C)) / (n_cap + len(tier_C) + n_b3) * 100
    print(f"AF>={af}            {r_b:>11.2f}% {r_A:>11.2f}% {r_B:>11.2f}% {r_C:>11.2f}%")
    out_rows.append(('af', af, n_b3, r_b, r_A, r_B, r_C))

print()
for dp in dp_thresholds:
    n_b3 = b3_count('dp', dp)
    r_b = n_cap / (n_cap + n_b3) * 100
    r_A = (n_cap + len(tier_A)) / (n_cap + len(tier_A) + n_b3) * 100
    r_B = (n_cap + len(tier_B)) / (n_cap + len(tier_B) + n_b3) * 100
    r_C = (n_cap + len(tier_C)) / (n_cap + len(tier_C) + n_b3) * 100
    print(f"alt_DP>={dp:<3}            {r_b:>11.2f}% {r_A:>11.2f}% {r_B:>11.2f}% {r_C:>11.2f}%")
    out_rows.append(('dp', dp, n_b3, r_b, r_A, r_B, r_C))

# Write TSV
with open(f"{RES}/06c_recall_by_tier.tsv", "w") as f:
    f.write("threshold_type\tthreshold\tb3_size\trecall_baseline_pct\trecall_tierA_pct\trecall_tierB_pct\trecall_tierC_pct\n")
    for tt, t, b3, rb, ra, rB, rC in out_rows:
        f.write(f"{tt}\t{t}\t{b3}\t{rb:.2f}\t{ra:.2f}\t{rB:.2f}\t{rC:.2f}\n")

# Save final numerator lists for each tier
for name, tier in [("A_strict", tier_A), ("B_moderate", tier_B), ("C_aggressive", tier_C)]:
    with open(f"{DATA}/salvage_tier_{name}.tsv", "w") as f:
        f.write("chrom\tpos\tref\talt\thits\tentropy\n")
        for k in tier:
            h, e = b4_meta[k]
            f.write(f"{k[0]}\t{k[1]}\t{k[2]}\t{k[3]}\t{h}\t{e:.2f}\n")

# Save summary
ts = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
with open(f"{RES}/06c_salvage_tiers_summary.txt", "w") as f:
    f.write(f"Salvage tier summary\nGenerated: {ts}\n\n")
    f.write(f"Tier A (strict, imputation says 0 ctrl ALT): {len(tier_A)}\n")
    f.write(f"Tier B (moderate, A + panel-novel): {len(tier_B)}\n")
    f.write(f"Tier C (aggressive, no imputation filter): {len(tier_C)}\n")
    f.write(f"Flagged (>=1 ctrl ALT in imputation, dropped from A/B): {len(robust) - len(tier_B)}\n\n")
    f.write(f"Baseline captured: {n_cap}\n")
    f.write(f"Expanded numerators: A={n_cap+len(tier_A)}, B={n_cap+len(tier_B)}, C={n_cap+len(tier_C)}\n\n")
    f.write(f"Headline recall @ alt_AF>=0.30:\n")
    n_b3_af30 = b3_count('af', 0.30)
    f.write(f"  baseline: {n_cap}/{n_cap+n_b3_af30} = {n_cap/(n_cap+n_b3_af30)*100:.2f}%\n")
    f.write(f"  tier A:   {n_cap+len(tier_A)}/{n_cap+len(tier_A)+n_b3_af30} = {(n_cap+len(tier_A))/(n_cap+len(tier_A)+n_b3_af30)*100:.2f}%\n")
    f.write(f"  tier B:   {n_cap+len(tier_B)}/{n_cap+len(tier_B)+n_b3_af30} = {(n_cap+len(tier_B))/(n_cap+len(tier_B)+n_b3_af30)*100:.2f}%\n")
    f.write(f"  tier C:   {n_cap+len(tier_C)}/{n_cap+len(tier_C)+n_b3_af30} = {(n_cap+len(tier_C))/(n_cap+len(tier_C)+n_b3_af30)*100:.2f}%\n")

print(f"\nSaved: {RES}/06c_recall_by_tier.tsv, {RES}/06c_salvage_tiers_summary.txt")
