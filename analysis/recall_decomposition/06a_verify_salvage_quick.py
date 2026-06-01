#!/usr/bin/env python3
"""
06a_verify_salvage_quick.py — Quick verification of 1,476 salvaged B4.
Outputs distributional + DNA-recurrence diagnostics (no external tool calls).
"""
import sys, gzip, os
from collections import Counter

BASE = "/home/yusanghyeon/RAT_project/PHMG_IT/results/recall_decomposition"
DATA = f"{BASE}/data"
RES  = f"{BASE}/results"
DT_DIR = "/home/yusanghyeon/RAT_project/PHMG_IT/results/dna_truth_coverage"

# Load robust B4 (KEPT only) — key (chrom, pos, ref, alt)
keep = set()
b4_records = {}  # key -> (hits, entropy, homopolymer)
with open(f"{DATA}/b4_robust.tsv") as f:
    next(f)
    for line in f:
        parts = line.rstrip("\n").split("\t")
        chrom, pos, ref, alt, hits, kmer, ent, hom, status, reason = parts
        if status == "KEPT":
            keep.add((chrom, pos, ref, alt))
            b4_records[(chrom, pos, ref, alt)] = (int(hits), float(ent), int(hom))

print(f"Robust B4 (KEPT) total: {len(keep)}")

# 1. Origin bucket (B1 vs B2)
print("\n=== (a) Origin bucket (B1 vs B2) ===")
origin = Counter()
with open(f"{DATA}/bucket_full.tsv") as f:
    next(f)
    for line in f:
        parts = line.rstrip("\n").split("\t")
        chrom, pos, ref, alt = parts[0], parts[1], parts[2], parts[3]
        b123 = parts[6]
        final = parts[8]
        if final == "B4":
            key = (chrom, pos, ref, alt)
            if key in keep:
                origin[b123] += 1
for k, n in origin.most_common():
    print(f"  Originally {k}: {n}  ({n*100/len(keep):.1f}%)")

# 2. Hit count / entropy distribution
print("\n=== (b) Hit count + entropy distribution among KEPT ===")
hits_list = sorted([v[0] for v in b4_records.values()])
ent_list  = sorted([v[1] for v in b4_records.values()])

def quantiles(lst, qs=(0.10, 0.25, 0.50, 0.75, 0.90)):
    return {q: lst[int(len(lst)*q)] for q in qs}

print("  Hits:", quantiles(hits_list), "  min/max:", hits_list[0], "/", hits_list[-1])
print(f"  Entropy: P10={ent_list[int(len(ent_list)*0.1)]:.2f}  P50={ent_list[int(len(ent_list)*0.5)]:.2f}  P90={ent_list[int(len(ent_list)*0.9)]:.2f}")

# 3. DNA-side treat_n recurrence among salvaged
print("\n=== (c) DNA-side treat_n recurrence ===")
dna_treat_n = {}
with gzip.open(f"{DT_DIR}/dna_truth_counts_OptA.tsv.gz", "rt") as f:
    next(f)
    for line in f:
        chrom, pos, ref, alt, c, t = line.rstrip("\n").split("\t")
        key = (chrom, pos, ref, alt)
        if key in keep:
            dna_treat_n[key] = int(t)

treat_dist = Counter(dna_treat_n.values())
print(f"  Salvaged positions found in DNA truth Option A: {len(dna_treat_n)}/{len(keep)}")
total = 0
for k in sorted(treat_dist):
    print(f"  treat_n={k}: {treat_dist[k]}")
    total += treat_dist[k]
mean_t = sum(k*v for k,v in treat_dist.items()) / total
print(f"  Mean treat_n among salvaged: {mean_t:.2f}")

# 4. AF in mapped reads at salvaged positions (should be 0 by B1/B2 definition; sanity check)
print("\n=== (d) Verify mapped-read alt evidence is indeed 0 (B1/B2 invariant) ===")
n_violations = 0
with open(f"{DATA}/bucket_full.tsv") as f:
    next(f)
    for line in f:
        parts = line.rstrip("\n").split("\t")
        chrom, pos, ref, alt = parts[0], parts[1], parts[2], parts[3]
        dp, alt_dp = int(parts[4]), int(parts[5])
        if (chrom, pos, ref, alt) in keep:
            if alt_dp > 0:
                n_violations += 1
print(f"  Violations (alt_dp > 0 among salvaged): {n_violations}  (should be 0)")

# 5. Chromosome distribution (any single-chrom artifact?)
print("\n=== (e) Chromosome distribution ===")
chrom_dist = Counter(k[0] for k in keep)
for c, n in chrom_dist.most_common(15):
    print(f"  {c}: {n}")

# Save summary
with open(f"{RES}/06a_verify_salvage_quick.txt", "w") as f:
    f.write(f"Robust B4 (KEPT): {len(keep)}\n\n")
    f.write("(a) Origin bucket:\n")
    for k, n in origin.most_common():
        f.write(f"  {k}: {n} ({n*100/len(keep):.1f}%)\n")
    f.write("\n(b) Hit count / entropy distribution:\n")
    f.write(f"  Hits: min={hits_list[0]}, P10={hits_list[int(len(hits_list)*0.1)]}, P50={hits_list[len(hits_list)//2]}, P90={hits_list[int(len(hits_list)*0.9)]}, max={hits_list[-1]}\n")
    f.write(f"  Entropy: P10={ent_list[int(len(ent_list)*0.1)]:.2f}, P50={ent_list[len(ent_list)//2]:.2f}, P90={ent_list[int(len(ent_list)*0.9)]:.2f}\n")
    f.write("\n(c) DNA-side treat_n recurrence:\n")
    for k in sorted(treat_dist):
        f.write(f"  treat_n={k}: {treat_dist[k]}\n")
    f.write(f"  Mean: {mean_t:.2f}\n")
    f.write(f"\n(d) Sanity: B1/B2 alt_dp violations among salvaged: {n_violations} (expected 0)\n")
    f.write("\n(e) Chromosome distribution (top 15):\n")
    for c, n in chrom_dist.most_common(15):
        f.write(f"  {c}: {n}\n")

print(f"\nSaved: {RES}/06a_verify_salvage_quick.txt")
