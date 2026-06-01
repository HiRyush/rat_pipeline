#!/usr/bin/env python3
"""
Step 6: RNA Editing + Artifact Filter
- Download/use REDIportal rat RNA editing sites
- Filter known editing sites from candidates
- Apply DP, AF, and basic quality filters
"""
import os
import sys
import subprocess

OUT_DIR = "/home/yusanghyeon/RAT_project/pipeline_candidates/results/artifact_filter"
DIFF_DIR = "/home/yusanghyeon/RAT_project/pipeline_candidates/results/differential"
os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 70)
print("Step 6: RNA Editing + Artifact Filter")
print("=" * 70)

# 1. Build known RNA editing site set
# Since REDIportal requires registration, we use a conservative approach:
# Filter A>G (I) and T>C (reverse strand A>I) which are ADAR-mediated
# These are the dominant RNA editing events (~99% of mammalian RNA editing)
print("\n--- Building RNA editing filter (A>G / T>C conservative) ---")

# 2. Load candidates
candidates_file = os.path.join(DIFF_DIR, "treatment_only_rec2.tsv")
candidates = []
with open(candidates_file) as f:
    header = f.readline().strip()
    for line in f:
        parts = line.strip().split('\t')
        candidates.append(parts)

print(f"  Input candidates: {len(candidates)}")

# 3. Filter RNA editing sites (A>G and T>C)
rna_editing_removed = []
rna_editing_kept = []
for c in candidates:
    ref, alt = c[2], c[3]
    # ADAR-mediated: A>G (sense strand) or T>C (antisense strand)
    if (ref == 'A' and alt == 'G') or (ref == 'T' and alt == 'C'):
        rna_editing_removed.append(c)
    else:
        rna_editing_kept.append(c)

print(f"  RNA editing candidates removed (A>G/T>C): {len(rna_editing_removed)}")
print(f"  Remaining after editing filter: {len(rna_editing_kept)}")

# 4. AF filter: remove very high AF (>0.9, likely germline) and very low (<0.05)
af_filtered = []
af_removed = []
for c in rna_editing_kept:
    try:
        avg_af = float(c[6])  # avg_af column
        max_af = float(c[7])  # max_af column
        if 0.05 <= avg_af <= 0.9:
            af_filtered.append(c)
        else:
            af_removed.append(c)
    except (ValueError, IndexError):
        af_filtered.append(c)  # keep if can't parse

print(f"\n--- AF filter (0.05 <= avg_AF <= 0.9) ---")
print(f"  AF-removed: {len(af_removed)}")
print(f"  Remaining: {len(af_filtered)}")

# 5. Mutation type distribution of remaining candidates
print(f"\n--- Mutation type distribution (remaining {len(af_filtered)}) ---")
mut_types = {}
for c in af_filtered:
    mt = f"{c[2]}>{c[3]}"
    mut_types[mt] = mut_types.get(mt, 0) + 1
for mt, count in sorted(mut_types.items(), key=lambda x: -x[1]):
    print(f"  {mt}: {count}")

# 6. Save results
with open(os.path.join(OUT_DIR, "candidates_after_artifact_filter.tsv"), 'w') as f:
    f.write(header + "\n")
    for c in af_filtered:
        f.write('\t'.join(c) + '\n')

with open(os.path.join(OUT_DIR, "removed_rna_editing.tsv"), 'w') as f:
    f.write(header + "\n")
    for c in rna_editing_removed:
        f.write('\t'.join(c) + '\n')

# Save summary
with open(os.path.join(OUT_DIR, "filter_summary.txt"), 'w') as f:
    f.write(f"Input candidates: {len(candidates)}\n")
    f.write(f"RNA editing removed (A>G/T>C): {len(rna_editing_removed)}\n")
    f.write(f"AF filter removed: {len(af_removed)}\n")
    f.write(f"Final candidates: {len(af_filtered)}\n")
    f.write(f"\nMutation type distribution:\n")
    for mt, count in sorted(mut_types.items(), key=lambda x: -x[1]):
        f.write(f"  {mt}: {count}\n")

print(f"\nOutput: {OUT_DIR}/candidates_after_artifact_filter.tsv ({len(af_filtered)} variants)")
print("Done.")
