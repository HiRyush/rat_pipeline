#!/usr/bin/env python3
"""
Step 4: Multi-sample Differential Filter
- Load all 15 samples' discordances
- Remove any discordance that appears in ANY Control sample
- Keep Treatment-only discordances with recurrence >= 2
"""
import os
import sys
from collections import defaultdict

DISC_DIR = "/home/yusanghyeon/RAT_project/pipeline_candidates/results/discordance"
OUT_DIR = "/home/yusanghyeon/RAT_project/pipeline_candidates/results/differential"
os.makedirs(OUT_DIR, exist_ok=True)

CONTROLS = ["C1", "C2", "C3", "C4", "C5"]
TREATMENTS = ["P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10"]

def load_discordances(sample):
    """Load discordance file: chrom, pos, ref, alt, obs_dp, obs_af, imp_gt, imp_ds"""
    variants = {}
    fpath = os.path.join(DISC_DIR, f"{sample}_discordances.tsv")
    with open(fpath) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                key = f"{parts[0]}:{parts[1]}:{parts[2]}:{parts[3]}"  # chrom:pos:ref:alt
                variants[key] = {
                    'chrom': parts[0],
                    'pos': parts[1],
                    'ref': parts[2],
                    'alt': parts[3],
                    'obs_dp': parts[4],
                    'obs_af': parts[5],
                    'imp_gt': parts[6] if len(parts) > 6 else '.',
                    'imp_ds': parts[7] if len(parts) > 7 else '.'
                }
    return variants

print("=" * 70)
print("Step 4: Multi-sample Differential Filter")
print("=" * 70)

# 1. Load all Control discordances → build Control set
print("\n--- Loading Control discordances ---")
control_keys = set()
for s in CONTROLS:
    disc = load_discordances(s)
    control_keys.update(disc.keys())
    print(f"  {s}: {len(disc)} discordances")
print(f"  Total unique Control discordance sites: {len(control_keys)}")

# 2. Load Treatment discordances
print("\n--- Loading Treatment discordances ---")
treatment_data = {}  # key -> list of (sample, data)
treatment_per_sample = {}
for s in TREATMENTS:
    disc = load_discordances(s)
    treatment_per_sample[s] = disc
    for key, data in disc.items():
        if key not in treatment_data:
            treatment_data[key] = []
        treatment_data[key].append((s, data))
    print(f"  {s}: {len(disc)} discordances")
print(f"  Total unique Treatment discordance sites: {len(treatment_data)}")

# 3. Filter: Treatment-only (not in any Control)
treatment_only = {k: v for k, v in treatment_data.items() if k not in control_keys}
print(f"\n--- After removing Control overlaps ---")
print(f"  Removed (in Control): {len(treatment_data) - len(treatment_only)}")
print(f"  Treatment-only sites: {len(treatment_only)}")

# 4. Recurrence filter
rec1 = {k: v for k, v in treatment_only.items() if len(v) >= 1}
rec2 = {k: v for k, v in treatment_only.items() if len(v) >= 2}
rec3 = {k: v for k, v in treatment_only.items() if len(v) >= 3}

print(f"\n--- Recurrence distribution ---")
print(f"  Recurrence >= 1: {len(rec1)}")
print(f"  Recurrence >= 2: {len(rec2)}")
print(f"  Recurrence >= 3: {len(rec3)}")

# 5. Output: recurrence >= 2 as main candidates
print(f"\n--- Writing output (recurrence >= 2) ---")
with open(os.path.join(OUT_DIR, "treatment_only_rec2.tsv"), 'w') as f:
    f.write("chrom\tpos\tref\talt\trecurrence\tsamples\tavg_af\tmax_af\n")
    for key in sorted(rec2.keys(), key=lambda x: (x.split(':')[0], int(x.split(':')[1]))):
        parts = key.split(':')
        samples_list = [s for s, d in rec2[key]]
        afs = []
        for s, d in rec2[key]:
            try:
                afs.append(float(d['obs_af']))
            except (ValueError, KeyError):
                pass
        avg_af = sum(afs) / len(afs) if afs else 0
        max_af = max(afs) if afs else 0
        f.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\t{parts[3]}\t"
                f"{len(samples_list)}\t{','.join(samples_list)}\t"
                f"{avg_af:.4f}\t{max_af:.4f}\n")

# Also output rec >= 1 for reference
with open(os.path.join(OUT_DIR, "treatment_only_rec1.tsv"), 'w') as f:
    f.write("chrom\tpos\tref\talt\trecurrence\tsamples\tavg_af\tmax_af\n")
    for key in sorted(rec1.keys(), key=lambda x: (x.split(':')[0], int(x.split(':')[1]))):
        parts = key.split(':')
        samples_list = [s for s, d in rec1[key]]
        afs = []
        for s, d in rec1[key]:
            try:
                afs.append(float(d['obs_af']))
            except (ValueError, KeyError):
                pass
        avg_af = sum(afs) / len(afs) if afs else 0
        max_af = max(afs) if afs else 0
        f.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\t{parts[3]}\t"
                f"{len(samples_list)}\t{','.join(samples_list)}\t"
                f"{avg_af:.4f}\t{max_af:.4f}\n")

# 6. Summary
print(f"\nOutput files:")
print(f"  treatment_only_rec1.tsv: {len(rec1)} variants")
print(f"  treatment_only_rec2.tsv: {len(rec2)} variants (main candidates)")
print(f"\nDone.")
