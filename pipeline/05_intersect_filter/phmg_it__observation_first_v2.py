#!/usr/bin/env python3
"""
Observation-first + Imputation-filter Pipeline v2
Uses pre-computed results:
  - MuTect2 Phase 9 differential: 65,375 SNPs
  - Imputed differential: 336,794 SNPs
  - Intersection → additional filters → DNA validation
"""
import os
import csv
import subprocess
import sys
from collections import defaultdict

# Force unbuffered output
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

BCFTOOLS = os.path.expanduser("~/RAT_project/miniforge3/envs/rnaseq/bin/bcftools")
BASE = "/home/yusanghyeon/RAT_project"
OUT_DIR = f"{BASE}/pipeline_candidates/results/observation_first"
os.makedirs(OUT_DIR, exist_ok=True)

CONTROLS = ["C1", "C2", "C3", "C4", "C5"]
TREATMENTS = ["P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10"]

def check_dna(chrom, pos, ref, alt, sample):
    dna_vcf = f"{BASE}/PHMG_IT/results/ground_truth/{sample}_dna_snps.vcf.gz"
    if not os.path.exists(dna_vcf):
        return "no_dna"
    result = subprocess.run(
        [BCFTOOLS, "view", "-r", f"{chrom}:{pos}-{pos}", dna_vcf],
        capture_output=True, text=True
    )
    for line in result.stdout.strip().split('\n'):
        if line.startswith('#') or not line.strip():
            continue
        fields = line.split('\t')
        if len(fields) >= 10 and fields[1] == str(pos):
            dna_alt = fields[4]
            gt = fields[9].split(':')[0]
            if alt in dna_alt.split(',') and gt in ['0/1','1/1','0|1','1|1','1/0','1|0']:
                return "dna_confirmed"
    return "dna_absent"

def validate_set(variants, label):
    true_somatic, germline_leak, rna_only = [], [], []
    total = len(variants)
    for i, v in enumerate(variants):
        if (i+1) % 100 == 0:
            print(f"    Validating {i+1}/{total}...")
        chrom, pos, ref, alt = v['chrom'], str(v['pos']), v['ref'], v['alt']
        treat_hits = sum(1 for s in TREATMENTS if check_dna(chrom, pos, ref, alt, s) == "dna_confirmed")
        ctrl_hits = sum(1 for s in CONTROLS if check_dna(chrom, pos, ref, alt, s) == "dna_confirmed")
        v['treat_dna'] = treat_hits
        v['ctrl_dna'] = ctrl_hits
        if treat_hits > 0 and ctrl_hits == 0:
            v['category'] = 'TRUE_SOMATIC'
            true_somatic.append(v)
        elif treat_hits > 0 and ctrl_hits > 0:
            v['category'] = 'GERMLINE_LEAK'
            germline_leak.append(v)
        else:
            v['category'] = 'RNA_ONLY'
            rna_only.append(v)

    print(f"\n  [{label}] Total: {total}")
    if total > 0:
        print(f"    TRUE SOMATIC:  {len(true_somatic)} ({100*len(true_somatic)/total:.1f}%)")
        print(f"    GERMLINE LEAK: {len(germline_leak)} ({100*len(germline_leak)/total:.1f}%)")
        print(f"    RNA ONLY:      {len(rna_only)} ({100*len(rna_only)/total:.1f}%)")
    return true_somatic, germline_leak, rna_only

print("=" * 70)
print("Observation-first + Imputation-filter Pipeline v2")
print("=" * 70)

# ======================================================================
# Step 2: Load MuTect2 Phase 9 differential
# ======================================================================
print("\n--- Step 2: Loading MuTect2 Phase 9 differential ---")
mutect2_keys = {}  # chrom:pos:alt -> variant data
with open(f"{BASE}/PHMG_IT/results/mutect2/differential/differential_snps.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        key = f"{row['chrom']}:{row['pos']}:{row['alt']}"
        mutect2_keys[key] = {
            'chrom': row['chrom'],
            'pos': row['pos'],
            'ref': row['ref'],
            'alt': row['alt'],
            'recurrence': int(row['recurrence']),
            'target_mean_af': float(row['target_mean_alt_ratio']),
            'fisher_pvalue': row['fisher_pvalue'],
        }
print(f"  MuTect2 differential SNPs: {len(mutect2_keys)}")

mutect2_rec2 = {k: v for k, v in mutect2_keys.items() if v['recurrence'] >= 2}
print(f"  After recurrence >= 2: {len(mutect2_rec2)}")

# ======================================================================
# Step 3+4: Load Imputed differential & Intersect
# ======================================================================
print("\n--- Step 4: Loading Imputed differential & Intersecting ---")
imputed_keys = set()
with open(f"{BASE}/PHMG_IT/results/imputed_differential/imputed_differential_snps.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        key = f"{row['chrom']}:{row['pos']}:{row['alt']}"
        imputed_keys.add(key)
print(f"  Imputed differential SNPs: {len(imputed_keys)}")

# Intersection
intersection_keys = set(mutect2_rec2.keys()) & imputed_keys
step4_variants = [mutect2_rec2[k] for k in intersection_keys]
step4_variants.sort(key=lambda v: (v['chrom'], int(v['pos'])))
print(f"  Intersection (MuTect2 rec>=2 ∩ Imputed): {len(step4_variants)}")

# Also check: what about MuTect2 ALL (no rec filter) ∩ Imputed?
all_intersection = set(mutect2_keys.keys()) & imputed_keys
print(f"  (Reference: MuTect2 ALL ∩ Imputed = {len(all_intersection)})")

# ======================================================================
# Step 5a: RNA editing filter
# ======================================================================
print("\n--- Step 5a: RNA editing filter (A>G / T>C) ---")
step5a = [v for v in step4_variants if not (
    (v['ref'] == 'A' and v['alt'] == 'G') or
    (v['ref'] == 'T' and v['alt'] == 'C')
)]
print(f"  Removed (A>G/T>C): {len(step4_variants) - len(step5a)}")
print(f"  Remaining: {len(step5a)}")

# ======================================================================
# Step 5b: Expression-aware filter
# ======================================================================
print("\n--- Step 5b: Expression-aware filter ---")
newly_expressed_file = f"{BASE}/pipeline_candidates/results/deseq2/newly_expressed_genes.txt"
newly_expressed = set()
if os.path.exists(newly_expressed_file):
    with open(newly_expressed_file) as f:
        newly_expressed = {l.strip() for l in f if l.strip()}
print(f"  Newly expressed genes: {len(newly_expressed)}")
# Gene-level filtering requires position→gene mapping; apply later
step5b = step5a
print(f"  (Positional gene mapping pending — {len(step5b)} kept)")

# ======================================================================
# Step 5c: Population AF filter (HRDP)
# ======================================================================
print("\n--- Step 5c: Population AF filter (HRDP) ---")
hrdp_variants = set()
for chri in range(1, 21):
    panel_vcf = f"{BASE}/PHMG_IT/imputation/ref_panel/phased/hrdp_chr{chri}_phased.vcf.gz"
    if os.path.exists(panel_vcf):
        result = subprocess.run(
            [BCFTOOLS, "query", "-f", "%CHROM:%POS:%ALT\n", panel_vcf],
            capture_output=True, text=True
        )
        for line in result.stdout.strip().split('\n'):
            if line:
                hrdp_variants.add(line)
print(f"  HRDP variants loaded: {len(hrdp_variants)}")

step5c = []
pop_removed = 0
for v in step5b:
    key = f"{v['chrom']}:{v['pos']}:{v['alt']}"
    if key in hrdp_variants:
        pop_removed += 1
    else:
        step5c.append(v)
print(f"  Removed (in HRDP): {pop_removed}")
print(f"  Remaining: {len(step5c)}")

# ======================================================================
# Step 5d: AF quality filter
# ======================================================================
print("\n--- Step 5d: AF quality filter (0.05-0.9) ---")
step5d = [v for v in step5c if 0.05 <= v['target_mean_af'] <= 0.9]
print(f"  Removed: {len(step5c) - len(step5d)}")
print(f"  Remaining: {len(step5d)}")

# ======================================================================
# Filter cascade summary
# ======================================================================
print("\n" + "=" * 70)
print("FILTER CASCADE SUMMARY")
print("=" * 70)
stages_summary = [
    ("MuTect2 differential (all)", len(mutect2_keys)),
    ("MuTect2 differential (rec>=2)", len(mutect2_rec2)),
    ("∩ Imputed differential", len(step4_variants)),
    ("- RNA editing (A>G/T>C)", len(step5a)),
    ("- HRDP population", len(step5c)),
    ("- AF filter (0.05-0.9)", len(step5d)),
]
for label, count in stages_summary:
    print(f"  {label:<40} {count:>8}")

# ======================================================================
# Step 6: DNA Ground Truth Validation
# ======================================================================
print("\n" + "=" * 70)
print("Step 6: DNA Ground Truth Validation")
print("=" * 70)

# Validate intersection (Step 4)
print("\n--- Validating Step 4: Intersection ---")
som4, germ4, rna4 = validate_set(step4_variants, "Step 4")

# Validate after RNA editing filter (Step 5a)
print("\n--- Validating Step 5a: After RNA editing filter ---")
som5a, germ5a, rna5a = validate_set(step5a, "Step 5a")

# Validate final (Step 5d)
print("\n--- Validating Step 5d: Final candidates ---")
som5d, germ5d, rna5d = validate_set(step5d, "Step 5d")

# ======================================================================
# Final results
# ======================================================================
print("\n" + "=" * 70)
print("FINAL RESULTS")
print("=" * 70)
print(f"\n{'Stage':<45} {'Total':>7} {'Somatic':>8} {'Germline':>9} {'RNA-only':>9} {'PPV':>7}")
print("-" * 87)
for label, som, germ, rna in [
    ("Step 4: MuTect2∩Imputed", som4, germ4, rna4),
    ("Step 5a: - RNA editing", som5a, germ5a, rna5a),
    ("Step 5d: Final (all filters)", som5d, germ5d, rna5d),
]:
    total = len(som) + len(germ) + len(rna)
    ppv = f"{100*len(som)/total:.1f}%" if total > 0 else "N/A"
    print(f"  {label:<43} {total:>7} {len(som):>8} {len(germ):>9} {len(rna):>9} {ppv:>7}")

# Save results
def save_tsv(variants, fpath):
    if not variants:
        with open(fpath, 'w') as f:
            f.write("(empty)\n")
        return
    with open(fpath, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=variants[0].keys(), delimiter='\t')
        w.writeheader()
        w.writerows(variants)

save_tsv(step4_variants, f"{OUT_DIR}/step4_intersection.tsv")
save_tsv(step5d, f"{OUT_DIR}/step5d_final_candidates.tsv")
save_tsv(som5d, f"{OUT_DIR}/true_somatic_final.tsv")
save_tsv(germ5d, f"{OUT_DIR}/germline_leak_final.tsv")
save_tsv(rna5d, f"{OUT_DIR}/rna_only_final.tsv")

# Mutation types of true somatic
if som5d:
    print(f"\n--- True Somatic Mutation Types ({len(som5d)}) ---")
    mt = defaultdict(int)
    for v in som5d:
        mt[f"{v['ref']}>{v['alt']}"] += 1
    for t, c in sorted(mt.items(), key=lambda x: -x[1]):
        print(f"  {t}: {c}")

# Comparison with Phase 11
print(f"\n--- Comparison with Phase 11 ---")
print(f"  Phase 11 (Phase7 ∩ Imputed): 977 → 597 true somatic (61.1% PPV)")
print(f"  This pipeline (MuTect2 ∩ Imputed + filters): {len(step5d)} → {len(som5d)} true somatic", end="")
if len(step5d) > 0:
    print(f" ({100*len(som5d)/len(step5d):.1f}% PPV)")
else:
    print(" (N/A)")

print(f"\nResults saved to: {OUT_DIR}")
print("Done.")
