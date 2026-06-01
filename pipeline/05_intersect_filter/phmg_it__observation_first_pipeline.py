#!/usr/bin/env python3
"""
Observation-first + Imputation-filter Pipeline
Step 2: Load MuTect2 Phase 9 differential (Treatment-only observed)
Step 4: Intersection with Imputation-based differential
Step 5: Additional filters (RNA editing, DESeq2, Population AF, GP, DP/AF)
Step 6: DNA Ground Truth Validation at each stage
"""
import os
import csv
import subprocess
from collections import defaultdict

BCFTOOLS = os.path.expanduser("~/RAT_project/miniforge3/envs/rnaseq/bin/bcftools")
BASE = "/home/yusanghyeon/RAT_project"
OUT_DIR = f"{BASE}/pipeline_candidates/results/observation_first"
os.makedirs(OUT_DIR, exist_ok=True)

CONTROLS = ["C1", "C2", "C3", "C4", "C5"]
TREATMENTS = ["P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10"]

# ======================================================================
# Helper: DNA ground truth validation
# ======================================================================
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
            if alt in dna_alt.split(',') and gt in ['0/1', '1/1', '0|1', '1|1', '1/0', '1|0']:
                return "dna_confirmed"
    return "dna_absent"

def validate_set(variants, label):
    """Validate a variant set against DNA. Returns categorized results."""
    true_somatic = []
    germline_leak = []
    rna_only = []

    for v in variants:
        chrom, pos, ref, alt = v['chrom'], v['pos'], v['ref'], v['alt']

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

    total = len(variants)
    print(f"\n  [{label}] Total: {total}")
    print(f"    TRUE SOMATIC:  {len(true_somatic)} ({100*len(true_somatic)/total:.1f}%)" if total else "")
    print(f"    GERMLINE LEAK: {len(germline_leak)} ({100*len(germline_leak)/total:.1f}%)" if total else "")
    print(f"    RNA ONLY:      {len(rna_only)} ({100*len(rna_only)/total:.1f}%)" if total else "")

    return true_somatic, germline_leak, rna_only

def save_variants(variants, filepath):
    if not variants:
        open(filepath, 'w').close()
        return
    with open(filepath, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=variants[0].keys(), delimiter='\t')
        writer.writeheader()
        writer.writerows(variants)

# ======================================================================
print("=" * 70)
print("Observation-first + Imputation-filter Pipeline")
print("=" * 70)

# ======================================================================
# Step 2: Load MuTect2 Phase 9 differential (Treatment-only)
# ======================================================================
print("\n--- Step 2: Loading MuTect2 differential variants ---")
diff_file = f"{BASE}/PHMG_IT/results/mutect2/differential/differential_snps.csv"
step2_variants = []
with open(diff_file) as f:
    reader = csv.DictReader(f)
    for row in reader:
        step2_variants.append({
            'chrom': row['chrom'],
            'pos': row['pos'],
            'ref': row['ref'],
            'alt': row['alt'],
            'recurrence': int(row['recurrence']),
            'target_mean_af': float(row['target_mean_alt_ratio']),
            'fisher_pvalue': row['fisher_pvalue'],
            'target_samples': row['target_samples_found'],
        })

print(f"  MuTect2 differential SNPs loaded: {len(step2_variants)}")

# Filter: recurrence >= 2
step2_rec2 = [v for v in step2_variants if v['recurrence'] >= 2]
print(f"  After recurrence >= 2: {len(step2_rec2)}")

# ======================================================================
# Step 4: Intersection with Imputed differential
# ======================================================================
print("\n--- Step 4: Building imputation-based differential set ---")

# Build per-sample imputed genotype lookup for Treatment and Control
# For each Treatment sample: collect sites where imputed genotype has ALT (0|1 or 1|1)
# For each Control sample: same
# Then find: Treatment-ALT but Control-REF in imputed

treatment_imputed = defaultdict(set)  # chrom:pos:alt -> set of samples
control_imputed = defaultdict(set)

for group, samples in [("treatment", TREATMENTS), ("control", CONTROLS)]:
    for s in samples:
        for chri in range(1, 21):
            imp_vcf = f"{BASE}/PHMG_IT/imputation/imputed/{s}_chr{chri}.vcf.gz"
            if not os.path.exists(imp_vcf):
                continue
            result = subprocess.run(
                [BCFTOOLS, "query", "-f", "%CHROM\t%POS\t%ALT\t[%GT]\n", imp_vcf],
                capture_output=True, text=True
            )
            for line in result.stdout.strip().split('\n'):
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) >= 4:
                    gt = parts[3]
                    if gt in ['0|1', '1|0', '1|1']:
                        key = f"{parts[0]}:{parts[1]}:{parts[2]}"
                        if group == "treatment":
                            treatment_imputed[key].add(s)
                        else:
                            control_imputed[key].add(s)
        print(f"    {s} imputed ALT loaded")

# Find imputed-differential: in Treatment but not in Control
imputed_diff_keys = set()
for key in treatment_imputed:
    if key not in control_imputed:
        imputed_diff_keys.add(key)
    elif len(treatment_imputed[key]) >= 2 and len(control_imputed[key]) == 0:
        imputed_diff_keys.add(key)

print(f"  Imputed Treatment-only ALT sites: {len(imputed_diff_keys)}")

# Intersection: Phase 9 observed ∩ Imputed differential
step4_variants = []
for v in step2_rec2:
    key = f"{v['chrom']}:{v['pos']}:{v['alt']}"
    if key in imputed_diff_keys:
        step4_variants.append(v)

print(f"  Intersection (Phase 9 ∩ Imputed): {len(step4_variants)}")

# ======================================================================
# Step 5a: RNA editing filter
# ======================================================================
print("\n--- Step 5a: RNA editing filter (A>G / T>C) ---")
step5a = [v for v in step4_variants if not (
    (v['ref'] == 'A' and v['alt'] == 'G') or
    (v['ref'] == 'T' and v['alt'] == 'C')
)]
removed_editing = len(step4_variants) - len(step5a)
print(f"  Removed (A>G/T>C): {removed_editing}")
print(f"  Remaining: {len(step5a)}")

# ======================================================================
# Step 5b: DESeq2 expression-aware filter
# ======================================================================
print("\n--- Step 5b: Expression-aware filter ---")
newly_expressed_file = f"{BASE}/pipeline_candidates/results/deseq2/newly_expressed_genes.txt"
newly_expressed_genes = set()
if os.path.exists(newly_expressed_file):
    with open(newly_expressed_file) as f:
        newly_expressed_genes = {line.strip() for line in f if line.strip()}
print(f"  Newly expressed genes in Treatment: {len(newly_expressed_genes)}")

# Map variant positions to genes using GTF
# For simplicity, we'll use bedtools-like approach via bcftools + GTF
# Skip this step if gene mapping is complex — flag for later
# For now, keep all (expression filter noted but not applied positionally)
step5b = step5a  # TODO: apply gene-level filter when gene mapping available
print(f"  (Gene-level mapping pending — all {len(step5b)} kept)")

# ======================================================================
# Step 5c: Population AF filter (HRDP)
# ======================================================================
print("\n--- Step 5c: Population AF filter (HRDP) ---")
# Load HRDP variants
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
print(f"  HRDP panel variants: {len(hrdp_variants)}")

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
print("\n--- Step 5d: AF quality filter ---")
step5d = [v for v in step5c if 0.05 <= v['target_mean_af'] <= 0.9]
print(f"  Removed (AF out of range): {len(step5c) - len(step5d)}")
print(f"  Remaining: {len(step5d)}")

# ======================================================================
# Summary of filter cascade
# ======================================================================
print("\n" + "=" * 70)
print("FILTER CASCADE SUMMARY")
print("=" * 70)
stages = {
    'Step 2: MuTect2 differential (rec>=2)': step2_rec2,
    'Step 4: ∩ Imputed differential': step4_variants,
    'Step 5a: RNA editing removed': step5a,
    'Step 5c: Population AF removed': step5c,
    'Step 5d: AF quality filter': step5d,
}
for label, vlist in stages.items():
    print(f"  {label}: {len(vlist)}")

# ======================================================================
# Step 6: DNA Ground Truth Validation (at key stages)
# ======================================================================
print("\n" + "=" * 70)
print("Step 6: DNA Ground Truth Validation")
print("=" * 70)

# Validate Step 4 (intersection)
print("\n--- Validating Step 4 (intersection) ---")
# Sample subset for speed if too many
if len(step4_variants) > 2000:
    print(f"  Too many variants ({len(step4_variants)}), validating first 2000...")
    validation_set = step4_variants[:2000]
else:
    validation_set = step4_variants

som4, germ4, rna4 = validate_set(validation_set, "Step 4: Intersection")
save_variants(validation_set, f"{OUT_DIR}/validation_step4_intersection.tsv")

# Validate after all filters (Step 5d)
print("\n--- Validating Step 5d (final candidates) ---")
som5, germ5, rna5 = validate_set(step5d, "Step 5d: Final")
save_variants(step5d, f"{OUT_DIR}/validation_step5d_final.tsv")

# ======================================================================
# Final Summary
# ======================================================================
print("\n" + "=" * 70)
print("FINAL RESULTS")
print("=" * 70)

print(f"\n{'Stage':<45} {'Total':>7} {'Somatic':>8} {'Germline':>9} {'RNA-only':>9} {'PPV':>6}")
print("-" * 85)

for label, som, germ, rna in [
    ("Step 4: Intersection", som4, germ4, rna4),
    ("Step 5d: Final (all filters)", som5, germ5, rna5),
]:
    total = len(som) + len(germ) + len(rna)
    ppv = f"{100*len(som)/total:.1f}%" if total > 0 else "N/A"
    print(f"  {label:<43} {total:>7} {len(som):>8} {len(germ):>9} {len(rna):>9} {ppv:>6}")

# Save final somatic candidates
save_variants(som5, f"{OUT_DIR}/true_somatic_final.tsv")
save_variants(germ5, f"{OUT_DIR}/germline_leak_final.tsv")
save_variants(rna5, f"{OUT_DIR}/rna_only_final.tsv")

# Mutation type distribution of true somatic
if som5:
    print(f"\n--- True Somatic Mutation Types ({len(som5)}) ---")
    mt_counts = defaultdict(int)
    for v in som5:
        if len(v['ref']) == 1 and len(v['alt']) == 1:
            mt_counts[f"{v['ref']}>{v['alt']}"] += 1
        else:
            mt_counts["INDEL"] += 1
    for mt, c in sorted(mt_counts.items(), key=lambda x: -x[1]):
        print(f"  {mt}: {c}")

print(f"\nResults saved to: {OUT_DIR}")
print("Done.")
