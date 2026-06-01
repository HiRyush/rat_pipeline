#!/usr/bin/env python3
"""
Step 7: Population AF + GP Confidence Filter
- Check if variant exists in HRDP population (germline → remove)
- Check imputation GP score for confidence
"""
import os
import subprocess

ARTIFACT_DIR = "/home/yusanghyeon/RAT_project/pipeline_candidates/results/artifact_filter"
IMPUTED_DIR = "/home/yusanghyeon/RAT_project/PHMG_IT/imputation/imputed"
PANEL_DIR = "/home/yusanghyeon/RAT_project/PHMG_IT/imputation/ref_panel/phased"
OUT_DIR = "/home/yusanghyeon/RAT_project/pipeline_candidates/results/pop_gp_filter"
os.makedirs(OUT_DIR, exist_ok=True)

BCFTOOLS = os.path.expanduser("~/RAT_project/miniforge3/envs/rnaseq/bin/bcftools")

print("=" * 70)
print("Step 7: Population AF + GP Confidence Filter")
print("=" * 70)

# 1. Load candidates from Step 6
candidates = []
with open(os.path.join(ARTIFACT_DIR, "candidates_after_artifact_filter.tsv")) as f:
    header = f.readline().strip()
    for line in f:
        parts = line.strip().split('\t')
        candidates.append(parts)

print(f"\nInput candidates: {len(candidates)}")

# 2. Build HRDP population variant set (all known germline positions)
print("\n--- Building HRDP population AF set ---")
hrdp_variants = set()
for chri in range(1, 21):
    panel_vcf = os.path.join(PANEL_DIR, f"hrdp_chr{chri}_phased.vcf.gz")
    if os.path.exists(panel_vcf):
        result = subprocess.run(
            [BCFTOOLS, "query", "-f", "%CHROM:%POS:%REF:%ALT\n", panel_vcf],
            capture_output=True, text=True
        )
        for line in result.stdout.strip().split('\n'):
            if line:
                hrdp_variants.add(line)

print(f"  HRDP panel variants loaded: {len(hrdp_variants)}")

# 3. Filter: remove variants found in HRDP (= known germline in rat population)
pop_removed = []
pop_kept = []
for c in candidates:
    key = f"{c[0]}:{c[1]}:{c[2]}:{c[3]}"
    if key in hrdp_variants:
        pop_removed.append(c)
    else:
        pop_kept.append(c)

print(f"  Removed (found in HRDP population): {len(pop_removed)}")
print(f"  Remaining: {len(pop_kept)}")

# 4. GP confidence filter
# Check imputation dosage score (DS) for each candidate
# DS close to 0 for 0|0 = high confidence REF → discordance is meaningful
# DS close to 1 = imputation was uncertain → might be imputation error
print("\n--- GP/DS confidence check ---")
# For candidates that pass population filter, check DS from imputed VCFs
# We need to check DS at each position in the Treatment samples that have the discordance

gp_kept = []
gp_removed = []
for c in pop_kept:
    chrom = c[0]
    pos = c[1]
    samples = c[5].split(',')
    chr_num = chrom.replace('chr', '')

    # Check DS in imputed VCFs for the first sample
    ds_values = []
    for s in samples[:3]:  # Check up to 3 samples for efficiency
        imp_vcf = os.path.join(IMPUTED_DIR, f"{s}_chr{chr_num}.vcf.gz")
        if os.path.exists(imp_vcf):
            result = subprocess.run(
                [BCFTOOLS, "query", "-r", f"{chrom}:{pos}-{pos}",
                 "-f", "[%DS]\n", imp_vcf],
                capture_output=True, text=True
            )
            if result.stdout.strip():
                try:
                    ds = float(result.stdout.strip().split('\n')[0])
                    ds_values.append(ds)
                except ValueError:
                    pass

    # DS close to 0 = confident REF imputation → keep
    # DS > 0.5 = uncertain → imputation might be wrong → remove
    if ds_values:
        avg_ds = sum(ds_values) / len(ds_values)
        if avg_ds <= 0.3:  # Confident REF imputation
            gp_kept.append(c + [f"{avg_ds:.3f}"])
        else:
            gp_removed.append(c + [f"{avg_ds:.3f}"])
    else:
        # No DS available — keep conservatively
        gp_kept.append(c + ["NA"])

print(f"  GP/DS removed (DS > 0.3, uncertain imputation): {len(gp_removed)}")
print(f"  Final candidates: {len(gp_kept)}")

# 5. Save results
with open(os.path.join(OUT_DIR, "candidates_after_pop_gp_filter.tsv"), 'w') as f:
    f.write(header + "\tavg_imp_ds\n")
    for c in gp_kept:
        f.write('\t'.join(c) + '\n')

with open(os.path.join(OUT_DIR, "removed_population.tsv"), 'w') as f:
    f.write(header + "\n")
    for c in pop_removed:
        f.write('\t'.join(c) + '\n')

with open(os.path.join(OUT_DIR, "filter_summary.txt"), 'w') as f:
    f.write(f"Input from artifact filter: {len(candidates)}\n")
    f.write(f"HRDP population removed: {len(pop_removed)}\n")
    f.write(f"GP/DS confidence removed: {len(gp_removed)}\n")
    f.write(f"Final candidates: {len(gp_kept)}\n")

# Mutation type distribution of final
print(f"\n--- Final mutation type distribution ({len(gp_kept)} variants) ---")
mut_types = {}
for c in gp_kept:
    ref, alt = c[2], c[3]
    if len(ref) == 1 and len(alt) == 1:
        mt = f"{ref}>{alt}"
    else:
        mt = "INDEL"
    mut_types[mt] = mut_types.get(mt, 0) + 1
for mt, count in sorted(mut_types.items(), key=lambda x: -x[1]):
    print(f"  {mt}: {count}")

print(f"\nOutput: {OUT_DIR}/candidates_after_pop_gp_filter.tsv")
print("Done.")
