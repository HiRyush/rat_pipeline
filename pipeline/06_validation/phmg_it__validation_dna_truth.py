#!/usr/bin/env python3
"""
Validation: DNA Ground Truth Check
- Check each pipeline candidate against matched DNA WGS
- Report: true somatic, germline leak, RNA-only artifact
- Also check broader candidate sets (each filter stage)
"""
import os
import subprocess

BCFTOOLS = os.path.expanduser("~/RAT_project/miniforge3/envs/rnaseq/bin/bcftools")
DNA_DIR = "/home/yusanghyeon/RAT_project/PHMG_IT/results/ground_truth"
OUT_DIR = "/home/yusanghyeon/RAT_project/pipeline_candidates/results/validation"
os.makedirs(OUT_DIR, exist_ok=True)

CONTROLS = ["C1", "C2", "C3", "C4", "C5"]
TREATMENTS = ["P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10"]

def check_dna(chrom, pos, ref, alt, sample):
    """Check if variant exists in DNA ground truth for given sample"""
    dna_vcf = os.path.join(DNA_DIR, f"{sample}_dna_snps.vcf.gz")
    if not os.path.exists(dna_vcf):
        return "no_dna"

    result = subprocess.run(
        [BCFTOOLS, "view", "-r", f"{chrom}:{pos}-{pos}", dna_vcf],
        capture_output=True, text=True
    )

    for line in result.stdout.strip().split('\n'):
        if line.startswith('#'):
            continue
        if not line.strip():
            continue
        fields = line.split('\t')
        if len(fields) >= 10:
            dna_pos = fields[1]
            dna_ref = fields[3]
            dna_alt = fields[4]
            gt_field = fields[9]
            gt = gt_field.split(':')[0]

            if dna_pos == pos:
                # Check if ALT matches
                if alt in dna_alt.split(','):
                    if gt in ['0/1', '1/1', '0|1', '1|1', '1/0', '1|0']:
                        return "dna_confirmed"
                    else:
                        return "dna_ref"  # DNA has the site but genotype is REF

    return "dna_absent"  # Position not in DNA VCF or no matching ALT

def validate_candidate_set(candidates_file, label, out_file):
    """Validate a set of candidates against DNA ground truth"""
    print(f"\n--- Validating: {label} ---")

    candidates = []
    with open(candidates_file) as f:
        header = f.readline().strip()
        for line in f:
            parts = line.strip().split('\t')
            candidates.append(parts)

    print(f"  Total candidates: {len(candidates)}")

    results = {
        'treatment_dna_only': [],  # In treatment DNA, NOT in control DNA → TRUE SOMATIC
        'control_dna_also': [],    # In control DNA too → germline leak
        'rna_only': [],            # Not in any DNA → RNA artifact or novel
        'mixed': []                # Partial matches
    }

    detailed = []

    for c in candidates:
        chrom, pos, ref, alt = c[0], c[1], c[2], c[3]
        samples_str = c[5] if len(c) > 5 else ""
        treatment_samples = samples_str.split(',') if samples_str else TREATMENTS

        # Check in Treatment DNA
        treatment_dna_hits = 0
        for s in TREATMENTS:
            status = check_dna(chrom, pos, ref, alt, s)
            if status == "dna_confirmed":
                treatment_dna_hits += 1

        # Check in Control DNA
        control_dna_hits = 0
        for s in CONTROLS:
            status = check_dna(chrom, pos, ref, alt, s)
            if status == "dna_confirmed":
                control_dna_hits += 1

        # Classify
        if treatment_dna_hits > 0 and control_dna_hits == 0:
            category = "TRUE_SOMATIC"
            results['treatment_dna_only'].append(c)
        elif treatment_dna_hits > 0 and control_dna_hits > 0:
            category = "GERMLINE_LEAK"
            results['control_dna_also'].append(c)
        elif treatment_dna_hits == 0:
            category = "RNA_ONLY"
            results['rna_only'].append(c)
        else:
            category = "MIXED"
            results['mixed'].append(c)

        detailed.append(c[:6] + [str(treatment_dna_hits), str(control_dna_hits), category])

    # Print summary
    total = len(candidates)
    true_som = len(results['treatment_dna_only'])
    germ = len(results['control_dna_also'])
    rna_only = len(results['rna_only'])

    print(f"  TRUE SOMATIC (Treatment DNA only): {true_som} ({100*true_som/total:.1f}%)" if total > 0 else "  No candidates")
    print(f"  GERMLINE LEAK (Control DNA also):  {germ} ({100*germ/total:.1f}%)" if total > 0 else "")
    print(f"  RNA ONLY (no DNA evidence):        {rna_only} ({100*rna_only/total:.1f}%)" if total > 0 else "")

    # Save detailed results
    with open(out_file, 'w') as f:
        f.write("chrom\tpos\tref\talt\trecurrence\tsamples\ttreatment_dna_hits\tcontrol_dna_hits\tcategory\n")
        for d in detailed:
            f.write('\t'.join(d) + '\n')

    return results

print("=" * 70)
print("DNA Ground Truth Validation")
print("=" * 70)

# Validate at each filter stage
stages = [
    ("/home/yusanghyeon/RAT_project/pipeline_candidates/results/differential/treatment_only_rec2.tsv",
     "Step 4: Treatment-only rec>=2 (117)",
     os.path.join(OUT_DIR, "validation_step4_rec2.tsv")),

    ("/home/yusanghyeon/RAT_project/pipeline_candidates/results/artifact_filter/candidates_after_artifact_filter.tsv",
     "Step 6: After artifact filter (79)",
     os.path.join(OUT_DIR, "validation_step6_artifact.tsv")),

    ("/home/yusanghyeon/RAT_project/pipeline_candidates/results/pop_gp_filter/candidates_after_pop_gp_filter.tsv",
     "Step 7: After pop+GP filter (5)",
     os.path.join(OUT_DIR, "validation_step7_final.tsv")),
]

all_results = {}
for fpath, label, out_file in stages:
    if os.path.exists(fpath):
        all_results[label] = validate_candidate_set(fpath, label, out_file)

# Summary table
print("\n" + "=" * 70)
print("VALIDATION SUMMARY")
print("=" * 70)
print(f"{'Stage':<45} {'Total':>6} {'Somatic':>8} {'Germline':>9} {'RNA-only':>9}")
print("-" * 70)
for label, res in all_results.items():
    total = sum(len(v) for v in res.values())
    som = len(res['treatment_dna_only'])
    germ = len(res['control_dna_also'])
    rna = len(res['rna_only'])
    print(f"{label:<45} {total:>6} {som:>8} {germ:>9} {rna:>9}")

# PPV calculation
print("\n--- Precision (PPV) ---")
for label, res in all_results.items():
    total = sum(len(v) for v in res.values())
    som = len(res['treatment_dna_only'])
    if total > 0:
        ppv = 100 * som / total
        print(f"  {label}: PPV = {ppv:.1f}% ({som}/{total})")

print(f"\nResults saved to: {OUT_DIR}")
print("Done.")
