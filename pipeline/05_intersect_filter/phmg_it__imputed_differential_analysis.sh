#!/bin/bash
# =============================================================================
# Imputed VCF Differential Analysis: Control vs Treatment
# Phase 7/9와 비교하기 위한 imputation 후 differential variant 분석
# =============================================================================
set -euo pipefail

BCFTOOLS="/home/yusanghyeon/RAT_project/miniforge3/envs/rnaseq/bin/bcftools"
IMPUTED_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/imputation/imputed"
OUTDIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/imputed_differential"
mkdir -p "$OUTDIR"

CONTROLS=(C1 C2 C3 C4 C5)
TREATMENTS=(P1 P2 P3 P4 P5 P6 P7 P8 P9 P10)

echo "============================================"
echo "Imputed VCF Differential Analysis"
echo "Started: $(date)"
echo "============================================"

# ---------------------------------------------------------
# Step 1: Merge all 15 imputed VCFs
# ---------------------------------------------------------
MERGED="$OUTDIR/all15_imputed_merged.vcf.gz"
if [ ! -f "$MERGED" ]; then
    echo "[Step 1] Merging 15 imputed VCFs..."

    # Create file list
    VCFLIST="$OUTDIR/vcf_list.txt"
    > "$VCFLIST"
    for s in "${CONTROLS[@]}" "${TREATMENTS[@]}"; do
        echo "$IMPUTED_DIR/${s}_imputed.vcf.gz" >> "$VCFLIST"
    done

    $BCFTOOLS merge -l "$VCFLIST" -Oz -o "$MERGED" --threads 4
    $BCFTOOLS index -t "$MERGED"
    echo "[Step 1] Done. $(date)"
else
    echo "[Step 1] Merged VCF already exists, skipping."
fi

# ---------------------------------------------------------
# Step 2: Extract differential variants using GT comparison
# ---------------------------------------------------------
echo "[Step 2] Extracting differential variants (GT-based)..."

python3 << 'PYEOF'
import gzip
import csv
import sys
from collections import defaultdict

MERGED = "/home/yusanghyeon/RAT_project/PHMG_IT/results/imputed_differential/all15_imputed_merged.vcf.gz"
OUTDIR = "/home/yusanghyeon/RAT_project/PHMG_IT/results/imputed_differential"

CONTROLS = ["C1", "C2", "C3", "C4", "C5"]
TREATMENTS = ["P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10"]

# Parameters matching Phase 7 criteria
MIN_RECURRENCE = 2  # min treatment samples with ALT

print("Parsing merged imputed VCF...")
sample_names = []
ctrl_idx = []
treat_idx = []

diff_snps = []
diff_indels = []
total_sites = 0
sites_with_treat_alt = 0

with gzip.open(MERGED, 'rt') as f:
    for line in f:
        if line.startswith('##'):
            continue
        if line.startswith('#CHROM'):
            fields = line.strip().split('\t')
            sample_names = fields[9:]
            for i, s in enumerate(sample_names):
                if s in CONTROLS:
                    ctrl_idx.append(i)
                elif s in TREATMENTS:
                    treat_idx.append(i)
            print(f"Samples: {len(sample_names)} (Control: {len(ctrl_idx)}, Treatment: {len(treat_idx)})")
            continue

        total_sites += 1
        if total_sites % 1000000 == 0:
            print(f"  Processed {total_sites:,} sites... (SNPs: {len(diff_snps):,}, INDELs: {len(diff_indels):,})")

        fields = line.strip().split('\t')
        chrom, pos, _, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
        info = fields[7]

        # Skip multi-allelic for simplicity
        if ',' in alt:
            continue

        # Parse FORMAT
        fmt = fields[8].split(':')
        gt_idx = fmt.index('GT') if 'GT' in fmt else 0
        ds_idx = fmt.index('DS') if 'DS' in fmt else -1

        # Extract genotypes
        samples = fields[9:]

        # Check control samples — any ALT?
        ctrl_has_alt = 0
        for i in ctrl_idx:
            gt = samples[i].split(':')[gt_idx]
            alleles = gt.replace('|', '/').split('/')
            if '1' in alleles:
                ctrl_has_alt += 1

        # Check treatment samples — how many have ALT?
        treat_with_alt = 0
        treat_dosages = []
        for i in treat_idx:
            sample_fields = samples[i].split(':')
            gt = sample_fields[gt_idx]
            alleles = gt.replace('|', '/').split('/')
            if '1' in alleles:
                treat_with_alt += 1
            if ds_idx >= 0 and len(sample_fields) > ds_idx:
                try:
                    treat_dosages.append(float(sample_fields[ds_idx]))
                except ValueError:
                    treat_dosages.append(0.0)

        if treat_with_alt == 0:
            continue

        sites_with_treat_alt += 1

        # Differential: Treatment has ALT, Control does NOT
        if ctrl_has_alt == 0 and treat_with_alt >= MIN_RECURRENCE:
            # Determine variant type
            is_snp = len(ref) == 1 and len(alt) == 1
            mean_dosage = sum(treat_dosages) / len(treat_dosages) if treat_dosages else 0

            record = {
                'chrom': chrom,
                'pos': int(pos),
                'ref': ref,
                'alt': alt,
                'recurrence': treat_with_alt,
                'ctrl_alt_count': ctrl_has_alt,
                'mean_dosage': round(mean_dosage, 4),
                'variant_type': 'SNP' if is_snp else 'INDEL'
            }

            if is_snp:
                diff_snps.append(record)
            else:
                diff_indels.append(record)

print(f"\nTotal sites: {total_sites:,}")
print(f"Sites with treatment ALT: {sites_with_treat_alt:,}")
print(f"Differential SNPs (ctrl=0, treat>={MIN_RECURRENCE}): {len(diff_snps):,}")
print(f"Differential INDELs (ctrl=0, treat>={MIN_RECURRENCE}): {len(diff_indels):,}")

# Recurrence distribution
print("\n=== Recurrence Distribution (SNPs) ===")
rec_dist = defaultdict(int)
for r in diff_snps:
    rec_dist[r['recurrence']] += 1
for k in sorted(rec_dist.keys(), reverse=True):
    print(f"  {k}/10: {rec_dist[k]:,}")

# Write results
for label, data in [('snps', diff_snps), ('indels', diff_indels)]:
    outpath = f"{OUTDIR}/imputed_differential_{label}.csv"
    if data:
        with open(outpath, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=data[0].keys())
            writer.writeheader()
            writer.writerows(sorted(data, key=lambda x: (-x['recurrence'], x['chrom'], x['pos'])))
    print(f"Written: {outpath} ({len(data):,} records)")

# Write summary
with open(f"{OUTDIR}/imputed_differential_summary.txt", 'w') as f:
    f.write("Imputed VCF Differential Analysis Summary\n")
    f.write(f"Date: $(date)\n")
    f.write(f"Total sites in merged VCF: {total_sites:,}\n")
    f.write(f"Sites with treatment ALT: {sites_with_treat_alt:,}\n")
    f.write(f"Differential SNPs (recurrence>={MIN_RECURRENCE}, ctrl=0): {len(diff_snps):,}\n")
    f.write(f"Differential INDELs (recurrence>={MIN_RECURRENCE}, ctrl=0): {len(diff_indels):,}\n")
    f.write(f"\nRecurrence distribution (SNPs):\n")
    for k in sorted(rec_dist.keys(), reverse=True):
        f.write(f"  {k}/10: {rec_dist[k]:,}\n")

print("\nDone!")
PYEOF

echo "[Step 2] Done. $(date)"
echo "============================================"
echo "Results in: $OUTDIR"
