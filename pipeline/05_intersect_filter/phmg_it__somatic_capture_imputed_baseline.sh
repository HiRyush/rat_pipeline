#!/bin/bash
# =============================================================================
# 방안 2: Imputed Haplotype 기반 Somatic Mutation Capture
#
# 원리: Imputed VCF(germline baseline) vs RNA-seq observed variant 비교
#   - RNA-seq에서 ALT 관찰 + Imputed에서 REF(0|0) = somatic candidate
#   - Imputed에서도 ALT = germline variant (제외)
#
# 입력: per-sample RNA-seq called VCF + per-sample imputed VCF
# 출력: per-sample somatic candidates + Control vs Treatment differential
# =============================================================================
set -euo pipefail

BCFTOOLS="/home/yusanghyeon/RAT_project/miniforge3/envs/rnaseq/bin/bcftools"
CALLED_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/differential/per_sample_vcf"
IMPUTED_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/imputation/imputed"
OUTDIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/somatic_capture"
mkdir -p "$OUTDIR/per_sample"

CONTROLS=(C1 C2 C3 C4 C5)
TREATMENTS=(P1 P2 P3 P4 P5 P6 P7 P8 P9 P10)
ALL_SAMPLES=("${CONTROLS[@]}" "${TREATMENTS[@]}")

# Minimum quality filters for RNA-seq calls
MIN_QUAL=5
MIN_DP=5
MIN_ALT_RATIO=0.15

echo "============================================"
echo "방안 2: Somatic Capture (Imputed Baseline)"
echo "Started: $(date)"
echo "============================================"

# ---------------------------------------------------------
# Step 1: Per-sample somatic extraction
# ---------------------------------------------------------
echo "[Step 1] Per-sample somatic candidate extraction..."

for SAMPLE in "${ALL_SAMPLES[@]}"; do
    CALLED="$CALLED_DIR/${SAMPLE}_calls.vcf.gz"
    IMPUTED="$IMPUTED_DIR/${SAMPLE}_imputed.vcf.gz"
    OUT="$OUTDIR/per_sample/${SAMPLE}_somatic.vcf.gz"

    if [ -f "$OUT" ]; then
        echo "  [SKIP] $SAMPLE — already exists"
        continue
    fi

    echo "  [RUN] $SAMPLE..."

    python3 << PYEOF
import gzip
import sys

CALLED = "$CALLED"
IMPUTED = "$IMPUTED"
OUT = "$OUTDIR/per_sample/${SAMPLE}_somatic.vcf.gz"
SAMPLE = "$SAMPLE"
MIN_QUAL = $MIN_QUAL
MIN_DP = $MIN_DP
MIN_ALT_RATIO = $MIN_ALT_RATIO

# Step A: Load imputed genotypes (germline baseline)
imputed_gt = {}  # "chrom_pos_ref_alt" -> GT
with gzip.open(IMPUTED, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        chrom, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]
        gt_field = fields[9].split(':')[0]
        alleles = gt_field.replace('|', '/').split('/')
        has_alt = '1' in alleles
        key = f"{chrom}_{pos}_{ref}_{alt}"
        imputed_gt[key] = has_alt

# Step B: Scan RNA-seq called VCF, find somatic candidates
somatic = []
germline_match = 0
total_called = 0
low_quality = 0

with gzip.open(CALLED, 'rt') as f:
    header_lines = []
    for line in f:
        if line.startswith('#'):
            header_lines.append(line)
            continue

        fields = line.strip().split('\t')
        chrom, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]

        # Quality filter
        try:
            qual = float(fields[5]) if fields[5] != '.' else 0.0
        except ValueError:
            qual = 0.0

        # Parse DP and AD from FORMAT
        info = fields[7]
        fmt = fields[8].split(':')
        sample_data = fields[9].split(':')

        dp = 0
        alt_ratio = 0.0
        try:
            dp_idx = fmt.index('DP')
            dp = int(sample_data[dp_idx])
        except (ValueError, IndexError):
            pass

        try:
            ad_idx = fmt.index('AD')
            ad_vals = sample_data[ad_idx].split(',')
            ref_count = int(ad_vals[0])
            alt_count = int(ad_vals[1]) if len(ad_vals) > 1 else 0
            if ref_count + alt_count > 0:
                alt_ratio = alt_count / (ref_count + alt_count)
        except (ValueError, IndexError):
            pass

        # Apply quality filters
        if qual < MIN_QUAL and fields[5] != '.':
            low_quality += 1
            continue
        if dp < MIN_DP:
            low_quality += 1
            continue
        if alt_ratio < MIN_ALT_RATIO:
            low_quality += 1
            continue

        total_called += 1

        # Check each ALT allele
        for a in alt.split(','):
            key = f"{chrom}_{pos}_{ref}_{a}"

            if key in imputed_gt:
                if imputed_gt[key]:
                    # Imputed says ALT too -> germline
                    germline_match += 1
                else:
                    # Imputed says REF, RNA-seq says ALT -> SOMATIC CANDIDATE
                    somatic.append(line)
            else:
                # Not in imputed panel -> could be novel somatic
                somatic.append(line)

# Write somatic VCF
import gzip as gz
with gz.open(OUT, 'wt') as f:
    for h in header_lines:
        f.write(h)
    for s in somatic:
        f.write(s)

print(f"    {SAMPLE}: called={total_called:,} germline={germline_match:,} somatic={len(somatic):,} filtered={low_quality:,}")
PYEOF

    # Index the output
    $BCFTOOLS index -t "$OUTDIR/per_sample/${SAMPLE}_somatic.vcf.gz" 2>/dev/null || true
done

echo "[Step 1] Done. $(date)"

# ---------------------------------------------------------
# Step 2: Summary statistics
# ---------------------------------------------------------
echo ""
echo "[Step 2] Summary..."

python3 << 'PYEOF'
import gzip
import os
from collections import defaultdict

OUTDIR = "/home/yusanghyeon/RAT_project/PHMG_IT/results/somatic_capture"
CONTROLS = ["C1", "C2", "C3", "C4", "C5"]
TREATMENTS = ["P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10"]

print("=" * 60)
print("Per-sample Somatic Candidate Summary")
print("=" * 60)
print(f"\n{'Sample':<8} {'Group':<10} {'Somatic Candidates':>20}")
print("-" * 40)

sample_counts = {}
sample_variants = {}  # sample -> set of "chrom_pos_ref_alt"

for sample in CONTROLS + TREATMENTS:
    group = "Control" if sample in CONTROLS else "PHMG"
    vcf_path = f"{OUTDIR}/per_sample/{sample}_somatic.vcf.gz"

    if not os.path.exists(vcf_path):
        print(f"  {sample:<8} {group:<10} {'MISSING':>20}")
        continue

    variants = set()
    with gzip.open(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            for alt in fields[4].split(','):
                variants.add(f"{fields[0]}_{fields[1]}_{fields[3]}_{alt}")

    sample_counts[sample] = len(variants)
    sample_variants[sample] = variants
    print(f"  {sample:<8} {group:<10} {len(variants):>20,}")

# Group averages
ctrl_avg = sum(sample_counts.get(s, 0) for s in CONTROLS) / len(CONTROLS)
treat_avg = sum(sample_counts.get(s, 0) for s in TREATMENTS) / len(TREATMENTS)
print(f"\n  {'Control avg':<18} {ctrl_avg:>20,.0f}")
print(f"  {'Treatment avg':<18} {treat_avg:>20,.0f}")

# Differential: Treatment-specific somatic
# Find variants in >=2 treatment samples AND 0 control samples
print(f"\n{'='*60}")
print("Differential Somatic Analysis (Treatment-specific)")
print(f"{'='*60}")

all_positions = set()
for s in TREATMENTS:
    all_positions |= sample_variants.get(s, set())

ctrl_positions = set()
for s in CONTROLS:
    ctrl_positions |= sample_variants.get(s, set())

diff_somatic = defaultdict(int)  # variant -> treatment recurrence
for var in all_positions:
    if var in ctrl_positions:
        continue
    rec = sum(1 for s in TREATMENTS if var in sample_variants.get(s, set()))
    if rec >= 2:
        diff_somatic[var] = rec

print(f"\nTotal treatment somatic positions: {len(all_positions):,}")
print(f"Also in control (germline, removed): {len(all_positions & ctrl_positions):,}")
print(f"Treatment-specific (recurrence>=2): {len(diff_somatic):,}")

# Recurrence distribution
rec_dist = defaultdict(int)
for var, rec in diff_somatic.items():
    rec_dist[rec] += 1
print(f"\nRecurrence distribution:")
for k in sorted(rec_dist.keys(), reverse=True):
    print(f"  {k}/10: {rec_dist[k]:,}")

# Save differential somatic
import csv
outpath = f"{OUTDIR}/differential_somatic_candidates.csv"
with open(outpath, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['chrom', 'pos', 'ref', 'alt', 'recurrence'])
    for var in sorted(diff_somatic.keys(), key=lambda x: (-diff_somatic[x], x)):
        chrom, pos, ref, alt = var.split('_')
        writer.writerow([chrom, pos, ref, alt, diff_somatic[var]])

# Also classify by variant type
snp_count = sum(1 for v in diff_somatic if len(v.split('_')[2]) == 1 and len(v.split('_')[3]) == 1)
indel_count = len(diff_somatic) - snp_count
print(f"\nDifferential somatic SNPs: {snp_count:,}")
print(f"Differential somatic INDELs: {indel_count:,}")
print(f"\nSaved to: {outpath}")
PYEOF

echo ""
echo "============================================"
echo "Done: $(date)"
echo "Results in: $OUTDIR"
echo "============================================"
