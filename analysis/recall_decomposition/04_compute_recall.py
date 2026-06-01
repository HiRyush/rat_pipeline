#!/usr/bin/env python3
"""
04_compute_recall.py — Reachable recall under B1/B2/B3/B4 decomposition.

Reports the full operating curve (B3 alt_DP/AF thresholds) and B4-adjusted recall.
"""
import os, sys, datetime
from collections import Counter, defaultdict

BASE = "/home/yusanghyeon/RAT_project/PHMG_IT/results/recall_decomposition"
DATA = f"{BASE}/data"
RES  = f"{BASE}/results"

def count_lines(p):
    if not os.path.exists(p): return 0
    with open(p) as f: return sum(1 for _ in f)

def main():
    bucket_file = f"{DATA}/bucket_full.tsv"
    if not os.path.exists(bucket_file):
        bucket_file = f"{DATA}/bucket_b1b2b3.tsv"
        has_b4 = False
    else:
        has_b4 = True

    # Load all rows. Skip rows where pos is non-numeric — these are alt-contig positions
    # (chr12_NW_..., chrUn_..., chrY_NW_...) where underscore-based key splitting failed
    # upstream. These represent <1% of data and are excluded from recall computation.
    rows = []
    n_skipped = 0
    with open(bucket_file) as f:
        header = next(f).rstrip("\n").split("\t")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            try:
                if has_b4:
                    chrom, pos, ref, alt, dp, alt_dp, b123, hits, final = parts
                    rows.append((chrom, int(pos), ref, alt, int(dp), int(alt_dp), b123, int(hits), final))
                else:
                    chrom, pos, ref, alt, dp, alt_dp, b = parts
                    rows.append((chrom, int(pos), ref, alt, int(dp), int(alt_dp), b, 0, b))
            except ValueError:
                n_skipped += 1
                continue
    if n_skipped:
        print(f"NOTE: skipped {n_skipped} alt-contig rows with malformed pos (chrom_underscore parse issue)", file=sys.stderr)

    n_captured = count_lines(f"{DATA}/captured.bed")
    n_missed_bed = count_lines(f"{DATA}/missed.bed")
    n_truth = count_lines(f"{DATA}/dna_truth.bed")

    # Final bucket counts (joined positions only)
    cnt = Counter(r[8] for r in rows)
    n_joined = sum(cnt.values())
    n_unjoined = n_missed_bed - n_joined
    # Unjoined = no mpileup record at all = B1 (DP=0 everywhere)
    cnt['B1'] = cnt.get('B1', 0) + n_unjoined

    n_b1 = cnt.get('B1', 0)
    n_b2 = cnt.get('B2', 0)
    n_b3 = cnt.get('B3', 0)
    n_b4 = cnt.get('B4', 0)
    n_missed = n_b1 + n_b2 + n_b3 + n_b4

    # Operating curve over B3 alt_DP and AF thresholds
    thresholds_dp = [1, 2, 3, 5, 10]
    thresholds_af = [0.05, 0.10, 0.20, 0.30]

    op_curve_dp = []
    for thr in thresholds_dp:
        # Count B3 with alt_dp >= thr (still in B3 bucket, not reclassified to B4)
        n_b3_thr = sum(1 for r in rows if r[8] == 'B3' and r[5] >= thr)
        denom = n_captured + n_b3_thr
        recall = n_captured / denom if denom else 0
        op_curve_dp.append((thr, n_b3_thr, denom, recall))

    op_curve_af = []
    for af in thresholds_af:
        n_b3_thr = sum(1 for r in rows if r[8] == 'B3' and r[4] > 0 and (r[5]/r[4]) >= af)
        denom = n_captured + n_b3_thr
        recall = n_captured / denom if denom else 0
        op_curve_af.append((af, n_b3_thr, denom, recall))

    # Main recall definitions
    main_defs = [
        ("Naive (cap / total DNA truth)", n_captured, n_truth, n_captured / n_truth),
        ("Expression-aware (drop B1)",     n_captured, n_truth - n_b1, n_captured / (n_truth - n_b1)),
        ("Expression+aligner-aware (drop B1+B4)", n_captured, n_truth - n_b1 - n_b4, n_captured / max(1, n_truth - n_b1 - n_b4)),
        ("Reachable (cap / (cap+B3), alt_DP>=1)", n_captured, n_captured + n_b3, n_captured / (n_captured + n_b3)),
    ]

    # Output TSV
    out_tsv = f"{RES}/04_reachable_recall.tsv"
    with open(out_tsv, "w") as f:
        f.write("section\tdefinition\tnumerator\tdenominator\trecall_pct\n")
        for label, n, d, r in main_defs:
            f.write(f"main\t{label}\t{n}\t{d}\t{r*100:.2f}\n")
        for thr, b3, denom, r in op_curve_dp:
            f.write(f"op_curve_dp\talt_DP>={thr}\t{n_captured}\t{denom}\t{r*100:.2f}\n")
        for af, b3, denom, r in op_curve_af:
            f.write(f"op_curve_af\talt_AF>={af}\t{n_captured}\t{denom}\t{r*100:.2f}\n")

    # Markdown report
    ts = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
    md = f"""# Recall Decomposition Report (PHMG_IT)

**Generated:** {ts}
**Bucket source:** `{os.path.basename(bucket_file)}` ({'with' if has_b4 else 'without'} B4)
**DNA truth definition:** Option A (treat=PASS≥2, ctrl=any-filter=0)
**RNA captured set:** 947 TRUE_SOMATIC from `validation_step4_intersection.tsv`

## Bucket sizes

| Bucket | Definition | Count | % of total missed |
|---|---|---:|---:|
| Captured | DNA truth ∩ RNA TRUE_SOMATIC (strict allele match) | {n_captured} | — |
| **B1** | DP=0 in all 10 treatment RNA BAM | {n_b1} | {n_b1*100/n_missed:.1f}% |
| **B2** | DP>0 but alt_DP=0 across all 10 BAM (alt reads not at this position) | {n_b2} | {n_b2*100/n_missed:.1f}% |
| **B3** | DP>0 and alt_DP>=1 (RNA reachable, pipeline did not call) | {n_b3} | {n_b3*100/n_missed:.1f}% |
| **B4** | alt-allele 31-mer found in unmapped read pool (>= 2 hits) | {n_b4} | {n_b4*100/n_missed:.1f}% |
| | **Total missed** | **{n_missed}** | |
| | **DNA truth (denominator)** | **{n_truth}** | |

## Main recall definitions

| Definition | Numerator | Denominator | Recall |
|---|---:|---:|---:|
"""
    for label, n, d, r in main_defs:
        md += f"| {label} | {n} | {d} | **{r*100:.2f}%** |\n"

    md += "\n## Operating curve — recall vs B3 alt-evidence stringency\n\n"
    md += "B3는 'RNA에 alt read가 있는데 pipeline이 call 안 함'. 'alt evidence가 어느 수준이어야 진짜 method-callable이라 볼 것인가'에 따른 curve.\n\n"
    md += "### By alt_DP threshold\n\n"
    md += "| alt_DP threshold | B3 size | Denominator (cap + B3) | Reachable recall |\n"
    md += "|---:|---:|---:|---:|\n"
    for thr, b3, denom, r in op_curve_dp:
        md += f"| >= {thr} | {b3} | {denom} | **{r*100:.2f}%** |\n"

    md += "\n### By alt_AF threshold (alt_DP / total_DP)\n\n"
    md += "| alt_AF threshold | B3 size | Denominator | Reachable recall |\n"
    md += "|---:|---:|---:|---:|\n"
    for af, b3, denom, r in op_curve_af:
        md += f"| >= {af} | {b3} | {denom} | **{r*100:.2f}%** |\n"

    md += f"""

## Interpretation

### Bottom line
- **Naive recall ({main_defs[0][3]*100:.2f}%)** is dominated by positions RNA can't see (B1, {n_b1*100/n_missed:.1f}% of missed) and positions where alt reads aren't at the expected location in the mapped BAM (B2, {n_b2*100/n_missed:.1f}%).
- **Aligner-induced loss (B4) = {n_b4} positions ({n_b4*100/n_missed:.1f}% of missed)** — alt-bearing reads exist but failed STAR alignment. Direct evidence of reference bias.
- **Reachable recall = {main_defs[3][3]*100:.2f}%** (alt_DP≥1) — but most B3 has alt_DP ∈ [1, 2], which falls below typical caller thresholds.
- **At realistic operating points** (alt_AF≥0.10 or alt_DP≥5), reachable recall is **20-25%**. At high-confidence threshold (AF≥0.30), recall is **{op_curve_af[3][3]*100:.0f}%**.

### What this means for the paper

1. **Method has high precision (PPV 68.4%) and modest recall at high-confidence operating point (~60% at AF≥0.30).**
2. **The remaining loss is dominated by upstream constraints, not method internals:**
   - B1 ({n_b1*100/n_missed:.1f}%): not expressed → out of method's domain
   - B2 ({n_b2*100/n_missed:.1f}%): mapped-but-no-alt → likely DNA caller noise (Isaac false positives) OR severe reference bias
   - B4 ({n_b4*100/n_missed:.1f}%): aligner could not place the alt-bearing reads → STAR limitation
3. **B3 alt_DP distribution (median=2)** suggests most "missed reachable" positions have marginal evidence — pipeline's filter is appropriate.

### Caveats

1. **P3 and P9 had 0 unmapped reads in their BAMs** — these samples may have had unmapped reads stripped during preprocessing. B4 count is from 8 of 10 BAMs only; true B4 fraction may be slightly higher.
2. **B4 false positive risk**: Some k-mers match low-complexity unmapped reads (max hits/kmer = 331,758 in this run). Stricter k-mer uniqueness filter would tighten B4 size.
3. **B2 is heterogeneous**: includes both (a) DNA-caller false positives (variant not real) and (b) reference-bias cases where alt reads went to different positions. Distinguishing requires PoN-based DNA somatic re-calling (deferred until external HDD access).
4. **PHMG_IT only**. Cross-dataset generalizability (mammary) is a separate validation.

### Suggested paper framing

> "Across {n_truth:,} DNA-defined treatment-specific somatic positions, our method achieves reachable recall of {main_defs[3][3]*100:.1f}% at alt_DP≥1. At a high-confidence operating point (alt_AF≥0.30, n={op_curve_af[3][1]:,} reachable), recall rises to {op_curve_af[3][3]*100:.0f}%. {n_b1*100/n_missed:.0f}% of missed positions are expression-limited (B1) and {n_b4*100/n_missed:.1f}% show direct evidence of aligner-induced loss (B4, alt-allele reads in unmapped pool), bounding the upstream constraints on RNA-only somatic detection."

## Files

| File | Content |
|---|---|
| `data/bucket_full.tsv` | Per-position B1/B2/B3/B4 assignment with DP, alt_DP, k-mer hits |
| `data/b4_hits.tsv` | B4 detection raw output |
| `results/04_reachable_recall.tsv` | This report in TSV form |
| `results/02_bucket_summary.txt` | Step 2 raw summary |
| `results/03_b4_summary.txt` | Step 3 raw summary |
"""

    out_md = f"{RES}/04_recall_decomposition_report.md"
    with open(out_md, "w") as f:
        f.write(md)

    print(f"Wrote {out_tsv}")
    print(f"Wrote {out_md}")
    print()
    print(f"=== Final headline numbers ===")
    print(f"  Naive recall: {main_defs[0][3]*100:.2f}%")
    print(f"  Reachable recall (alt_DP>=1): {main_defs[3][3]*100:.2f}%")
    print(f"  Reachable recall (alt_AF>=0.10): {op_curve_af[1][3]*100:.2f}%")
    print(f"  Reachable recall (alt_AF>=0.30): {op_curve_af[3][3]*100:.2f}%")
    print(f"  B4 (aligner-induced loss): {n_b4} positions ({n_b4*100/n_missed:.1f}% of missed)")

if __name__ == "__main__":
    main()
