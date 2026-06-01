#!/usr/bin/env python3
"""
05_b4_salvage_to_numerator.py — Filter B4 for quality and promote robust B4 to numerator.

Quality filters applied to B4 positions:
  - Hit count in range [MIN_HITS, MAX_HITS]: exclude low-complexity k-mers that hit huge.
  - K-mer Shannon entropy >= MIN_ENTROPY: exclude homopolymer / repetitive k-mers.
  - K-mer (fwd) not appearing in reference at >1 location: skip non-unique k-mers (deferred,
    requires BWA index lookup; for now relies on entropy + complexity).

After filtering, robust B4 positions get added to numerator and recall curve is recomputed.

Outputs:
  data/b4_robust.tsv — filtered B4 with kept/dropped + reason
  data/numerator_v2.tsv — expanded numerator (captured + robust B4)
  results/05_salvage_summary.txt
  results/05_recall_with_salvage.tsv
"""
import os, sys, math, datetime
from collections import Counter, defaultdict

BASE = "/home/yusanghyeon/RAT_project/PHMG_IT/results/recall_decomposition"
DATA = f"{BASE}/data"
RES  = f"{BASE}/results"

KMER_LEN = 31
FLANK = 15

# Filter thresholds
MIN_HITS = 2
MAX_HITS = 200             # exclude crazy-high (likely low-complexity k-mers)
MIN_ENTROPY = 1.5          # Shannon entropy of nucleotides in k-mer; full random is ~2.0
MAX_HOMOPOLYMER = 8        # disallow runs of same base >= 8 in k-mer

def revcomp(s):
    return s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]

def shannon_entropy(seq):
    if not seq: return 0
    c = Counter(seq)
    L = len(seq)
    return -sum((n/L) * math.log2(n/L) for n in c.values() if n)

def max_homopolymer(seq):
    best = cur = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i-1]:
            cur += 1
            if cur > best: best = cur
        else:
            cur = 1
    return best

# 1. Rebuild k-mer per position from flanks file (need k-mer sequence for entropy check)
print("[1/4] Rebuilding alt k-mers from flanks...", file=sys.stderr)
position_kmer = {}  # (chrom, pos, ref, alt) -> kmer_fwd
with open(f"{DATA}/b1b2_flanks.tsv") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        name = parts[0].split("::")[0]
        seq = parts[1].upper()
        if len(seq) != KMER_LEN: continue
        toks = name.split("_")
        if len(toks) < 5: continue
        chrom = "_".join(toks[:-4])
        pos = toks[-4]; refb = toks[-3]; altb = toks[-2]; bucket = toks[-1]
        if seq[FLANK] != refb: continue
        if altb not in "ACGT" or len(altb) != 1: continue
        kmer = seq[:FLANK] + altb + seq[FLANK+1:]
        if "N" in kmer: continue
        position_kmer[(chrom, pos, refb, altb)] = kmer
print(f"  Loaded {len(position_kmer)} k-mers", file=sys.stderr)

# 2. Apply quality filters to B4 positions in bucket_full.tsv
print("[2/4] Filtering B4 quality...", file=sys.stderr)
b4_raw = []   # (chrom, pos, ref, alt, hits)
with open(f"{DATA}/bucket_full.tsv") as f:
    next(f)
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if parts[8] != "B4": continue
        chrom, pos, ref, alt, dp, alt_dp, b123, hits, final = parts
        b4_raw.append((chrom, pos, ref, alt, int(hits)))
print(f"  Raw B4: {len(b4_raw)}", file=sys.stderr)

robust_b4 = []
drop_reason_counts = Counter()
with open(f"{DATA}/b4_robust.tsv", "w") as out:
    out.write("chrom\tpos\tref\talt\thits\tkmer\tentropy\tmax_homopolymer\tstatus\treason\n")
    for chrom, pos, ref, alt, hits in b4_raw:
        kmer = position_kmer.get((chrom, pos, ref, alt))
        if kmer is None:
            ent = 0.0; hom = 0
            status = "DROPPED"; reason = "no_kmer_found"
        else:
            ent = shannon_entropy(kmer)
            hom = max_homopolymer(kmer)

            if hits < MIN_HITS:
                status, reason = "DROPPED", f"hits<{MIN_HITS}"
            elif hits > MAX_HITS:
                status, reason = "DROPPED", f"hits>{MAX_HITS}_low_complexity"
            elif ent < MIN_ENTROPY:
                status, reason = "DROPPED", f"entropy<{MIN_ENTROPY:.1f}"
            elif hom > MAX_HOMOPOLYMER:
                status, reason = "DROPPED", f"homopolymer>{MAX_HOMOPOLYMER}"
            else:
                status, reason = "KEPT", "ok"
                robust_b4.append((chrom, pos, ref, alt, hits))

        drop_reason_counts[reason] += 1
        out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{hits}\t{kmer or 'NA'}\t{ent:.2f}\t{hom}\t{status}\t{reason}\n")

print(f"  Robust B4 (KEPT): {len(robust_b4)}", file=sys.stderr)
print("  Filter breakdown:", file=sys.stderr)
for reason, n in drop_reason_counts.most_common():
    print(f"    {reason}: {n}", file=sys.stderr)

# 3. Build expanded numerator (captured + robust B4)
print("[3/4] Building expanded numerator...", file=sys.stderr)
captured = []
with open(f"{DATA}/captured.bed") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 5:
            captured.append((parts[0], parts[2], parts[3], parts[4]))  # chrom, pos, ref, alt
        elif len(parts) >= 4:
            # alternate format
            captured.append(tuple(parts[:4]))

with open(f"{DATA}/numerator_v2.tsv", "w") as f:
    f.write("chrom\tpos\tref\talt\tsource\n")
    for c, p, r, a in captured:
        f.write(f"{c}\t{p}\t{r}\t{a}\tcaptured\n")
    for c, p, r, a, h in robust_b4:
        f.write(f"{c}\t{p}\t{r}\t{a}\tb4_salvage\n")

n_cap = len(captured)
n_b4r = len(robust_b4)
n_num_v2 = n_cap + n_b4r
print(f"  Captured: {n_cap}", file=sys.stderr)
print(f"  Robust B4 salvage: {n_b4r}", file=sys.stderr)
print(f"  Numerator v2: {n_num_v2}", file=sys.stderr)

# 4. Recompute recall curve with new numerator
print("[4/4] Recomputing recall curve with expanded numerator...", file=sys.stderr)

# Load bucket_full to get B3 by AF threshold
rows = []
with open(f"{DATA}/bucket_full.tsv") as f:
    next(f)
    for line in f:
        parts = line.rstrip("\n").split("\t")
        try:
            chrom, pos, ref, alt = parts[0], int(parts[1]), parts[2], parts[3]
            dp, alt_dp = int(parts[4]), int(parts[5])
            final = parts[8]
            rows.append((chrom, pos, ref, alt, dp, alt_dp, final))
        except ValueError:
            continue

# Numerator universes
def recall_at(n_num, n_b3_thr):
    denom = n_num + n_b3_thr
    return n_num, denom, (n_num / denom * 100 if denom else 0)

dp_thresholds = [1, 2, 3, 5, 10]
af_thresholds = [0.05, 0.10, 0.20, 0.30]

with open(f"{RES}/05_recall_with_salvage.tsv", "w") as f:
    f.write("section\tdefinition\tnumerator\tdenominator\trecall_pct\n")
    # baseline (captured only)
    for thr in dp_thresholds:
        n_b3 = sum(1 for r in rows if r[6] == 'B3' and r[5] >= thr)
        n, d, r = recall_at(n_cap, n_b3)
        f.write(f"baseline_dp\talt_DP>={thr}\t{n}\t{d}\t{r:.2f}\n")
    for af in af_thresholds:
        n_b3 = sum(1 for r in rows if r[6] == 'B3' and r[4] > 0 and (r[5]/r[4]) >= af)
        n, d, r = recall_at(n_cap, n_b3)
        f.write(f"baseline_af\talt_AF>={af}\t{n}\t{d}\t{r:.2f}\n")
    # salvage
    for thr in dp_thresholds:
        n_b3 = sum(1 for r in rows if r[6] == 'B3' and r[5] >= thr)
        n, d, r = recall_at(n_num_v2, n_b3)
        f.write(f"salvage_dp\talt_DP>={thr}\t{n}\t{d}\t{r:.2f}\n")
    for af in af_thresholds:
        n_b3 = sum(1 for r in rows if r[6] == 'B3' and r[4] > 0 and (r[5]/r[4]) >= af)
        n, d, r = recall_at(n_num_v2, n_b3)
        f.write(f"salvage_af\talt_AF>={af}\t{n}\t{d}\t{r:.2f}\n")

# Summary
ts = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
with open(f"{RES}/05_salvage_summary.txt", "w") as f:
    f.write(f"B4 salvage to numerator — summary\n")
    f.write(f"Generated: {ts}\n\n")
    f.write(f"Filter parameters:\n")
    f.write(f"  Hit count range: [{MIN_HITS}, {MAX_HITS}]\n")
    f.write(f"  Min Shannon entropy: {MIN_ENTROPY}\n")
    f.write(f"  Max homopolymer run: {MAX_HOMOPOLYMER}\n\n")
    f.write(f"Raw B4: {len(b4_raw)}\n")
    f.write(f"Robust B4 (KEPT): {n_b4r}\n")
    f.write(f"Drop reasons:\n")
    for reason, n in drop_reason_counts.most_common():
        f.write(f"  {reason}: {n}\n")
    f.write(f"\nNumerator: {n_cap} (baseline) → {n_num_v2} (with B4 salvage) (+{n_b4r})\n\n")
    f.write(f"Recall comparison at AF>=0.30:\n")
    n_b3_af30 = sum(1 for r in rows if r[6] == 'B3' and r[4] > 0 and (r[5]/r[4]) >= 0.30)
    baseline = n_cap / (n_cap + n_b3_af30) * 100
    salvage  = n_num_v2 / (n_num_v2 + n_b3_af30) * 100
    f.write(f"  Baseline:  {n_cap}/{n_cap + n_b3_af30} = {baseline:.2f}%\n")
    f.write(f"  Salvage:   {n_num_v2}/{n_num_v2 + n_b3_af30} = {salvage:.2f}%  (delta +{salvage-baseline:.2f} pp)\n")

print()
print("=== Headline ===")
print(f"  Raw B4: {len(b4_raw)}")
print(f"  Robust B4 (kept): {n_b4r}")
print(f"  Numerator: {n_cap} → {n_num_v2}")
n_b3_af30 = sum(1 for r in rows if r[6] == 'B3' and r[4] > 0 and (r[5]/r[4]) >= 0.30)
print(f"  Recall@AF>=0.30: {n_cap/(n_cap+n_b3_af30)*100:.2f}% → {n_num_v2/(n_num_v2+n_b3_af30)*100:.2f}%")
n_b3_af10 = sum(1 for r in rows if r[6] == 'B3' and r[4] > 0 and (r[5]/r[4]) >= 0.10)
print(f"  Recall@AF>=0.10: {n_cap/(n_cap+n_b3_af10)*100:.2f}% → {n_num_v2/(n_num_v2+n_b3_af10)*100:.2f}%")
n_b3_af05 = sum(1 for r in rows if r[6] == 'B3' and r[4] > 0 and (r[5]/r[4]) >= 0.05)
print(f"  Recall@AF>=0.05: {n_cap/(n_cap+n_b3_af05)*100:.2f}% → {n_num_v2/(n_num_v2+n_b3_af05)*100:.2f}%")
