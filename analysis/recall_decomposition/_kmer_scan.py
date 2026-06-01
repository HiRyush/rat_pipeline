import sys, os, subprocess
from collections import defaultdict

FLANK = int(os.environ['FLANK'])
KMER = 2*FLANK + 1
ROOT = os.environ['ROOT']
DATA = os.environ['DATA']
LOGS = os.environ['LOGS']
TREAT = os.environ['TREAT'].split()

# 1. Build alt 31-mer set + mapping kmer → position_key
def revcomp(s):
    return s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]

kmer_to_keys = defaultdict(set)  # kmer -> set of position keys
position_info = {}                # key -> (chrom, pos, ref, alt, orig_bucket)

skipped = 0; built = 0
with open(f"{DATA}/b1b2_flanks.tsv") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        # name col format: chrom_pos_ref_alt_bucket  followed by ::chrom:start-end (bedtools -name appends)
        name = parts[0].split("::")[0]
        seq = parts[1].upper()
        if len(seq) != KMER:
            skipped += 1; continue
        toks = name.split("_")
        if len(toks) < 5:
            skipped += 1; continue
        chrom = "_".join(toks[:-4])
        pos = toks[-4]; refb = toks[-3]; altb = toks[-2]; bucket = toks[-1]
        if seq[FLANK] != refb:
            skipped += 1; continue
        if altb not in "ACGT" or len(altb) != 1:
            skipped += 1; continue
        kmer_fwd = seq[:FLANK] + altb + seq[FLANK+1:]
        if "N" in kmer_fwd:
            skipped += 1; continue
        kmer_rc = revcomp(kmer_fwd)
        key = f"{chrom}:{pos}:{refb}:{altb}"
        position_info[key] = (chrom, pos, refb, altb, bucket)
        kmer_to_keys[kmer_fwd].add(key)
        kmer_to_keys[kmer_rc].add(key)
        built += 1
print(f"  Built {built} positions, {len(kmer_to_keys)} unique kmers, skipped {skipped}", file=sys.stderr)

# 2. Stream unmapped reads from each treat BAM, scan for k-mer hits
hits = defaultdict(int)  # position_key -> total read-hit count
total_reads = 0

samtools = "samtools"
for sample in TREAT:
    bam = f"{ROOT}/results/mutect2/markdup/{sample}.dedup.bam"
    print(f"  Scanning unmapped reads from {sample}...", file=sys.stderr)
    proc = subprocess.Popen(
        [samtools, "view", "-f", "4", "-F", "256", bam],
        stdout=subprocess.PIPE, text=True
    )
    n_reads_sample = 0
    for line in proc.stdout:
        cols = line.split("\t", 11)
        if len(cols) < 10: continue
        seq = cols[9]
        if len(seq) < KMER: continue
        n_reads_sample += 1
        # Slide window
        for i in range(len(seq) - KMER + 1):
            kmer = seq[i:i+KMER]
            if kmer in kmer_to_keys:
                for key in kmer_to_keys[kmer]:
                    hits[key] += 1
    proc.wait()
    total_reads += n_reads_sample
    print(f"    {sample}: {n_reads_sample} reads scanned, cumulative hits to date: {sum(hits.values())}", file=sys.stderr)

print(f"  Total unmapped reads scanned: {total_reads}", file=sys.stderr)
print(f"  Total positions with hits: {len(hits)}", file=sys.stderr)

# 3. Emit per-position hit table + final bucket reclassification
MIN_HITS = int(os.environ.get('MIN_HITS', '2'))
out_hits = open(f"{DATA}/b4_hits.tsv", "w")
out_hits.write("chrom\tpos\tref\talt\torig_bucket\tunmapped_kmer_hits\tcalled_b4\n")
for key, info in position_info.items():
    chrom, pos, refb, altb, bucket = info
    n = hits.get(key, 0)
    is_b4 = (n >= MIN_HITS)
    out_hits.write(f"{chrom}\t{pos}\t{refb}\t{altb}\t{bucket}\t{n}\t{int(is_b4)}\n")
out_hits.close()
print("  Wrote b4_hits.tsv", file=sys.stderr)
