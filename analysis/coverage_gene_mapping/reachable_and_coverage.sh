#!/usr/bin/env bash
# reachable_and_coverage.sh — self-contained reachable recall + gene coverage 재계산.
#
# recall_decomposition 프레임워크(04_compute_recall.py)와 동일 정의를 corrected candidate(C=1,356)에 적용.
#   recall = captured / (captured + B3_reachable)
#     captured     = candidate ∩ DNA truth (position)
#     B3_reachable = missed DNA-truth 중 mpileup으로 RNA 증거 충분(alt_AF≥0.30 또는 alt_DP≥10)
#   reachable    = captured ∪ B3_reachable
#   gene recall  = |genes(captured) ∩ genes(reachable)| / |genes(reachable)|
#
# 입력 (전부 committed 스크립트/소스로 재현되는 것):
#   - DNA truth      : recall_decomposition/data/dna_truth.bed     (ground_truth VCF → 01_build_inputs.sh; treat≥2 & ctrl=0)
#   - missed pileup  : recall_decomposition/data/bucket_b1b2b3.tsv  (treatment BAM mpileup → 02_classify_b1b2b3.sh)
#   - gene 주석      : dna_truth_coverage/genes.bed                 (GTF 유래)
#   - candidate      : differential ∩ imputed (= C 1,356) , -RNA editing (= D 961)  [reproduce_ablation.py와 동일]
# 키는 모두 탭 구분(chrom<TAB>pos) — chrom에 underscore 있는 alt-contig 안전.
set -euo pipefail
BASE=/home/yusanghyeon/RAT_project
BT=$BASE/miniforge3/envs/rnaseq/bin/bedtools
RD=$BASE/PHMG_IT/results/recall_decomposition/data
GENES=$BASE/PHMG_IT/results/dna_truth_coverage/genes.bed
DIFF=$BASE/PHMG_IT/results/mutect2/differential/differential_snps.csv
IMP=$BASE/PHMG_IT/results/imputed_differential/imputed_differential_snps.csv
OUT=$(dirname "$(readlink -f "$0")")
T=$(mktemp -d)

# --- candidate (chrom\tpos\tref\talt) ---
awk -F, 'NR>1{print $1"\t"$2"\t"$3"\t"$4}' "$DIFF" | sort -u > "$T/B.tsv"
awk -F, 'NR>1{print $1"\t"$2"\t"$3"\t"$4}' "$IMP"  | sort -u > "$T/I.tsv"
comm -12 "$T/B.tsv" "$T/I.tsv" > "$T/C.tsv"                                  # 1,356
awk -F'\t' '!(($3=="A"&&$4=="G")||($3=="T"&&$4=="C"))' "$T/C.tsv" > "$T/D.tsv" # 961

# --- DNA truth / B3 reachable positions (chrom\tpos) ---
awk -F'\t' '{print $1"\t"$3}' "$RD/dna_truth.bed" | sort -u > "$T/truth.pos"            # col3 = 1-based pos
awk -F'\t' 'NR>1 && $5>0 && ($6/$5)>=0.30 {print $1"\t"$2}' "$RD/bucket_b1b2b3.tsv" | sort -u > "$T/b3_af30.pos"
awk -F'\t' 'NR>1 && $6>=10               {print $1"\t"$2}' "$RD/bucket_b1b2b3.tsv" | sort -u > "$T/b3_dp10.pos"

pos2bed(){ awk -F'\t' '$2 ~ /^[0-9]+$/ {print $1"\t"$2-1"\t"$2}' "$1" | sort -k1,1 -k2,2n; }
genes_of(){ pos2bed "$1" | "$BT" intersect -a - -b "$GENES" -wb 2>/dev/null | awk '{print $NF}' | sort -u; }

SUM="$OUT/coverage_summary.tsv"
printf "candidate_set\tcaptured\treach_af30\trecall_af30\treach_dp10\trecall_dp10\tgenes_cap\tgenes_reach_af30\tgene_recall_af30\tgenes_reach_dp10\tgene_recall_dp10\n" > "$SUM"

run(){
  local name=$1 cand=$2
  cut -f1,2 "$cand" | sort -u > "$T/$name.cand.pos"
  comm -12 "$T/$name.cand.pos" "$T/truth.pos" > "$T/$name.cap.pos"
  local ncap; ncap=$(wc -l < "$T/$name.cap.pos")
  sort -u "$T/$name.cap.pos" "$T/b3_af30.pos" > "$T/$name.reach_af.pos"
  sort -u "$T/$name.cap.pos" "$T/b3_dp10.pos" > "$T/$name.reach_dp.pos"
  local raf rdp; raf=$(wc -l < "$T/$name.reach_af.pos"); rdp=$(wc -l < "$T/$name.reach_dp.pos")
  local prec_af prec_dp
  prec_af=$(awk "BEGIN{printf \"%.2f\",100*$ncap/$raf}")
  prec_dp=$(awk "BEGIN{printf \"%.2f\",100*$ncap/$rdp}")
  genes_of "$T/$name.cap.pos"      > "$T/$name.gcap"
  genes_of "$T/$name.reach_af.pos" > "$T/$name.graf"
  genes_of "$T/$name.reach_dp.pos" > "$T/$name.grdp"
  local gcap graf grdp; gcap=$(wc -l<"$T/$name.gcap"); graf=$(wc -l<"$T/$name.graf"); grdp=$(wc -l<"$T/$name.grdp")
  local gi_af gi_dp; gi_af=$(comm -12 "$T/$name.gcap" "$T/$name.graf"|wc -l); gi_dp=$(comm -12 "$T/$name.gcap" "$T/$name.grdp"|wc -l)
  local grec_af grec_dp
  grec_af=$(awk "BEGIN{printf \"%.2f\",100*$gi_af/$graf}")
  grec_dp=$(awk "BEGIN{printf \"%.2f\",100*$gi_dp/$grdp}")
  echo "=== $name : captured=$ncap ==="
  echo "  position recall: AF≥0.30 ${ncap}/${raf}=${prec_af}%   alt_DP≥10 ${ncap}/${rdp}=${prec_dp}%"
  echo "  gene recall    : AF≥0.30 ${gi_af}/${graf}=${grec_af}%   alt_DP≥10 ${gi_dp}/${grdp}=${grec_dp}%   (captured genes=${gcap})"
  echo
  printf "%s\t%s\t%s\t%s%%\t%s\t%s%%\t%s\t%s\t%s%%\t%s\t%s%%\n" \
    "$name" "$ncap" "$raf" "$prec_af" "$rdp" "$prec_dp" "$gcap" "$graf" "$grec_af" "$grdp" "$grec_dp" >> "$SUM"
}

run "C_1356_main" "$T/C.tsv"
run "D_961_alt"   "$T/D.tsv"
echo "요약: $SUM"
rm -rf "$T"
