#!/usr/bin/env bash
# recompute_coverage.sh — 재현 numerator(961 D / 1,356 C)로 coverage·gene mapping 재계산.
#
# 폐기된 1,392 기반 headline(position 70.59% / gene 74.44%)을 대체한다.
# Metric 정의(통일): recall = |후보 ∩ reachable| / |reachable|  (gene-level은 gene 집합 기준).
#   gene-level @ AF≥0.30 = |genes(후보) ∩ genes(reachable)| / |genes(reachable)|  ← 원래 198/266=74.44% 정의와 동일.
#
# 입력 = 현재 디스크 산출물:
#   - 후보: differential ∩ imputed (= C, 1,356) , RNA editing 제거 (= D, 961)  ← reproduce_ablation.py와 동일 정의
#   - reachable DNA truth: gene_mapping/reachable_{af30,dp10}.bed  (사전 산출 artifact)
#   - gene 주석: gene_mapping/genes.sorted.bed
#   - gene 분모: gene_mapping/genes_reachable_{af30,dp10}.txt, genes_dna_truth.txt
set -euo pipefail
BASE=/home/yusanghyeon/RAT_project
BT=$BASE/miniforge3/envs/rnaseq/bin/bedtools
GM=$BASE/pipeline_candidates/results/observation_first/gene_mapping
DIFF=$BASE/PHMG_IT/results/mutect2/differential/differential_snps.csv
IMP=$BASE/PHMG_IT/results/imputed_differential/imputed_differential_snps.csv
OUT=$(dirname "$(readlink -f "$0")")
T=$(mktemp -d)

# --- 후보 재구축 (C, D) ---
awk -F, 'NR>1{print $1"_"$2"_"$3">"$4}' "$DIFF" | sort -u > "$T/B.keys"
awk -F, 'NR>1{print $1"_"$2"_"$3">"$4}' "$IMP"  | sort -u > "$T/I.keys"
comm -12 "$T/B.keys" "$T/I.keys" > "$T/C.keys"                                   # 1,356
awk -F'[_>]' '!(($3=="A"&&$4=="G")||($3=="T"&&$4=="C"))' "$T/C.keys" > "$T/D.keys" # 961

key2bed(){ awk -F'[_>]' '{print $1"\t"$2-1"\t"$2"\t"$3">"$4}' "$1" | sort -k1,1 -k2,2n; }
genes_of(){ "$BT" intersect -a "$1" -b "$GM/genes.sorted.bed" -wb 2>/dev/null | awk '{print $NF}' | sort -u; }
poskey(){ awk '{print $1"_"$3}' "$1" | sort -u; }   # bed -> chr_pos(1-based)

REACH_AF30_POS=$(poskey "$GM/reachable_af30.bed" | wc -l)
REACH_DP10_POS=$(poskey "$GM/reachable_dp10.bed" | wc -l)
G_REACH_AF30=$(sort -u "$GM/genes_reachable_af30.txt" | wc -l)
G_REACH_DP10=$(sort -u "$GM/genes_reachable_dp10.txt" | wc -l)
G_DNA=$(sort -u "$GM/genes_dna_truth.txt" | wc -l)
DNA_POS=$(poskey "$GM/dna_truth.sorted.bed" | wc -l)

printf "reachable: AF30_pos=%s DP10_pos=%s | genes: reach_af30=%s reach_dp10=%s dna_all=%s | dna_pos=%s\n\n" \
  "$REACH_AF30_POS" "$REACH_DP10_POS" "$G_REACH_AF30" "$G_REACH_DP10" "$G_DNA" "$DNA_POS"

SUM="$OUT/coverage_summary.tsv"
printf "candidate_set\tn_cand\tn_genes\tpos_recall_af30\tpos_recall_dp10\tgene_recall_af30\tgene_recall_dp10\tgene_recall_dnaall\tcaptured_strict\n" > "$SUM"

run_set(){
  local name=$1 keys=$2
  key2bed "$keys" > "$T/$name.bed"
  local n; n=$(wc -l < "$keys")
  poskey "$T/$name.bed" > "$T/$name.pos"
  local cap_af30 cap_dp10 cap_dna
  cap_af30=$(comm -12 "$T/$name.pos" <(poskey "$GM/reachable_af30.bed") | wc -l)
  cap_dp10=$(comm -12 "$T/$name.pos" <(poskey "$GM/reachable_dp10.bed") | wc -l)
  cap_dna=$(comm -12 "$T/$name.pos" <(poskey "$GM/dna_truth.sorted.bed") | wc -l)
  genes_of "$T/$name.bed" > "$T/$name.genes"
  local ng; ng=$(wc -l < "$T/$name.genes")
  local gaf gdp gall
  gaf=$(comm -12 "$T/$name.genes" <(sort -u "$GM/genes_reachable_af30.txt") | wc -l)
  gdp=$(comm -12 "$T/$name.genes" <(sort -u "$GM/genes_reachable_dp10.txt") | wc -l)
  gall=$(comm -12 "$T/$name.genes" <(sort -u "$GM/genes_dna_truth.txt") | wc -l)
  local pr_af pr_dp gr_af gr_dp gr_all
  pr_af=$(awk "BEGIN{printf \"%.2f\",100*$cap_af30/$REACH_AF30_POS}")
  pr_dp=$(awk "BEGIN{printf \"%.2f\",100*$cap_dp10/$REACH_DP10_POS}")
  gr_af=$(awk "BEGIN{printf \"%.2f\",100*$gaf/$G_REACH_AF30}")
  gr_dp=$(awk "BEGIN{printf \"%.2f\",100*$gdp/$G_REACH_DP10}")
  gr_all=$(awk "BEGIN{printf \"%.2f\",100*$gall/$G_DNA}")
  echo "=== $name : 후보 $n / 매핑 gene $ng ==="
  echo "  position recall: AF≥0.30 ${cap_af30}/${REACH_AF30_POS}=${pr_af}%   DP≥10 ${cap_dp10}/${REACH_DP10_POS}=${pr_dp}%"
  echo "  gene recall    : AF≥0.30 ${gaf}/${G_REACH_AF30}=${gr_af}%   DP≥10 ${gdp}/${G_REACH_DP10}=${gr_dp}%   DNA-all ${gall}/${G_DNA}=${gr_all}%"
  echo "  strict captured (∩ DNA truth, position): ${cap_dna}"
  echo
  printf "%s\t%s\t%s\t%s%%\t%s%%\t%s%%\t%s%%\t%s%%\t%s\n" "$name" "$n" "$ng" "$pr_af" "$pr_dp" "$gr_af" "$gr_dp" "$gr_all" "$cap_dna" >> "$SUM"
}

run_set "C_1356_pre_editing" "$T/C.keys"
run_set "D_961_final"        "$T/D.keys"
echo "요약 저장: $SUM"
rm -rf "$T"
