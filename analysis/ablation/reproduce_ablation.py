#!/usr/bin/env python3
"""
reproduce_ablation.py — Observation-first + Imputation-as-Filter ablation 캐스케이드 재현.

발표용 headline 숫자를 "현재 디스크의 중간산출물"에서 결정론적으로 재현한다.
2026-05-28 인라인 ablation(휘발됨)을 대체하는 committed·reproducible 버전.

캐스케이드:
  A. MuTect2 union sites                         (raw observation)
  B. + Fisher differential (treatment vs control)
  C. + Imputation-as-filter (∩ imputed differential)   ← 핵심
  D. + RNA editing filter (A>G / T>C 제거)              ← 최종 candidate set

각 단계에서 DNA ground truth(15 samples)로 TP/PPV 산출:
  TRUE_SOMATIC  = treatment DNA에 confirmed & control DNA엔 없음
  GERMLINE_LEAK = treatment & control 둘 다 confirmed
  RNA_ONLY      = 어느 DNA에도 없음
PPV = TRUE_SOMATIC / total.

키 형식: chrom_pos_REF>ALT  (SNV 기준; 4-tuple strict).
실행: conda activate rnaseq 후 python reproduce_ablation.py
"""
import os, subprocess, csv, sys

BASE = "/home/yusanghyeon/RAT_project"
BCF  = f"{BASE}/miniforge3/envs/rnaseq/bin/bcftools"
DIFF = f"{BASE}/PHMG_IT/results/mutect2/differential/differential_snps.csv"
IMP  = f"{BASE}/PHMG_IT/results/imputed_differential/imputed_differential_snps.csv"
UNION= f"{BASE}/PHMG_IT/results/mutect2/union_sites.vcf.gz"
GT   = f"{BASE}/PHMG_IT/results/ground_truth"
OUT  = os.path.dirname(os.path.abspath(__file__))

CONTROLS   = ["C1","C2","C3","C4","C5"]
TREATMENTS = ["P1","P2","P3","P4","P5","P6","P7","P8","P9","P10"]

def key(c,p,r,a): return f"{c}_{p}_{r}>{a}"

def load_csv_keys(path):
    ks=set(); meta={}
    with open(path) as f:
        for row in csv.DictReader(f):
            k=key(row['chrom'],row['pos'],row['ref'],row['alt'])
            ks.add(k); meta[k]=row
    return ks, meta

def sample_dna_keys(sample):
    """샘플 DNA VCF에서 alt allele 보유(het/hom) 변이 키 집합 (1회 호출)."""
    vcf=f"{GT}/{sample}_dna_snps.vcf.gz"
    if not os.path.exists(vcf): return set()
    r=subprocess.run([BCF,"query","-i",'GT[*]="alt"',
                      "-f","%CHROM\\_%POS\\_%REF>%ALT\\n",vcf],
                     capture_output=True,text=True)
    return {l for l in r.stdout.split("\n") if l}

print("="*70); print("Ablation 캐스케이드 재현 (현재 디스크 입력)"); print("="*70)

# A
nA = int(subprocess.run([BCF,"view","-H",UNION],capture_output=True,text=True).stdout.count("\n")) \
     if os.path.exists(UNION) else 0
# B
B, Bmeta = load_csv_keys(DIFF)
# imputed
IMPk, _ = load_csv_keys(IMP)
# C = B ∩ imputed
C = B & IMPk
# D = C - RNA editing
def is_edit(k):
    ra=k.split("_")[-1]  # REF>ALT
    r,a=ra.split(">"); return (r=="A" and a=="G") or (r=="T" and a=="C")
D = {k for k in C if not is_edit(k)}

print(f"\n  A. MuTect2 union sites            {nA:>10}")
print(f"  B. + Fisher differential          {len(B):>10}")
print(f"  C. + Imputation-as-filter (∩)     {len(C):>10}")
print(f"  D. + RNA editing filter (final)   {len(D):>10}")

# DNA truth per-sample sets (15 calls)
print("\n  DNA ground truth 로딩(15 samples)...", flush=True)
dna={s:sample_dna_keys(s) for s in CONTROLS+TREATMENTS}

def classify(keys):
    som=germ=rna=0
    for k in keys:
        t=sum(1 for s in TREATMENTS if k in dna[s])
        c=sum(1 for s in CONTROLS   if k in dna[s])
        if t>0 and c==0: som+=1
        elif t>0 and c>0: germ+=1
        else: rna+=1
    return som,germ,rna

print("\n"+"="*70); print("DNA validation (TP / PPV)"); print("="*70)
print(f"\n{'Stage':<34}{'Total':>8}{'Somatic':>9}{'Germline':>10}{'RNA-only':>10}{'PPV':>8}")
print("-"*79)
rows=[]
for label,ks in [("C. ∩ Imputation",C),("D. + RNA editing (final)",D)]:
    s,g,r=classify(ks); tot=s+g+r
    ppv=f"{100*s/tot:.1f}%" if tot else "N/A"
    print(f"  {label:<32}{tot:>8}{s:>9}{g:>10}{r:>10}{ppv:>8}")
    rows.append((label,tot,s,g,r,ppv))

# strict 4-tuple match (모든 샘플 DNA union)
dna_union=set().union(*dna.values())
print(f"\n  [참고] strict 4-tuple match (전체 DNA union):")
print(f"    C ∩ DNA = {len(C & dna_union)} / {len(C)}")
print(f"    D ∩ DNA = {len(D & dna_union)} / {len(D)}")

# summary 저장 (작은 파일 → commit 가능)
with open(f"{OUT}/ablation_summary.tsv","w") as f:
    f.write("stage\ttotal\tsomatic\tgermline\trna_only\tppv\n")
    f.write(f"A_mutect2_union\t{nA}\t\t\t\t\n")
    f.write(f"B_fisher_differential\t{len(B)}\t\t\t\t\n")
    for label,tot,s,g,r,ppv in rows:
        f.write(f"{label}\t{tot}\t{s}\t{g}\t{r}\t{ppv}\n")
print(f"\n요약 저장: {OUT}/ablation_summary.tsv")
print("Done.")
