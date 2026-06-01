# Secondary Pipeline Candidate: RNA-seq Somatic Mutation Detection

**작성일:** 2026-05-06  
**목적:** DNA 없이 RNA-seq만으로 somatic mutation을 탐지하는 method 개발  
**핵심 claim:** 72% gene-level coverage 달성 + treatment-induced somatic mutation detection

---

## 연구 구조

```
Validation datasets (DNA + RNA matched):
  - PHMG_IT: 15 samples (5 Control, 10 Treatment), matched DNA WGS
  - mammary_cancer: 23 samples, matched WES (NMU-induced)

Application dataset (RNA only):
  - PHMG_IH: RNA-seq only, Control + Treatment 구조
```

---

## Core Pipeline (RNA-seq only로 작동)

### Step 1: RNA-seq Variant Calling

**하는 일:** RNA-seq reads에서 reference genome과 다른 position을 찾아냄

**원리:**
- Reference와 비교하여 ALT allele이 관찰되는 위치를 기록
- 이 시점에서는 germline/somatic 구분 없이 "관찰된 variant" 확보
- Tool: bcftools mpileup / GATK HaplotypeCaller / MuTect2

**한계:** RNA-seq이므로 발현되는 유전자 영역에서만 관찰 가능 (~56% gene-level)

---

### Step 2: Genotype Imputation (GLIMPSE2 + HRDP panel)

**하는 일:** 각 sample의 "expected germline genotype"을 population reference panel 기반으로 추정

**원리 (Linkage Disequilibrium):**
- 유전체에서 물리적으로 가까운 variant들은 함께 유전됨 (haplotype block)
- 일부 site만 관찰해도 LD pattern으로 나머지를 추정 가능
- HRDP 75 rat strains의 germline haplotype을 reference로 사용

**핵심 가정:** Reference panel = germline → imputation 결과 = per-sample germline estimation

**효과:** 
- Gene-level coverage: 56% → **72%** 상승
- Low-depth (DP 1-4) site에서도 GL 기반 imputation 가능
- 각 site에 Genotype Probability (GP) score 부여 → confidence 정량화

**Imputation의 역할 정의:**
- "somatic을 찾는 도구"가 아님
- "germline이 뭔지 추정하는 도구" → somatic은 그 반대로 정의

---

### Step 3: Discordance Detection (Unbiased)

**하는 일:** Step 1 (실제 관찰)과 Step 2 (germline 추정)를 비교

**원리:**
```
Imputed germline = 0|0 (REF)  +  RNA-seq observed = ALT reads
→ Discordance: "germline은 REF인데 왜 ALT가 보이지?"
→ 가능한 원인: somatic mutation / imputation error / RNA editing / artifact
```

**적용 범위:** Imputed 72% gene coverage 전체에서 수행 (unbiased, 특정 DB 한정 아님)

---

### Step 4: Multi-sample Differential Filter

**하는 일:** Discordance의 원인을 분류 — somatic vs imputation error 구분

**핵심 논리:**
```
Imputation error → Control과 Treatment에서 비슷한 빈도로 발생 (systematic)
Somatic mutation → Treatment에서만 발생 (treatment-specific)
```

**적용:**
- Control에서도 같은 discordance → imputation error → 제거
- Treatment에서만 discordance → somatic candidate 유지
- Recurrence ≥2 in Treatment → random error가 아닌 systematic event

---

### Step 5: Expression-aware Filter (DESeq2 기반)

**하는 일:** 발현 변화로 인해 "새로 보이는" germline variant를 제거

**원리:**
```
문제: Gene X가 Control에서 비발현 → Treatment에서 발현 시작
     → Treatment에서만 variant 관찰 → "somatic?"
     → 실은 germline variant가 발현 변화로 새로 관찰된 것

해결: DESeq2로 Treatment에서 newly expressed된 유전자 식별
     → 해당 유전자의 "Treatment-only" variant 제거 또는 flag
     → Imputed germline과 일치하면 "germline이 새로 보인 것" 확정
```

---

### Step 6: Artifact Filtering

**하는 일:** Biological/technical artifact 제거

| Filter | 원리 |
|--------|------|
| RNA editing DB (RADAR/REDIportal) | DNA 변이 없이 RNA에서만 A→I(G) 변환. 알려진 site 제거 |
| Strand bias | 진짜 variant면 양 strand에서 보여야 함. 한 strand만 → artifact |
| DP filter | Too low depth → sequencing error와 구별 불가 |
| AF filter | Too high (>0.9) → germline 가능성. Too low (<0.05) → noise 가능성 |

---

### Step 7: Population AF + GP Confidence Filter

**하는 일:** 개체 간 germline 차이를 somatic으로 오인하는 것 방지

**적용:**
- HRDP 75 strain 중 1개라도 해당 variant 보유 → population에 존재하는 germline → 제거
- GP(0|0) < threshold (e.g., 0.9) → imputation 불확실 → 제거 또는 low-confidence flag
- DNA validation set에서 GP threshold 최적화

---

### Step 8: Known Somatic DB Annotation (Post-hoc)

**하는 일:** Detection 후 surviving candidates를 cancer DB와 대조

**원리:**
- Detection은 unbiased (Step 3-7에서 DB 참조 없이 수행)
- Annotation은 post-hoc: 발견된 것의 biological meaning 부여
- Cancer type enrichment: Fisher's exact test → OR, p-value

**적용:**
- COSMIC/ClinVar/TCGA 대조 → 어떤 cancer type과 overlap?
- Enrichment test (DB size normalize): "PHMG somatic are enriched in lung cancer sites"
- DB hit = confidence 가중치 ↑ (imputation error가 우연히 COSMIC site에 떨어질 확률 극히 낮음)

---

### Step 9: Mutational Signature Analysis

**하는 일:** PHMG-specific mutational signature 추출 및 matching

**원리:**
- 모든 mutagen은 고유한 trinucleotide context 패턴을 남김 (5'-N[mutation]N-3')
- DNA-confirmed somatic에서 PHMG signature 추출 (최초 규명)
- RNA pipeline candidates와 signature 일치 여부 확인

**적용:**
- Tool: SigProfiler / MutationalPatterns (R)
- PHMG signature와 일치 → PHMG-induced confidence ↑
- NMU signature (mammary_cancer)와 비교 → mutagen-specific 구별 가능

---

## Validation Layer (DNA matched datasets 활용)

### A. DNA Ground Truth Validation

```
PHMG_IT:
  DNA WGS에서 직접 somatic calling (Treatment DNA vs Control DNA)
  → "DNA-confirmed somatic list" 확보
  → RNA pipeline 결과와 비교 → Sensitivity, PPV 산출

mammary_cancer:
  WES에서 somatic calling
  → Independent validation → same metrics
```

### B. Permutation Test

```
Treatment/Control label 무작위 섞기 (1000회)
→ 매번 pipeline 실행 → null distribution 생성
→ Observed vs Permuted 비교 → p-value
→ "Treatment label에 의존하는 진짜 signal 존재" 증명
```

### C. Cross-dataset Validation

```
Strategy A: PHMG_IT에서 threshold 최적화 → mammary_cancer에서 test
Strategy B: 반대 방향
→ 양방향 성능 유지 = overfitting 아님, method generalizable
→ 다른 mutagen (PHMG vs NMU), 다른 tissue에서도 작동
```

### D. Signature Extraction

```
DNA-confirmed somatic mutations에서 mutational signature 추출
→ PHMG signature 최초 규명
→ RNA-only candidates에 signature matching 적용 시 성능 향상 확인
```

---

## Confidence Scoring (최종 output)

```
High confidence:   Discordance + Treatment-only + COSMIC hit + Signature match
Medium confidence: Discordance + Treatment-only + (COSMIC OR Signature)
Low confidence:    Discordance + Treatment-only only
Novel candidate:   Discordance + Treatment-only + No DB hit (PHMG-specific novel mutation)
```

---

## Reviewer 방어 구조

| 공격 포인트 | 방어 | 해결 Step |
|------------|------|-----------|
| Imputation error → somatic 오인 | Treatment-only filter + GP score | Step 4, 7 |
| RNA editing → somatic 오인 | RADAR/REDIportal DB 제거 | Step 6 |
| Expression 변화 → 새 variant 관찰 | DESeq2 filter | Step 5 |
| 개체 간 germline 차이 | Population AF + imputation GP | Step 7 |
| 순환 논리 (DB-guided detection) | Unbiased detection → post-hoc annotation | Step 3-7 vs Step 8 |
| Coverage 주장의 근거 | DNA ground truth로 검증 | Validation A |
| Pipeline FP rate 보증 | Permutation + cross-dataset | Validation B, C |
| PHMG vs background 구별 | Mutational signature | Step 9, Validation D |
| Low-depth site 신뢰도 | Multi-evidence + depth tier 보고 | Step 8-9 |

---

## 기술적 한계 (Limitation, 논리적 결함 아님)

1. **비발현 유전자** — RNA-seq 본질적 한계, 어떤 방법으로도 관찰 불가
2. **Very low AF somatic (<0.05)** — Read depth 한계로 탐지 어려움
3. **Recurrence ≥2 filter** — Sample-private somatic mutation 놓침 (sensitivity 손실)
4. **Human→Rat ortholog mapping** — COSMIC annotation 시 불완전한 mapping
5. **Imputation precision (~75-82%)** — GP threshold로 통제하나 완벽하지 않음

---

## 논문 구조 매핑

| Section | 내용 |
|---------|------|
| Introduction | RNA-seq based somatic detection의 필요성, 기존 한계 |
| Method | Core Pipeline (Step 1-9) |
| Results - Validation 1 | PHMG_IT: sensitivity/PPV, signature, permutation |
| Results - Validation 2 | mammary_cancer: independent validation, cross-dataset |
| Results - Application | PHMG_IH: RNA-only 적용 결과 |
| Discussion | COSMIC enrichment, cancer type correlation, PHMG signature 해석 |
| Limitation | 기술적 한계 명시 |

---

## Pipeline Claim (한 문장)

> "RNA-seq 단독으로 72% gene-level coverage를 달성하면서, imputation-based germline estimation과 multi-sample differential analysis, mutational signature matching을 결합하여 treatment-induced somatic mutation을 탐지하고, known cancer somatic database와의 enrichment를 통해 화학물질의 발암 메커니즘을 추론하는 통합 pipeline"

---

## 기존 방법과의 차이

| 항목 | Phase 7 (기존) | Phase 12 (기존) | Secondary Pipeline (본 문서) |
|------|---------------|----------------|--------------------------|
| Coverage | 8% (variant site) | 72% (gene) | 72% (gene) |
| Somatic logic | Treatment-only observed | Imputed vs observed | Discordance + Treatment-only + multi-evidence |
| FP control | 없음 | 약함 (85% FP) | 5-layer filtering |
| Validation | DNA 대조 | DNA 대조 | DNA + permutation + cross-dataset |
| Signature | 없음 | 없음 | PHMG signature 추출 + matching |
| DB 활용 | SnpEff annotation | 없음 | Post-hoc enrichment (non-circular) |
| Reviewer 방어 | 약함 | 약함 | 강함 (모든 공격 대비됨) |

---

## PHMG_IT 실행 결과 (2026-05-06)

### 실행 환경
- Dataset: PHMG_IT 15 samples (C1-C5 Control, P1-P10 Treatment)
- 기존 자산 재활용: MuTect2 PASS VCF, BEAGLE imputed VCF, DNA ground truth VCF, HRDP panel
- 스크립트: `/home/yusanghyeon/RAT_project/pipeline_candidates/scripts/`
- 결과: `/home/yusanghyeon/RAT_project/pipeline_candidates/results/`

### Step 3: Discordance Detection 결과

| Sample | Observed PASS (DP≥5, AF≥0.05) | Imputed REF (0\|0) | Discordances |
|--------|:---:|:---:|:---:|
| C1 | 86,151 | 6,912,320 | 344 |
| C2 | 107,234 | 7,029,011 | 335 |
| C3 | 101,039 | 6,969,481 | 293 |
| C4 | 104,103 | 7,046,796 | 283 |
| C5 | 92,328 | 7,051,148 | 275 |
| P1 | 80,443 | 7,056,045 | 372 |
| P2 | 98,999 | 7,054,814 | 230 |
| P3 | 72,920 | 7,049,055 | 190 |
| P4 | 69,632 | 6,938,573 | 276 |
| P5 | 82,867 | 6,975,569 | 253 |
| P6 | 70,842 | 7,001,742 | 203 |
| P7 | 91,589 | 7,027,639 | 392 |
| P8 | 55,943 | 7,072,692 | 224 |
| P9 | 92,704 | 7,048,993 | 486 |
| P10 | 80,881 | 7,072,671 | 305 |
| **Control 평균** | 98,171 | 7,001,751 | **306** |
| **Treatment 평균** | 79,682 | 7,030,379 | **293** |

**관찰:** Control과 Treatment 간 discordance 수가 유사 → discordance 대부분은 imputation error (group 비특이적)

### Step 4: Multi-sample Differential Filter 결과

| 단계 | Count |
|------|:---:|
| Total Treatment discordance sites (unique) | 2,734 |
| Control overlap 제거 | -161 |
| Treatment-only sites | 2,573 |
| **Recurrence ≥ 1** | **2,573** |
| **Recurrence ≥ 2** | **117** |
| **Recurrence ≥ 3** | **12** |

### Step 5: DESeq2 Expression Analysis 결과

| 항목 | 수치 |
|------|:---:|
| 분석 유전자 수 (count ≥ 10) | 25,143 |
| Significantly DE genes (padj < 0.05) | 5,048 |
| Newly expressed in Treatment (Control<5, Treatment≥10) | 894 |

*Note: DESeq2 expression filter는 5개 최종 후보에 대해서는 적용 전이나, 전체 파이프라인 프레임워크에서 활용 가능*

### Step 6: Artifact Filter 결과

| Filter | Input | Removed | Output |
|--------|:---:|:---:|:---:|
| RNA editing (A>G / T>C) | 117 | 37 | 80 |
| AF filter (0.05-0.9) | 80 | 1 | **79** |

Mutation type distribution (79개):
- G>A: 24, C>T: 20, G>T: 6, C>A: 6, G>C: 6, A>T: 5, A>C: 4, T>A: 2, C>G: 1, INDEL: 5

### Step 7: Population AF + GP Confidence Filter 결과

| Filter | Input | Removed | Output |
|--------|:---:|:---:|:---:|
| HRDP population (11.3M variants) | 79 | 73 | 6 |
| GP/DS confidence (DS > 0.3) | 6 | 1 | **5** |

**최종 5개 후보:**

| chrom | pos | ref | alt | recurrence | samples | avg_AF | avg_DS |
|-------|-----|-----|-----|:---:|---------|:---:|:---:|
| chr13 | 23242612 | G | GT | 3 | P7,P9,P10 | 0.201 | 0.000 |
| chr15 | 20626488 | long INDEL | A | 2 | P9,P10 | 0.077 | 0.020 |
| chr19 | 22961825 | G | T | 3 | P1,P8,P10 | 0.384 | 0.047 |
| chr20 | 3315342 | TAAAA | T | 3 | P2,P3,P7 | 0.392 | 0.053 |
| chr7 | 107110535 | G | A | 2 | P3,P5 | 0.484 | 0.090 |

### DNA Ground Truth Validation 결과

| Stage | Total | True Somatic | Germline Leak | RNA-only | PPV |
|-------|:---:|:---:|:---:|:---:|:---:|
| Step 4 (rec≥2) | 117 | **0 (0%)** | 98 (83.8%) | 19 (16.2%) | **0%** |
| Step 6 (artifact) | 79 | **0 (0%)** | 66 (83.5%) | 13 (16.5%) | **0%** |
| Step 7 (pop+GP) | 5 | **0 (0%)** | 0 (0%) | 5 (100%) | **0%** |

### 결론: Pipeline 실패 원인 분석

**True somatic = 0. 이 pipeline은 somatic mutation을 탐지하지 못했다.**

#### 실패 원인

1. **Imputation precision 한계 (75-82%)**
   - 18-25%의 위치에서 imputation이 틀림
   - 이 error가 개체마다 다른 위치에서 발생 → 우연히 Treatment에서만 error가 발생하는 site 존재
   - "Treatment-only discordance"의 83.8%가 실제로는 germline leak (imputation이 Control에서는 맞추고 Treatment에서는 틀린 것)

2. **Discordance 기반 접근의 근본적 한계**
   - Imputation error rate > Somatic mutation rate
   - 진짜 somatic signal이 imputation error noise에 묻힘
   - Multi-sample differential filter로도 분리 불가 (error도 개체별로 다르므로)

3. **Coverage 72%는 달성했으나 precision이 핵심**
   - Coverage가 높아도 imputation이 틀리면 discordance가 somatic인지 error인지 구분 불가
   - Precision이 ~100%에 가까워야 이 접근이 작동

#### 시사점

Imputation-discordance 기반 somatic detection은 **현재 imputation precision 수준에서는 작동하지 않음.**
대안:
- Phase 7/9의 direct observation 기반 differential analysis가 이미 더 나은 결과 (597 true somatic)
- Imputation의 역할을 "somatic detection"이 아닌 "germline characterization"으로 한정해야 할 수 있음
- 또는 imputation precision을 95%+ 로 올릴 수 있는 방법 필요 (GLIMPSE2, larger panel 등)
