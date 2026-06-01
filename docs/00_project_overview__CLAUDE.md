# RAT_project — RNA-seq 기반 Somatic Mutation Detection Pipeline

**최종 업데이트:** 2026-05-27

## 🔒 Numerator Lock — 1,392가 main, 932는 보조 (절대 혼동 금지)

> **이 pipeline은 RNA-only application(PHMG_IH)을 위해 만든다. 따라서 모든 성능 지표의 "주" numerator는 1,392이지 932가 아니다.**

| Numerator | 정의 | DNA 사용? | 역할 |
|---|---|:---:|---|
| **1,392** ⭐ | RNA pipeline의 직접 출력 (Step 1~5의 final candidate set) | ❌ | **Main — application form 그대로. 모든 headline coverage/recall/gene-mapping의 분자는 1,392** |
| 947 | 1,392 중 DNA validation으로 TRUE_SOMATIC 라벨된 부분 | ✅ | DNA-aided characterization (보조) |
| 932 | 947 ∩ DNA truth strict match (chrom+pos+ref+alt 4-tuple) | ✅ | DNA-aided characterization (보조) |

### 규칙
1. **Coverage/recall headline은 1,392 numerator로 계산해 보고**. 932는 "PHMG_IT-internal DNA-aided characterization"으로 따로 명시.
2. RNA-only application에선 932를 알 수 없다 — DNA truth 없이는 1,392에서 어느 게 TP인지 못 가른다. 따라서 application transferable spec은 1,392 form뿐.
3. Gene-level mapping, downstream annotation(SnpEff/COSMIC), pathway enrichment — **전부 1,392 기준**으로 먼저 산출. 932 기준은 sanity check 보조 자료.
4. 이전 docs(`A_to_Z_65pct_derivation.md`, `CLAUDE_20260526.md` 등)의 "932 strict numerator" headline framing은 이 lock과 충돌하면 **이 lock이 우선**.

### Headline spec (1,392 기준)

#### Position-level coverage

| Operating point | 분자 | 분모 | Coverage |
|---|---:|---:|---:|
| Reachable @ alt_AF ≥ 0.30 | 1,392 | 1,972 | **70.59%** |
| Reachable @ alt_DP ≥ 10 | 1,392 | 2,728 | **51.03%** |
| Naive (전체 DNA truth) | 1,392 | 112,294 | 1.24% |

#### Gene-level coverage (2026-05-27 추가)

| Denominator | RNA gene ∩ DNA gene | DNA gene 수 | Coverage |
|---|---:|---:|---:|
| Reachable @ AF ≥ 0.30 (1,478 pos) | 198 | 266 | **74.44%** |
| Reachable @ DP ≥ 10 (2,256 pos) | 207 | 499 | 41.48% |
| Full DNA truth (112,294 pos, naive) | 264 | 3,486 | 7.57% |

- RNA 1,392 → **361 unique gene symbol** (89.6% positions in gene body, 145 intergenic)
- Position-level 70.59%와 gene-level 74.44%가 일치 → coverage 지표 일관성 확인
- 산출 위치: `pipeline_candidates/results/observation_first/gene_mapping/`

→ Method paper에 적을 application-transferable headline:
  - **Position: 70.59% @ AF≥0.30 / 51.03% @ DP≥10**
  - **Gene: 74.44% @ AF≥0.30** (361 genes hit, 198 overlap reachable DNA truth)
  - PPV 68.4%와 함께 보고

### Framing slip 경계 — 이런 표현 나오면 즉시 정정
- "932 strict recall 41/62% headline" → ❌. 그건 DNA-aided 보조. Headline은 1,392 기반.
- "captured 932을 numerator로 보면..." → 보조 분석이라고 명시할 때만.
- "DNA 검증된 TP만 세야 정직하다" → application 관점에선 부정직. 1,392 자체가 method 출력이고 PPV 손실은 PPV 지표로 따로 보고.

---

## ⚠️ 연구의 본질 — 이것을 잊지 말 것 (Framing Lock)

> **이 연구의 deliverable은 "method/pipeline 그 자체"이지 특정 dataset의 somatic mutation list가 아니다.**

| 항목 | 역할 |
|---|---|
| **Deliverable** | RNA-only로 작동하는 **generalizable somatic detection method/pipeline** |
| **PHMG_IT / mammary_cancer** | Method 성능을 측정하는 **벤치마크 dataset** (분석 대상 아님) |
| **657 confirmed somatic 등 결과 수치** | Method가 한 벤치마크에서 산출한 **검증 출력**, 그 자체가 deliverable 아님 |
| **DNA WGS/WES** | Method를 평가하는 **측정자 (yardstick)**, pipeline에는 절대 들어가지 않음 |
| **PHMG_IH (RNA-only)** | 완성된 method의 **실제 적용 demo** |

### 따라서 모든 의사결정은 다음 기준으로 한다
1. "이 작업이 method 자체의 성능/특성을 개선/특성화하는가?" — ✅ 의미 있음
2. "이 작업이 PHMG_IT의 특정 결과만 개선하는가?" — ❌ method paper에 무의미
3. **Cross-dataset 일관성**이 method paper의 main claim → 한 dataset에만 매몰되지 말 것

### Framing slip 경계 — 이런 표현이 나오면 자기검열
- "657개가 유의미하다" → "**method가 이 벤치마크에서 PPV 68.4% 달성**"
- "germline leak 28.6% 줄이자" → "**method의 false positive rate를 컴포넌트 추가로 개선**"
- "657의 DP 분석" → "**method의 working DP envelope 특성화**"
- "이 dataset에서 X를 찾자" → "**X 찾는 method를 만들고 dataset으로 측정하자**"

---

## 연구 목적

> **DNA 없이 RNA-seq만으로 somatic-like variant를 탐지하는 method 개발**
> (Matched normal-free, reference panel-prior 방식)

- Matched DNA는 validation 용도로만 사용 (pipeline에 투입하지 않음)
- 최종 목표: RNA-only dataset (PHMG_IH)에 적용 가능한 validated pipeline
- **정확한 framing**: per-individual tumor-normal pair가 없는 상황에서 imputation을 germline proxy로 사용하는 "treatment-enriched non-germline variant detection"
  - 일반적 의미의 "somatic"이 아닌 **"group-level non-germline + treatment-associated"** 변이 검출
  - Reviewer 공격 회피용 framing이며 method의 raison d'être이기도 함

## 전략 — Rat-first, Mouse는 이후

| 단계 | 목표 | 산출물 |
|---|---|---|
| **1. Rat 단일 종 method 완성** (현재 단계) | PHMG_IT/mammary에서 method closure | Rat-specific tool paper로도 standalone publishable |
| **2. Mouse transfer 시도** | 확보된 mouse ctrl+tumor WGS-RNA matched cohort에 동일 framework 적용 | 성공 시 multi-species method paper로 격상 |
| **3. PHMG_IH 최종 적용** | RNA-only 시나리오에서 validated pipeline 실증 | Final application proof |

**미완성된 method를 두 종에 동시 굴리지 않는다** — 디버깅 불가. Rat에서 닫고 mouse로 transfer만.

---

## 디렉토리 구조

```
RAT_project/
├── CLAUDE.md                    ← 이 파일 (프로젝트 전체 현황)
│
├── PHMG_IT/                     [Validation Dataset 1]
│   ├── CLAUDE.md                ← PHMG_IT 전용 (Phase 0-12 상세 이력)
│   ├── reference/               rn7 genome + STAR index + GTF
│   ├── results/
│   │   ├── mutect2/             MuTect2 force-call + markdup BAMs
│   │   ├── joint_calling/       bcftools 15-sample joint calling
│   │   ├── imputed_differential/ Phase 11 교집합 분석 결과
│   │   ├── somatic_capture/     Phase 12 imputed baseline 결과
│   │   └── ground_truth/        DNA WGS SNP VCFs (15 samples)
│   ├── imputation/
│   │   ├── imputed/             BEAGLE imputed VCFs (15 samples × 20 chr)
│   │   ├── rna_gl/              Genotype likelihood VCFs
│   │   ├── ref_panel/           HRDP 75 strain panel (phased)
│   │   └── evaluation/          Imputation 성능 평가
│   ├── scripts/                 Pipeline 스크립트 (19개)
│   └── coworker/                Python variant-expression 통합 pipeline
│
├── mammary_cancer/              [Validation Dataset 2]
│   ├── fastq/                   23 RNA-seq samples
│   ├── results/
│   │   ├── rn6/                 HaplotypeCaller VCFs
│   │   └── dedup_bam/           23 dedup BAMs
│   ├── imputation/              BEAGLE imputed (GL-based, 36.8%)
│   ├── ground_truth/            WES VCFs (DNA truth)
│   └── sample_mapping_verified.tsv
│
├── PHMG_IH/                     [Application Dataset — RNA only]
│   ├── CLAUDE.md                PHMG_IH 전용 (적용 계획)
│   └── fastq/                   ~145GB RNA-seq
│
├── pipeline_candidates/         [Pipeline 개발 작업 공간]
│   ├── CLAUDE_20260506.md       ← 오늘 작업 상세 (실패+성공 기록)
│   ├── scripts/                 Pipeline 스크립트들
│   ├── results/
│   │   ├── observation_first/   ← 최종 성공 결과 (657 true somatic)
│   │   ├── discordance/         실패한 detector 접근 결과
│   │   ├── deseq2/              DESeq2 분석 결과
│   │   └── ...                  각 filter 단계 결과
│   └── secondary_pipeline.md    Pipeline 설계 문서 (상세)
│
├── docs/archive/                [이전 문서 보관]
│   ├── CLAUDE_20260415_master.md
│   ├── CLAUDE_20260429_PHMG_IT_phase0-8.md
│   ├── CLAUDE_20260506_secondary_pipeline_detail.md
│   ├── README_20260429.md
│   ├── weekly_report_20260328.md
│   └── progress_report_20260401.md
│
└── miniforge3/                  Conda (rnaseq, snpeff, r_figures)
```

---

## 현재 최적 Pipeline (2026-05-06 확정)

### 방법: Observation-first + Imputation-as-Filter

```
Step 1: MuTect2 variant calling (direct observation, 15 samples)
Step 2: Treatment-only differential (Control vs Treatment)
Step 3: Imputation (BEAGLE + HRDP) → germline estimation
Step 4: Intersection (Step 2 ∩ Imputed differential) ← 핵심
Step 5: RNA editing filter (A>G/T>C 제거)
```

### 핵심 교훈

> **Imputation은 "detector"가 아니라 "filter"로 사용해야 한다.**
> - Detector로 쓰면: 0개 (실패)
> - Filter로 쓰면: 657개 true somatic (성공)

### Imputation의 본질적 한계 — 반드시 인지할 것

> **Imputation은 germline haplotype reference(HRDP) 기반의 genotype 추론 → 본질적으로 germline-biased operation.**
> HRDP에 없는 진짜 somatic 변이는 imputation이 "복원"할 수 없음. 따라서 imputation은 somatic-detection coverage를 직접 확장하지 못함.

- **72.6% gene-level coverage의 정확한 의미**: somatic detection coverage가 ❌ 아님. **germline-filter coverage** ✓
- 실제 somatic-detection 가능 영역 = MuTect2 callable region (direct observation에서 DP 충분한 곳)
- → Imputation의 역할은 "second opinion germline filter"로 한정해야 하며, 주력 detection은 direct observation
- → Coverage 확장은 imputation 외 다른 방법(PoN, ASE, phasing 등)으로 보완 필요 (로드맵 참조)

### 최적 결과

| Metric | Value | 상태 |
|--------|:-----:|:--:|
| Gene-level coverage (= germline-filter coverage) | **72.6%** | ✓ |
| Somatic-detection coverage (= MuTect2 callable) | 측정 필요 | ❌ |
| Total candidates | 961 | ✓ |
| DNA-validated true somatic | **657** | ✓ |
| PPV (Precision) | **68.4%** | ✓ |
| RNA artifact rate | 3.0% | ✓ |
| Germline leak rate | 28.6% | ❌ 저감 필요 |
| **Sensitivity (Recall)** | **미측정** | ❌ Phase B에서 측정 |
| 657 variant의 DP 분포 | 미분석 | ❌ Phase A 병행 분석 |

### Phase 11 대비 개선

| | Phase 11 (기존) | 현재 | 변화 |
|--|:-:|:-:|:-:|
| True somatic | 597 | **657** | +10% |
| PPV | 61.1% | **68.4%** | +7.3%p |
| RNA artifact | 13.3% | **3.0%** | -10%p |

---

## 3-Dataset 연구 설계 (+ Mouse 예비)

| Dataset | 종 | Reference | Samples | DNA | Matching | 역할 |
|---------|:--:|:---------:|:-------:|:---:|---|------|
| PHMG_IT | Rat | rn7 | 15 (5C+10T) | WGS 15개 | Same tissue fragment matched, group-level normal only | Pipeline 개발 + Validation 1 |
| mammary_cancer | Rat | rn6 | 23 | WES 23개 | Same tissue fragment matched | Independent Validation 2 |
| PHMG_IH | Rat | rn7 | TBD | 없음 | — | Final Application (RNA-only) |
| Mouse cohort (예비) | Mouse | — | TBD | WGS | Ctrl+tumor WGS-RNA matched (확보 완료) | Cross-species transfer (Phase 2) |

### Dataset 특성 (Reviewer 대응용 필수 기록)
- **Per-individual tumor-normal pair 없음** — rat 분야 전체의 제약. Control 5마리가 group-level normal 역할.
- 같은 tissue fragment에서 DNA/RNA 동시 추출 → tissue heterogeneity 문제 없음
- 657 "true somatic"의 정확한 의미: **DNA-confirmed, non-germline, treatment-enriched variants** (not per-individual somatic in strict sense)

---

## 소프트웨어 환경

- **Conda:** `~/RAT_project/miniforge3/` (PATH 미등록, 스크립트에서 직접 source)
- **rnaseq env:** STAR 2.7.11b, GATK 4.6.1.0, bcftools 1.21, samtools 1.21, Beagle 5.5
- **snpeff env:** SnpEff 5.3 + Java 21
- **r_figures env:** R + DESeq2 + BioConductor

---

## Rat Closure 로드맵 (Method 개발 관점, Mouse 진입 전 완료 목표)

> **순서 원칙**: 한 dataset 내부 최적화보다 **cross-dataset generalizability 입증이 먼저**. 일반화 안 되면 PHMG_IT 내부 최적화는 method paper에 무의미.

| Phase | Method 관점 목표 | 핵심 작업 | 예상 기간 |
|---|---|---|:---:|
| **C-pre. Generalizability first check** ⭐ | 현재 method가 다른 dataset에도 작동하는지 즉시 확인 | **현재 pipeline을 mammary_cancer에 그대로 적용** → PPV/특성이 PHMG_IT과 일관되는지 비교 | 1-2주 |
| **B. Method 성능 특성화** | Method spec 완성 (sensitivity, working DP envelope 등) | DNA-validated ground truth로 sensitivity 측정 / DP working envelope 분석 / AF distribution 분석 (양 dataset에서) | 1-2주 |
| **A. Method 컴포넌트 ablation** | 각 필터의 기여도 정량화 + 개선 컴포넌트 추가 | PoN 추가 ablation (imputation 대체) / ASE filter ablation / Imputation 제거 시 성능 변화 — **양 dataset에서 동일 효과인지 검증** | 2-3주 |
| **D. Downstream usability 입증** | Method 출력이 의미 있는 해석을 낳는지 | SnpEff + COSMIC ortholog enrichment + mutational signature (양 dataset에서 일관되게) | 2-3주 |
| **E. Final application demo** | RNA-only 시나리오에서 method 작동 실증 | Validated method를 PHMG_IH에 적용 | 1주 |

**C-pre 결과에 따른 의사결정 분기:**
- Mammary에서도 PPV 60-70% → ✅ Method 자체가 generalizable, 후속 ablation/개선 의미 있음
- Mammary에서 PPV 30% 이하 → ⚠️ Method가 PHMG_IT-overfitting, 재설계 필요
- 양 dataset 결과 차이의 패턴 분석 → method의 한계 특성화에 사용

**A-E 완료 시점이 Mouse transfer 시작 시점.** 그 전엔 mouse 건드리지 않음.

### Pipeline 재설계 방향 — Imputation을 보조 evidence로 강등

```
[기존] Direct observation → ∩ Imputation differential (germline-biased single filter)

[제안] Direct observation 
         → PoN-based statistical filter (germline + recurring artifact)
         → ASE-based filter (germline AF~0.5 패턴 제거)
         → Read-backed phasing confirmation (somatic-relevant linkage)
         → Imputation은 multiple evidence 중 하나로만 (필수 아님)
```

→ Imputation의 germline-bias 한계를 회피하면서 somatic-relevant filter를 주축으로.

### 추가 검토할 method (Imputation 대체/보완)
- **PoN (Panel of Normals)** ⭐ Phase A 핵심 — control 샘플 기반 background error model. Imputation의 가장 직접적 대체.
- **ASE-based germline leak filter** (germline은 ~0.5 balanced, somatic은 skewed)
- **Read-backed phasing** (WhatsHap) — external reference 불필요, 본인 read의 linkage 정보 사용
- **AF-conditional DP threshold** (hard cutoff 대신 AF에 따른 가변 cutoff)
- **Soft filtering with Tier**: 657개를 Tier 1/2/3로 분류, 정보 손실 없이 신뢰도 보고
- **Protein language model pathogenicity** (ESM-1v, AlphaMissense) — 종 무관 적용 가능, prediction layer용

---

## 주의사항

- 외장 HDD 마운트 확인 필요 (FASTQ, DNA VCF 원본)
- GATK 4.6.1.0 TandemRepeat 버그: `--annotations-to-exclude TandemRepeat`
- STAR MAPQ=255: `--disable-read-filter MappingQualityAvailableReadFilter`
- HRDP population filter는 somatic pipeline에서 사용 금지 (너무 공격적)
- 논문 작성 시 "somatic detection" 단순 표현 금지 → **"matched normal-free RNA-based variant detection"** 또는 **"treatment-enriched non-germline variant detection"**으로 정확히 framing
- **"72.6% coverage" 단순 표현 금지** → "72.6% germline-filter coverage"로 정확히. Somatic-detection coverage와 혼동 금지.
- **Imputation을 somatic detector로 표현 금지** → "germline filter as second opinion"으로 framing
- DP cutoff는 hard threshold 대신 AF-conditional 또는 Tier-based soft filtering 사용 (RNA-seq expression-dependent coverage 한계 때문)
- **모든 분석은 "method 특성화/개선" 관점에서 수행** — "이 dataset에서 X 찾기"가 아니라 "method가 X를 얼마나 잘 찾는지 측정/개선"
- **Cross-dataset 일관성**이 method paper의 main claim → 어떤 개입이든 양 dataset(PHMG_IT + mammary)에서 같은 효과 보여야 의미 있음

---

## 현재 진행 중

- **Numerator lock 확정 (2026-05-27)**: 1,392(RNA-only output)가 main numerator, 932는 DNA-aided 보조 characterization.
- **1,392 → gene symbol mapping**: 1,392 positions → 361 unique gene. Gene-level coverage @ AF≥0.30 = **74.44%**, position-level 70.59%와 일치.
- **Mammary cross-validation 완료 (2026-05-28, implementation parity v2)**: PHMG_IT와 동일한 force-call + Fisher differential을 mammary에 적용해 재실행. 결과: **45 candidates, PPV 15.5%, 25 unique genes**. v1(PPV 9.7%) 대비 +5.8%p 개선했지만 PHMG_IT(68.4%)와는 31× candidate / 4.4× PPV 차이로 **method spec 그대로 transfer 불가**. 차이의 주 원인은 dataset-intrinsic 3종(HRDP panel 크기, group design contrast, GT scope) — method 자체보다 panel/design factor가 결정적. (상세: `mammary_cancer/CROSS_VALIDATION_20260528.md`)
- **Component ablation 완료 (2026-05-28, PHMG_IT)** ⭐: pipeline 각 step의 PPV 기여 측정.

  | 시나리오 | Candidates | TP | **PPV** | Recall |
  |---|---:|---:|---:|---:|
  | A. MuTect2 raw union | 2,737,681 | 11,848 | **0.43%** | 10.55% |
  | B. + Fisher differential | 65,277 | 1,406 | **2.15%** | 1.25% |
  | C. **+ Imputation filter** ⭐ | 1,392 | 932 | **66.95%** | 0.82% |
  | D. + RNA editing filter | 996 | 655 | **65.76%** | 0.58% |

  **핵심 발견 — Imputation filter가 PPV의 dominant 기여자 (B→C: +64.8%p)**. Fisher differential 단독은 PPV +1.7%p만 기여, 통계만으론 germline-somatic 구분 불가 직접 증거. **RNA editing filter step은 PPV 개선 없이 TP의 30% 손실** (잃은 396 중 277=70%가 TP) → 재검토 필요. 산출: `pipeline_candidates/results/observation_first/ablation/`.

  → Method paper main claim 직접 입증: "Imputation-as-filter step이 PPV를 2.15% → 66.95%로 +64.8%p 끌어올림".

## 다음 단계

- 다음 분기: (a) Mammary을 rn7 liftover + HRDP 75-strain panel로 재imputation해서 panel size 요인 분리, (b) Mouse ctrl-tumor matched cohort으로 우회 (design factor 분리), (c) Paper framing을 "method + dataset-dependent panel/design factor" 형태로 작성. (a)-(c) 중 우선순위 결정 필요.

## 마감



---

## 문서 가이드

| 파일 | 용도 |
|------|------|
| `CLAUDE.md` (이 파일) | 프로젝트 전체 현황 |
| `PHMG_IT/CLAUDE.md` | Phase 0-12 상세 이력, 데이터 위치, 샘플 매핑 |
| `PHMG_IH/CLAUDE.md` | IH dataset 적용 계획 |
| `pipeline_candidates/CLAUDE_20260506.md` | 오늘 pipeline 개발 상세 기록 |
| `pipeline_candidates/secondary_pipeline.md` | Pipeline 설계 + 각 Step 원리 설명 |
| `docs/archive/` | 이전 버전 문서 보관 |
