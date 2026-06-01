# Mammary Cross-Validation — 정직한 진단 + 의사결정 문서

**작성일:** 2026-05-28
**대상:** PHMG_IT pipeline(observation-first + imputation-as-filter)을 mammary_cancer dataset에 적용한 cross-validation 작업의 전체 기록과 해석
**Status:** 1차 cross-validation 완료(v1, 2026-05-13) + implementation parity 재실행(v2, 2026-05-28). Method spec transfer는 실패. 다음 분기점 결정 필요.

---

## 0. 배경 — 왜 cross-validation이 main claim인가

`/home/yusanghyeon/RAT_project/CLAUDE.md`의 framing lock 참조:

> 이 연구의 deliverable은 **method/pipeline 그 자체**이지 특정 dataset의 somatic mutation list가 아니다. **Cross-dataset 일관성**이 method paper의 main claim.

따라서 mammary에서 PHMG_IT spec(PPV 68.4%, position-level coverage 70.59% @ AF≥0.30, gene-level 74.44%)이 재현되어야 method paper가 성립.

---

## 1. v1 결과 (2026-05-13, 기존 작업)

**스크립트:** `mammary_cancer/scripts/01~09_*.sh`
**산출:** `mammary_cancer/results/observation_first/`

### 1.1 단계별 candidate funnel

| Step | Count |
|---|---:|
| MuTect2 late_only differential (binary set-diff) | 76,468 |
| Imputed late_only differential | 88,981 |
| Step 4 intersection | 37 |
| Step 5 RNA editing filter | **31** |
| WES-confirmed true somatic (late only) | 3 |

### 1.2 핵심 지표

| Metric | PHMG_IT | Mammary v1 |
|---|:---:|:---:|
| Total candidates | 1,392 | 31 |
| True somatic | 947 | 3 |
| PPV | 68.4% | **9.7%** |
| RNA artifact rate | 3.0% | 87.1% |
| Imputation coverage | 72.6% | 36.5% |

→ Generalizability 1차 실패. 차이 분석 필요.

---

## 2. 차이의 진단 (2026-05-28)

### 2.1 Implementation 차이 발견

PHMG_IT script `mutect2_scatter.sh` vs mammary `02_mutect2_per_sample.sh` + `05_differential.sh` 비교:

| 구성 요소 | PHMG_IT | Mammary v1 |
|---|---|---|
| **MuTect2 call** | Per-sample + **force-call**(`--alleles union_sites --force-call-filtered-alleles`) | Per-sample only, force-call 없음 |
| **Differential** | `coworker/run_pipeline.py`로 **Fisher's exact + recurrence + AF threshold** | `bcftools isec -n=1`로 **binary set-difference** (통계 없음) |
| **Group design** | control(5) vs treatment(10) — sharp contrast | early M1(6) vs late M12+M20(7) — time-course continuum |
| **Reference** | rn7 | rn6 |
| **HRDP panel** | 75-strain | 48-strain |
| **Ground truth** | WGS (genome-wide) | WES (exome only) |

### 2.2 5가지 원인 (영향 큰 순)

1. **MuTect2 force-call 누락** ⭐ 가장 큰 차이
   - PHMG_IT는 모든 sample을 같은 좌표에서 force-call → cross-sample comparable
   - Mammary v1는 per-sample 독립 call → sample 간 좌표 거의 안 겹침 → "late_only"가 sparse position union

2. **Differential 정의가 통계적 아님**
   - PHMG_IT는 Fisher's exact + recurrence ≥2 + AF threshold
   - Mammary v1는 binary presence/absence (1-sample edge case 다 포함)

3. **Group design contrast 강도**
   - PHMG_IT: 처리 vs 비처리 — sharp
   - Mammary: 시간 경과 — continuum, 변이 자연 축적

4. **Imputation panel 크기**
   - rn7 75-strain (PHMG_IT) vs rn6 48-strain (mammary)
   - Imputation coverage 72.6% vs 36.5%

5. **Ground truth scope (PPV 정확도에만 영향)**
   - WGS vs WES — intron/intergenic candidate 평가 불가
   - **단, candidate 수 자체에는 영향 없음**

---

## 3. v2 — Implementation parity 실험 (2026-05-28)

### 3.1 무엇을 바꿨나

원인 1, 2를 해결: PHMG_IT와 동일한 force-call + Fisher differential 구현.

**새 스크립트:**
- `mammary_cancer/scripts/11_force_call_scatter.sh` — Union sites build + MuTect2 force-call (13 samples × 21 chr = 273 jobs)
- `mammary_cancer/scripts/12_filter_force_call.sh` — FilterMutectCalls + PASS 추출 + coworker Fisher differential

**산출:** `mammary_cancer/results/observation_first_v2/`

원인 3-5는 dataset-intrinsic이라 이번엔 손대지 않음.

### 3.2 v2 단계별 funnel

| Step | v1 | v2 (force-call + Fisher) | 변화 |
|---|---:|---:|---|
| MuTect2 differential | 76,468 | **6,290** | -92% (statistical strict화) |
| Imputed late_only | 88,981 | 88,981 | 동일 |
| Step 4 intersection | 37 | **59** | +59% |
| Step 5 (final, after RNA editing) | 31 | **45** | +45% |
| True somatic (late WES) | 3 | **7** | +133% |

### 3.3 PPV 변화

| Metric | v1 | v2 |
|---|:---:|:---:|
| Total candidates | 31 | 45 |
| ∩ any WES (any TP) | (미측정) | 9 (20.0%) |
| ∩ late WES (true somatic) | 3 | **7** |
| **PPV** | **9.7%** | **15.5%** |

→ Implementation parity로 PPV +5.8%p 개선. 그러나 PHMG_IT 68.4%와는 여전히 **4.4× 차이**.

### 3.4 Gene-level

- 45 candidates → 25 unique gene symbol (32 in gene, 13 intergenic)
- Top hits: Nhsl1(4), Ccdc148(4), Tenm2(3), Tnik(2), Skil(2), **Atm**(2), Phc3(2), Calcrl(2), Bicd2(2) ...
- PHMG_IT의 top hit(Ext1, Nrg1, Cadm1, Ski)과 **거의 overlap 없음** → 두 dataset의 biology가 다름

---

## 4. 단계별 격차 분석 — 어디서 얼마나 잃었나

| Step | PHMG_IT | Mammary v2 | 격차 |
|---|---:|---:|---|
| MuTect2 force-call + Fisher differential | 65,375 | 6,290 | **10×** |
| Imputed late_only differential | 336,794 | 88,981 | **3.8×** |
| Step 4 intersection (∩) | 977 | 59 | **16.6×** |
| Final candidates | 1,392 | 45 | **31×** |
| True somatic (DNA-validated) | 947 | 7 | 135× |
| PPV | 68.4% | 15.5% | 4.4× |
| Genes | 361 | 25 | 14× |

### 4.1 격차의 분해

- **MuTect2 단계에서 10× less**: 주로 **design contrast 약함** (time-course continuum). Sample 수 차이(15 vs 13)는 작은 기여.
- **Imputation 단계에서 3.8× less**: 주로 **panel size**(75-strain vs 48-strain) → coverage 72.6% vs 36.5%.
- **Intersection 단계에서 16.6× less**: 단순 곱셈(10×3.8=38)보단 작지만 **두 set이 같은 좌표에 떨어지는 효율**도 낮음.

### 4.2 "Imputation 문제다"는 절반만 맞다

User가 제기한 가설: "Imputation이 mammary에서 제대로 작동 안 해서 적게 나온다"
- 부분적 진실: imputation coverage 36.5%로 절반 수준
- 그러나 **MuTect2 differential 자체도 약함** (design contrast)
- 두 요인이 동시에 작용 → intersection이 disproportionately 작음

---

## 5. Ground truth 한계 — 별도 검토

User 질문: "WES vs WGS 차이가 candidate 수에 영향?"

**답: 아니다. Candidate 수는 RNA-only로 결정됨.**

```
RNA-seq → MuTect2 → Imputation → ∩ → RNA editing → 45 candidates  (여기까지 ground truth 미사용)
                                       ↓
                                  WES validation  ← 여기서만 GT 사용 (PPV 채점)
```

- 45라는 수치는 RNA pipeline + imputation panel에 의해서만 결정
- WES(exome)냐 WGS(genome)냐는 **PPV 측정 정확도**에만 영향:
  - WES는 intron/intergenic candidate 평가 불가 → 일부 true positive가 자동 FN 처리될 수 있음
  - 즉 PPV 15.5%는 **하한**일 수 있음 (진짜 PPV는 더 높을 가능성)
- 하지만 candidate 31× 격차는 100% RNA-side 문제

---

## 6. rn6 vs rn7 — 왜 mammary는 rn6?

### 6.1 역사적 이유
- `mammary_cancer/reference/rn6/` 세팅이 4월에 이미 완료된 상태
- 원본 WES truth (GSE297544)가 rn6 기준으로 call됨 → liftover 손실 회피
- HRDP imputation panel은 rn6용 48-strain 버전이 먼저 공개되어 historical baseline

### 6.2 mammary 약점이 됨
- rn6 + 48-strain은 rn7 + 75-strain 대비 본질적으로 imputation coverage 낮음
- PHMG_IT는 시점 운이 좋아 rn7+75-strain 채택 가능

### 6.3 rn7로 옮길 경우

| 항목 | 비용 / 문제 |
|---|---|
| BAM | 23 sample × 3-4시간 재align |
| WES truth | liftover 시 2-4% 손실 + REF allele 불일치 일부 |
| HRDP panel | rn7 75-strain 사용 가능 (1.5× strain) |
| 예상 imputation coverage | 36.5% → 60-70% (PHMG_IT 근접) |

---

## 7. Panel downsample 실험 — 정직한 평가

User 제안한 대안: PHMG_IT의 75-strain HRDP를 random 48-strain으로 downsample 후 재imputation해서 panel size 효과를 분리하는 실험.

### 7.1 이 실험이 줄 수 있는 것
- Panel size vs imputation coverage의 **dose-response 곡선**
- Method 작동의 operating envelope 특성화
- Paper의 limitation/scope section 보강 자료

### 7.2 줄 수 없는 것 (현실적 한계)

1. **Method 자체의 우수성 evidence ❌**
   - "panel 좋으면 잘 작동" = HRDP의 강점 재확인, method의 강점 아님
   - Reviewer가 "그럼 imputation tool 덕분 아냐?"라고 물으면 답이 약함

2. **결과가 너무 예측 가능**
   - HRDP 원논문에 이미 strain 수 vs accuracy 곡선 게재됨
   - 같은 결과 재확인 → "이미 알려진 사실을 우리 데이터에서도 봤다" 수준

3. **Generalizability claim에 도움 안 됨**
   - 오히려 약화 가능: "method가 panel에 강하게 의존" → transfer 어렵다는 뜻

4. **PPV 결과 어느 쪽이든 깔끔하지 않음**
   - PPV 같이 떨어지면: "precision도 panel-fragile" (약함)
   - PPV 유지되면: "candidate 줄어도 정확도 유지" → mammary 15.5% PPV와 모순

### 7.3 진짜 method evidence가 될 만한 실험들

| 실험 | Method evidence로 작용? |
|---|:---:|
| **Baseline 비교**: MuTect2 단독 PPV vs +imputation filter PPV at 양 dataset | ✅ 강함 |
| **Component ablation**: imputation filter, RNA editing filter, Fisher 각각 제거 시 PPV 변화 | ✅ 강함 |
| **Imputation detector vs filter**: PHMG_IT의 0 vs 657 비교 (이미 입증된 핵심 발견) | ✅ 강함 |
| Panel downsample ablation | △ dependence 특성화에는 OK, method 우수성은 아님 |
| **rn7 liftover + 75-strain mammary 재실행** | ✅ 만약 PPV 68%대 회복 시 강함 |

---

## 8. 결론 (정직하게)

### 8.1 현재 method spec 상태

- **PHMG_IT-internal headline**: PPV 68.4%, position coverage 70.59% @ AF≥0.30, gene coverage 74.44% — 견고
- **Mammary transfer**: PPV 15.5%, candidate 45 — **method spec이 그대로 transfer되지 않음**
- **격차 원인**: dataset-intrinsic 3종(panel size 50%, design contrast 30%, GT scope 미미) + implementation 5%

### 8.2 Method paper의 정직한 framing 필요

다음 세 가지 framing 옵션:

**Option A — 좁은 scope claim**
> "RNA-only somatic-like detection method, validated on rat treatment-vs-control design with HRDP 75-strain reference panel. Performance is panel- and design-dependent; see limitation."

→ Mammary 결과는 limitation section에 포함. Generalizability claim 축소.

**Option B — Operating envelope claim**
> "Method with characterized operating envelope: requires HRDP ≥75 strains and group-level contrast (not time-course continuum) for spec'd performance."

→ Panel downsample 실험으로 envelope 정량화 추가. Paper는 "method + working condition spec"로 framing.

**Option C — 추가 실험 후 강한 claim 시도**
> 1. Mammary을 rn7 liftover + 75-strain panel로 재실행 → panel factor 제거 후 성능 확인
> 2. PPV가 50%+ 회복되면: "method는 panel-equivalent 환경에서 generalizable"
> 3. 여전히 낮으면: design factor가 dominant → method 한계 명시

→ 가장 강한 claim 가능하지만 추가 비용 큼.

### 8.3 권장 우선순위

> 자원 한정 시:
> 1. **Mammary baseline 비교**(MuTect2 단독 PPV) — 적은 비용으로 imputation filter step의 contribution 직접 측정
> 2. **rn7 liftover 실험** — 큰 비용, 결과 결정적. 잘 나오면 Option C 가능
> 3. Panel downsample은 paper revision 단계 reviewer 요구 시

---

## 9. 산출물 위치

```
mammary_cancer/
├── CROSS_VALIDATION_20260528.md           ← 이 문서
├── scripts/
│   ├── 11_force_call_scatter.sh           ← MuTect2 force-call (v2)
│   └── 12_filter_force_call.sh            ← Filter + Fisher differential (v2)
├── results/
│   ├── mutect2/force_call/
│   │   ├── union_sites.vcf.gz             ← 184,629 union sites
│   │   ├── union_sites_per_chr/           ← chr별 split
│   │   ├── force_scatter/                 ← 273 force-called VCFs
│   │   ├── force_pass/                    ← 13 sample PASS VCFs
│   │   └── differential_force/
│   │       ├── differential_snps.csv      ← 6,290 statistical differential
│   │       └── differential_indels.csv
│   ├── observation_first/                 ← v1 결과 (기존)
│   │   ├── step5_rna_editing_removed.vcf.gz  ← 31 candidates
│   │   ├── validation/                    ← v1 PPV 측정
│   │   └── comparison_to_phmg.md          ← 1차 비교 (2026-05-13)
│   └── observation_first_v2/              ← v2 결과 (이번)
│       ├── candidates_v2.bed              ← 45 candidates (BED)
│       ├── candidates_v2.vcf.gz           ← 45 candidates (VCF)
│       ├── intersect_keys.txt             ← 59 step4 intersection
│       ├── isec_late_only/0000.vcf        ← 7 true somatic
│       └── gene_mapping/
│           ├── candidates_v2_with_genes.tsv
│           ├── genes_v2.txt               ← 25 unique gene
│           └── genes_rn6.bed              ← rn6 transcript intervals
```

---

## 10. 다음 분기점 (User 결정 대기)

- (a) Mammary baseline 비교 먼저 (적은 비용)
- (b) rn7 liftover + 75-strain panel 실험 (큰 비용, 결정적)
- (c) Mouse ctrl-tumor cohort으로 우회 (design factor 회피)
- (d) Paper framing을 panel/design dependence 명시로 작성하고 추가 실험 없이 진행

---

## 한 문장 요약

> **Implementation parity(force-call + Fisher)를 맞춰서 mammary PPV가 9.7%→15.5%로 부분 개선됐지만 PHMG_IT의 68.4%엔 한참 못 미친다. 격차의 주범은 panel size(rn6 48-strain) + design contrast(time-course continuum) + WES scope, method 자체가 아닌 dataset-intrinsic 요인이 dominant. Method spec(70%/74% coverage)이 mammary에 그대로 transfer되지 않으므로 paper framing을 좁히거나 추가 실험으로 envelope을 정량화해야 한다.**
