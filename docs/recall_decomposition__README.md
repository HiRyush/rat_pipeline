# Recall Decomposition Analysis (PHMG_IT)

**작성일:** 2026-05-17
**분석 단위:** PHMG_IT only (mammary는 generalizability 검증 단계에서 별도)
**상태:** 진행 중

---

## Motivation — 왜 이 분석이 필요한가

### 기존 coverage metric의 한계

기존 분석에서 사용한 "72.6% gene-level coverage"는 본질적으로 **germline-filter coverage** (HRDP imputation으로 germline 추정 가능한 gene 비율) 이지, **somatic detection recall** 이 아니다.

Root `CLAUDE.md` (framing-lock 조항)에서도 명시:
> "72.6% gene-level coverage" 단순 표현 금지 → "72.6% germline-filter coverage"로 정확히. Somatic-detection coverage와 혼동 금지.

### 새 metric 정의 — DNA-defined somatic recall

```
Recall = | RNA_captured_somatic ∩ DNA_treatment_specific_somatic | / | DNA_treatment_specific_somatic |
```

- 분자: 현재 pipeline output의 TRUE_SOMATIC (947개 in `validation_step4_intersection.tsv`)
- 분모: DNA WGS 기반 treatment-specific somatic set (Isaac caller VCF에서 ctrl=0 & treat≥k 정의)

### 초기 측정 결과 (naive denominator)

| Universe | treat≥2 | treat≥3 |
|---|---|---|
| Whole genome | 932/164,907 = **0.57%** | 591/65,073 = **0.91%** |
| MuTect2-evaluated | 932/14,555 = **6.40%** | 591/7,458 = **7.93%** |
| MuTect2 ∩ Transcript | 932/11,732 = **7.94%** | 591/6,018 = **9.82%** |

→ 한 자릿수 recall. 이는 method의 진짜 성능이 아니라, **분모에 "RNA가 원리상 볼 수 없는 위치"가 다수 포함된 결과**.

---

## 본 분석의 핵심 가설 — 분모 정제로 reachable recall 산정

분모(미캡처 DNA truth 위치, ~110K개)를 다음 4 bucket으로 분해:

| Bucket | 정의 | 의미 | 분모 처리 |
|---|---|---|---|
| **B1** | 모든 treatment RNA sample에서 mapped DP=0 | 발현 없음 (gene not expressed) | **분모에서 제거** (method's domain 밖) |
| **B2** | DP>0 but alt_count=0 in all treat samples | Reference로만 read 존재, alt-bearing read는 misalign/filter cut | **분모에서 제거** (정상 alignment된 alt 정보 자체가 없음) |
| **B3** | DP>0, alt_count>0 in ≥1 treat sample, but pipeline rejected | 진짜 method-level miss (필터가 컷한 것) | **분모에 포함** (개선 가능한 영역) |
| **B4** | Unmapped read pool에 alt-allele k-mer 존재 | Read 자체가 STAR alignment 실패 → unmapped로 빠짐 (reference bias의 가장 강한 케이스) | **분모에서 제거** (aligner-induced loss) |

### Reachable recall

```
Reachable recall = Captured / (Captured + B3)
                 = 947 / (947 + |B3|)
```

이는 "method가 원리상 detect 가능했던 위치 중 실제로 detect한 비율"이며, paper의 진짜 method recall claim이 된다.

### 부수 메트릭 — Aligner bias 정량

```
Aligner-induced loss rate = |B4| / |Total missed|
Expression-limited loss rate = |B1| / |Total missed|
```

이건 method 자체의 한계가 아니라 **upstream aligner (STAR)의 reference bias**가 가져오는 손실로, method paper에서 별도 figure로 보고 가능.

---

## 데이터 입력

| 파일 | 위치 | 용도 |
|---|---|---|
| DNA truth counts | `../dna_truth_coverage/dna_truth_counts.tsv.gz` | 분모 set 구축 |
| DNA truth (Option A) | `../dna_truth_coverage/dna_truth_counts_OptA.tsv.gz` | 분모 set (treat=PASS, ctrl=any) |
| RNA captured | `~/RAT_project/pipeline_candidates/results/observation_first/validation_step4_intersection.tsv` | 분자 |
| RNA markdup BAMs | `../mutect2/markdup/{C1-C5,P1-P10}.dedup.bam` | mpileup + unmapped 추출 |
| MuTect2 union sites | `../mutect2/union_sites.vcf.gz` | universe restriction |
| Reference genome | `../../reference/rn7.fa` | k-mer 합성 |
| Gene BED | `../dna_truth_coverage/genes.bed` | gene-level recall 보조 metric |

---

## 분석 단계

### Step 1 — Build inputs
`scripts/01_build_inputs.sh`
- DNA truth BED (treat≥2 & ctrl=0, Option A)
- RNA captured BED
- Missed set = DNA truth − RNA captured

### Step 2 — B1/B2/B3 classification
`scripts/02_classify_b1b2b3.sh`
- 10개 treatment BAM에 대해 missed 위치 `samtools mpileup -l missed.bed`
- 위치별 sum: total_dp, total_alt_dp
- 분류:
  - B1: total_dp = 0
  - B2: total_dp > 0 & total_alt_dp = 0
  - B3: total_dp > 0 & total_alt_dp > 0

### Step 3 — B4 detection via unmapped k-mer
`scripts/03_b4_kmer_scan.sh`
- Non-B3 missed 위치마다 reference에서 alt-flanking 31-mer 합성 (forward + revcomp)
- 10개 treatment BAM에서 `samtools view -f 4` 로 unmapped read 추출
- K-mer scan (KMC 또는 grep 기반)
- K-mer present in unmapped pool ≥ threshold → B4

### Step 4 — Reachable recall 산정
`scripts/04_compute_recall.py`
- Bucket sizes 집계
- Reachable recall 계산 (여러 universe와 treat threshold 조합으로)
- 최종 표 + JSON 출력

---

## 결과 파일 (예정)

| 파일 | 내용 |
|---|---|
| `data/missed_truth.bed` | DNA truth minus captured |
| `data/bucket_b1b2b3.tsv` | per-position B1/B2/B3 분류 + DP/alt_DP 값 |
| `data/bucket_b4.tsv` | B4 후보 k-mer match 정보 |
| `results/bucket_summary.tsv` | bucket count 요약 |
| `results/reachable_recall.tsv` | 최종 recall 매트릭스 |
| `results/recall_decomposition_report.md` | 분석 결과 narrative |

---

## 대화 요약 — 분석 setup에 이르기까지

1. **요청 출발점:** "IT layer에서 coverage를 다르게 계산하자. 우리가 capture한 RNA somatic이 DNA somatic을 얼마나 cover하는지가 핵심."

2. **분모 정의 합의:** DNA에서 ctrl=0 & treat≥k인 set이 분모로 가야 한다는 데 동의. naive 정의로 측정 → recall 1-10% (universe에 따라).

3. **분모 strictness 논의:** Option A (treat=PASS, ctrl=any-filter)가 가장 정직. Option C (5 control 모두 DP≥5 confirmed REF)는 gVCF가 필요한데 현재 local에 없음 (외장 HDD, 사용자 원격 작업 중).

4. **Reference bias 가설 (사용자 제기):**
   > "RNA에서 align이 안되서 버려지는 것들 중 mutation이 있는 애들이 있을 수도 있는거 아님?"

   → 이게 분석의 핵심 motivation. Naive recall이 낮은 이유 중 큰 비중이 method 한계가 아니라 aligner bias 일 가능성.

5. **DNA 재처리 검토:** 현재 SNP VCF만 local에 있고 gVCF/BAM/FASTQ는 외장 HDD. DNA-side 개선은 외장 HDD 마운트 필요 — 사용자가 연구실 복귀 후로 미루기로 결정.

6. **현 자원으로 가능한 가장 강한 lever:** RNA BAM 기반 B1/B2/B3 분해 + unmapped 풀 k-mer scan으로 B4 측정. 이 4-bucket 분해가 "reachable recall" 정의를 가능케 함.

7. **사용자 결정:** 새 폴더에서 분석 1-4 단계를 모두 진행. 본 README가 그 기록.

---

## 주의사항

- **Framing 유지:** Reachable recall은 "method recall under aligner constraint" 임. "DNA somatic 전체에 대한 sensitivity"가 아니라는 점 paper에 명시 필수.
- **Bucket 정의의 정직성:** B2를 분모에서 빼는 것에 대한 reviewer 공격 가능 — "alt allele이 충분 evidence 있는데도 caller가 못 잡으면 그건 method miss 아니냐" 반론. 따라서 B2 정의에 alt depth = 0 strict하게 적용 (1 read만 있어도 B3로 분류).
- **B4 false positive 우려:** Unmapped read에 k-mer가 우연히 매칭될 가능성. K-mer 길이를 31bp 이상, hit count ≥ 2 등 컷 적용.
- **결과는 PHMG_IT only.** Cross-dataset generalizability는 mammary에 같은 framework 적용 후 별도 검증.
