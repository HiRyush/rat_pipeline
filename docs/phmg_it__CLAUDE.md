# PHMG_IT — Validation Dataset 1

> **현재 상태 (2026-05-06):** Phase 0-12 완료. 최적 pipeline 확정 (Observation-first + Imputation-as-filter). 657 true somatic, PPV 68.4%.
> **전체 프로젝트 현황:** 루트 `CLAUDE.md` 참조
> **이전 상세 이력:** `docs/archive/CLAUDE_20260429_PHMG_IT_phase0-8.md`

## 프로젝트 목표

RNA-seq 기반 variant calling 파이프라인 개발. Sample-matched DNA WGS를 ground truth로 사용하여 coverage와 정확도를 체계적으로 최적화한다.

## 데이터 정보

- **종:** Rat (Rattus norvegicus), Reference: UCSC RN7 (mRatBN7.2)
- **RNA-seq:** Paired-end FASTQ, 15 samples (C1-C5 Control, P1-P10 PHMG-treated)
- **DNA WGS (Ground Truth):** 15 matched gVCF (IsaacVariantCaller v2.0.13)

### 데이터 위치

| 데이터 | 경로 |
|--------|------|
| RNA FASTQ | `/media/yusanghyeon/30B4E366B4E32D52/01_Projects/korea_PHMG/RNA_fastq/{sample}_{1,2}.fastq.gz` |
| DNA gVCF | `/media/yusanghyeon/30B4E366B4E32D52/01_Projects/korea_PHMG/DNA_vcf/{dna_sample}_sorted.genome.vcf.gz` |
| Reference genome | `/home/yusanghyeon/RAT_project/PHMG_IT/reference/rn7.fa` (+.fai, .dict) |
| STAR index | `/home/yusanghyeon/RAT_project/PHMG_IT/reference/star_index/` |
| GTF annotation | `/home/yusanghyeon/RAT_project/PHMG_IT/reference/rn7.ncbiRefSeq.gtf` |
| 동료 코드 | `/home/yusanghyeon/RAT_project/PHMG_IT/coworker/` |
| 외장 HDD | `/media/yusanghyeon/30B4E366B4E32D52/01_Projects/korea_PHMG/` |

### 샘플 매핑 (DNA ↔ RNA)

| DNA (WGS) | RNA-seq | Group |
|-----------|---------|-------|
| P5S40W-2 | C1 | Control |
| P5S54W-2 | C2 | Control |
| P5S54W-3 | C3 | Control |
| P5S54W-6 | C4 | Control |
| P5S54W-8 | C5 | Control |
| P5H35W-11 | P1 | PHMG |
| P5H40W-2 | P2 | PHMG |
| P5H40W-16 | P3 | PHMG |
| P5H48W-12 | P4 | PHMG |
| P5H49W-9 | P5 | PHMG |
| P5H49W-18 | P6 | PHMG |
| P5H52W-8 | P7 | PHMG |
| P5H54W-2 | P8 | PHMG |
| P5H54W-19 | P9 | PHMG |
| P5M44W-11 | P10 | PHMG |

## 소프트웨어 환경

- **Conda 경로:** `/home/yusanghyeon/RAT_project/miniforge3/`
- **Conda 환경:** `rnaseq` (STAR 2.7.11b, GATK 4.6.1.0, bcftools 1.21, samtools 1.21, Beagle 5.5)
- **Conda 환경:** `snpeff` (SnpEff 5.3 + Java 21, annotation 전용)
- **시스템:** Ubuntu native, NVMe SSD 916GB
- conda가 PATH에 미등록 → 스크립트에서 직접 경로 사용

---

## Coverage 정의

본 프로젝트에서 "coverage"는 두 가지 의미:

### 1. Genome Coverage (전통적)
RNA-seq read가 genome의 몇 %를 덮는지. DP≥5 기준 **~41.5%** (RNA-seq 발현 영역 한계)

### 2. DNA SNP Coverage (핵심 지표)
**DNA WGS 전체 SNP 중 RNA-seq에서 탐지(복원) 가능한 비율.**

| 단계 | DNA SNP Coverage | 방법 |
|------|:---:|------|
| Phase 4 (bcftools) | ~6.5% | variant calling만 |
| Phase 9 (MuTect2 force-call) | ~8% | 전체 genome 대비 |
| **Phase 10 (Imputation)** | **78.6%** | HRDP panel + Beagle |

---

## 수행 이력 요약

> Phase 0~8의 상세 내용은 `CLAUDE_old_20260429.md` 참조

| Phase | 기간 | 핵심 결과 |
|-------|------|-----------|
| 0 | ~03-23 | HISAT2+bcftools 기반 초기 시도, Coverage ~21%로 한계 |
| 1 | 03-24 | STAR aggressive (Arm C) 발견: DP≥5 **45.6%** (+59%), SplitNCigar 제거 결정 |
| 2 | 03-26~30 | 15샘플 전체 Arm C 적용, 범용성 확인 |
| 3 | 03-30 | Multi-sample callable region 분석, K≥10에서도 개별 DP=0 문제 발견 |
| 4 | 03-31~04-01 | bcftools joint calling, Sens 6.5%, Prec 81.4% |
| 5 | 04-01 | 필터 최적화, DP≥5 한정 Sens 24%가 상한, 병목은 caller |
| 6 | 04-01~06 | FreeBayes 비교 → bcftools가 RNA-seq에서 우수 |
| 7 | 04-08 | Differential analysis: 48,470 SNPs + 7,974 INDELs (PHMG-specific) |
| 8 | 04-08 | SnpEff annotation: 189 HIGH impact, 449 MODERATE |

---

## Phase 9: MuTect2 + Force-Calling (2026-04-17 ~ 04-18, 완료)

**목적:** MuTect2 + force-calling으로 differential variant 수 증가

**핵심 수정사항 (GATK 4.6.1.0 + RNA-seq 호환):**
- `--disable-read-filter MappingQualityAvailableReadFilter --minimum-mapping-quality 0` (STAR MAPQ=255)
- `--min-median-base-quality 0` (RNA-seq base quality)
- `--annotations-to-exclude TandemRepeat` (GATK 4.6.1.0 버그 workaround)
- Force-call에서 VCF 생성 후 .stats 누락 시 initial scatter stats 복사

#### 결과

| 항목 | Phase 7 (bcftools) | Phase 9 (MuTect2) | 변화 |
|------|:---:|:---:|:---:|
| Differential SNPs | 48,470 | **65,375** | +34.8% |
| Differential INDELs | 7,974 | **19,498** | +144.5% |

**스크립트:** `scripts/mutect2_scatter.sh`
**결과:** `results/mutect2/differential/`

---

## Phase 10: Genotype Imputation (2026-04-22 ~ 04-24, 완료)

**목적:** RNA-seq anchor SNP + HRDP reference panel + Beagle imputation으로 DNA SNP Coverage 상향

**방법:**
1. HRDP 75 strain VCF merge (per-chr) + missing→0/0 + Beagle phasing
2. RNA-seq BAM에서 panel SNP 위치의 GL 추출 (`bcftools mpileup -T`)
3. Beagle 5.5 imputation
4. DNA WGS ground truth 대비 평가

#### 결과: DNA SNP Coverage (Before → After)

| Sample | DNA SNPs | Before% | **After%** |
|--------|:--------:|:-------:|:---:|
| C1 | 4,889K | 8.1 | **74.5** |
| C2 | 4,861K | 8.3 | **73.9** |
| C3 | 4,887K | 8.3 | **74.1** |
| C4 | 3,915K | 9.3 | **82.2** |
| C5 | 4,025K | 8.7 | **81.4** |
| P1 | 4,895K | 6.9 | **72.3** |
| P2 | 4,798K | 7.9 | **72.7** |
| P3 | 4,840K | 6.4 | **73.9** |
| P4 | 4,007K | 7.3 | **83.7** |
| P5 | 4,027K | 8.6 | **82.3** |
| P6 | 4,025K | 6.8 | **81.5** |
| P7 | 4,028K | 8.6 | **80.9** |
| P8 | 3,985K | 6.2 | **80.2** |
| P9 | 3,971K | 8.0 | **82.9** |
| P10 | 4,005K | 7.7 | **80.5** |
| **평균** | | **7.8** | **78.6** |

**Precision:** 74~82% (imputed ALT 중 DNA에도 있는 비율)

#### 결과: 유전자 캡처 수

| | Before (RNA-seq) | After (Imputed) | 신규 유전자 |
|--|:---:|:---:|:---:|
| **평균** | 18,775 (56.4%) | **24,165 (72.6%)** | **+7,588개** |
| 전체 유전자 | 33,287개 | | |

**스크립트:** `imputation/scripts/01_download_hrdp.sh`, `02_run_imputation_v2.sh`, `03_phase_and_impute.sh`
**결과:** `imputation/evaluation/imputation_results.tsv`, `imputation/evaluation/gene_capture_analysis.md`
**학술 근거:** BMC Genomics 2025, Beagle (Browning), HRDP (Cell Genomics 2024)
**가이드 문서:** `imputation/RNA-seq_Genotype_Imputation_Guide.md`

---

## Phase 11: Imputation 전후 Differential Analysis + Somatic Annotation (2026-04-29 ~ 05-04, 완료)

**목적:** Imputed VCF로 differential analysis를 재수행하여 imputation이 somatic detection에 미치는 영향 확인 + human somatic DB 교차 검증

#### Step 1: Imputed Differential Analysis

15개 imputed VCF merge → GT 기반 Control vs Treatment 비교

| Method | Differential SNPs |
|--------|:---:|
| Phase 7 (bcftools) | 48,470 |
| Phase 9 (MuTect2) | 65,375 |
| **Imputed (GT-based)** | **336,794** |

- Phase 7 ∩ Imputed = 977개 (2.0%) — imputation 편향을 뚫고 살아남은 variant
- Phase 7의 98%가 imputed에서 소실 (germline으로 덮임)
- Imputed only 335K는 개체 간 germline haplotype 차이 (noise)

#### Step 2: DNA Ground Truth 검증 (977개)

| Category | Count | % |
|----------|:---:|:---:|
| **Treatment DNA only (>=2 samples)** | **597** | **61.1%** |
| Treatment DNA only (1 sample) | 18 | 1.8% |
| Control DNA에도 존재 (germline) | 232 | 23.7% |
| DNA에 아예 없음 (RNA-only) | 130 | 13.3% |

#### Step 3: SnpEff + Human Somatic DB Annotation (597개)

**SnpEff:** HIGH 1, MODERATE 8, LOW 17, MODIFIER 571

**COSMIC CGC 매칭 (HIGH/MODERATE 중):**

| Rat Gene | Human | Identity | Effect | COSMIC |
|----------|:---:|:---:|---|:---:|
| Cdk9 | CDK9 | 98.7% | missense x2 | Yes |
| Mad1l1 | MAD1L1 | 82.9% | missense | Yes |
| Tango6 | HAS3 | 97.6% | missense | Yes |

**ClinVar somatic:** 8개 gene에서 somatic variant 보고됨 (CDH1: 19건+12 oncogenic 최다)
**Protein-level exact match:** 0개 — rat-specific somatic mutation

**결론 (방안 1):** Imputation은 differential/somatic analysis에 직접 적용 불가하나, Phase 7 결과와의 교집합(977개) 중 DNA 검증된 597개를 high-confidence somatic candidate로 확보. Gene-level에서 CDK9/MAD1L1/HAS3 등 human cancer gene 교차 검증 유의미.

**스크립트:** `scripts/imputed_differential_analysis.sh`
**결과:** `results/imputed_differential/`

---

## Phase 12: Somatic Capture — Imputed Baseline (2026-05-04, 완료)

**목적:** 방안 2 실행 — Imputed VCF를 per-sample germline baseline으로 사용, RNA-seq observed와의 차이에서 somatic capture

**방법:**
1. Per-sample: RNA-seq called VCF (DP>=5, QUAL>=5, alt_ratio>=0.15)
2. Imputed VCF와 비교: RNA ALT + Imputed REF(0|0) = somatic candidate
3. Control vs Treatment differential (recurrence>=2)

#### Per-sample 결과

| Group | Somatic candidates 평균 |
|-------|:---:|
| Control | 31,382 |
| Treatment | 26,452 |

#### Differential + DNA 검증

| 단계 | Count |
|------|:---:|
| Treatment-specific (rec>=2) | 15,861 |
| **DNA-validated true somatic** | **123** |
| Germline leak | 2,140 (13.5%) |
| RNA-only artifact | 13,490 (85.1%) |

#### 방안 1 + 방안 2 통합

| | 방안 1 | 방안 2 | 합계 |
|--|:---:|:---:|:---:|
| DNA-validated somatic | 597 | 123 | **720** |
| Overlap | — | — | **0** |

- 두 방안이 완전히 다른 somatic signal을 포착 (상호 보완적)
- 방안 1: imputation 편향을 뚫은 robust signal
- 방안 2: imputed germline baseline과의 차이 (noise 많으나 새로운 후보 발굴)

**스크립트:** `scripts/somatic_capture_imputed_baseline.sh`
**결과:** `results/somatic_capture/`

---

## 현재 진행 중

### mammary_cancer GL-based Imputation (2026-04-29 ~)
- **상태:** STAR + MarkDup + GL 추출 진행 중 (BATCH_SIZE=1, 순차 처리)
- **로그:** `~/mammary_gl_final.out`
- **dedup BAM 보존:** `mammary_cancer/results/dedup_bam/`
- **이전 hard call 결과:** 평균 Coverage 3.2% → 36.5% (GL 기반으로 추가 상향 예상)
- **HRDP rn6 panel:** 48 strain joint VCF (5.5GB), phased 완료

### mammary_cancer 완료된 작업
- RNA-seq variant calling (23샘플 HaplotypeCaller) — 완료
- Hard call 기반 imputation (Beagle) — 완료, 평균 36.5%
- Sample mapping SRA API로 검증 — 완료 (기존 mapping 전부 틀렸음, 수정 완료)

---

## 다음 단계

1. **mammary_cancer GL-based imputation 완료** → hard call(36.5%) vs GL 비교
2. **720개 통합 somatic candidate set 분석** — SnpEff/pathway 통합, 방안1+2 합산
3. **방안 2 noise 저감** — RNA-only 85% 문제 (필터 강화, RNA editing DB 제거 등)
4. **IH layer 적용** — 검증된 파이프라인으로 RNA-only 데이터 분석
5. DESeq2 differential expression 분석
6. Pathway analysis (somatic + germline + expression 통합)
7. COSMIC 라이선스 확보 시 position-level 매칭 재시도

---

## 주의사항

- 데이터(FASTQ, DNA VCF)는 외장 HDD에 위치 — 작업 전 마운트 확인
- DNA VCF는 gVCF 형식 (BLOCKAVG) — SNP 추출 시 `GT="alt"` 필터 필요
- GATK 4.6.1.0 TandemRepeat 버그 — `--annotations-to-exclude TandemRepeat` 필수
- STAR MAPQ=255 — GATK에서 `--disable-read-filter MappingQualityAvailableReadFilter` 필수
- 이전 상세 이력: `CLAUDE_old_20260429.md` 참조
