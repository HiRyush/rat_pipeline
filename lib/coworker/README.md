# Mouse Variant-Expression Integration Pipeline

RNA-seq 데이터의 VCF 파일을 이용하여 **Control vs Disease** 또는 **Normal vs Tumor** 조건 간
variant 차이와 유전자 발현 차이를 통합 분석하는 Python 파이프라인입니다.

Zebrafish RNAmapper (Miller et al., 2013)의 VCF 파싱 개념을 기반으로,
멘델리안 linkage mapping 대신 **조건 간 differential variant + DESeq2 통합** 전략으로 재설계하였습니다.

---

## 핵심 특징

| 기능 | 설명 |
|------|------|
| **Dual VCF 지원** | bcftools mpileup/call 과 GATK HaplotypeCaller 포맷 자동 감지 |
| **SNP + INDEL 분리 분석** | INDEL을 제거하지 않고 별도 트랙으로 분석 (frameshift 등 고영향 variant 보존) |
| **두 가지 분석 모드** | `independent` (disease vs WT) / `paired` (tumor vs normal) |
| **Variant-Expression 통합** | Differential variant + DESeq2 결과를 Tier별로 통합하여 candidate gene list 도출 |
| **Pure Python** | 외부 R 의존성 없음 (scipy만 optional) |

---

## 설치 및 의존성

```bash
# 기본 (필수)
Python >= 3.8

# 선택 (Fisher's exact test 사용 시)
pip install scipy
```

---

## 파이프라인 구조

```
RNA-seq BAM files
       │
       ├──── bcftools mpileup/call ────┐
       │                                ├──► VCF files (input)
       └──── GATK HaplotypeCaller ─────┘
                                          │
                    ┌─────────────────────┘
                    │
    ┌───────────────▼───────────────┐
    │  Step 1: VCF Parsing          │  vcf_parser.py
    │  - 자동 포맷 감지             │  bcftools DP4 / GATK AD 통합 파싱
    │  - SNP / INDEL 분류           │  frameshift 자동 판별
    │  - near_indel flagging        │  INDEL 근처 SNP에 flag (제거 안함)
    └───────────────┬───────────────┘
                    │
    ┌───────────────▼───────────────┐
    │  Step 2: Differential Variant │  variant_analyzer.py
    │  - Mode A: Independent        │  recurrence + reference 제거
    │  - Mode B: Paired             │  개체 내 somatic variant 검출
    │  - Fisher's exact test        │  통계적 유의성 검정
    │  → SNP/INDEL 별도 CSV 출력    │
    │  → VEP 입력용 VCF 출력        │
    └───────────────┬───────────────┘
                    │
           [VEP / SnpEff 실행]  ← 외부 도구
                    │
    ┌───────────────▼───────────────┐
    │  Step 3: Integration          │  integration.py
    │  - VEP annotation 파싱        │
    │  - DESeq2 결과 로드           │
    │  - Tier별 candidate 분류      │
    │  → candidate_genes.csv        │
    └───────────────────────────────┘
```

---

## 사용법

### 기본 실행 (Step 1-2)

```bash
# Disease model vs WT (independent)
python run_pipeline.py --mode independent \
    --target vcf/disease_1.vcf vcf/disease_2.vcf vcf/disease_3.vcf \
    --reference vcf/ctrl_1.vcf vcf/ctrl_2.vcf vcf/ctrl_3.vcf \
    --output results/

# Tumor vs Normal (paired) — target과 reference 순서가 pair를 결정
python run_pipeline.py --mode paired \
    --target vcf/tumor_1.vcf vcf/tumor_2.vcf vcf/tumor_3.vcf \
    --reference vcf/normal_1.vcf vcf/normal_2.vcf vcf/normal_3.vcf \
    --output results/
```

### VEP 실행 (외부)

```bash
vep -i results/for_vep.vcf \
    -o results/vep_output.txt \
    --species mus_musculus \
    --cache --sift b --polyphen b --symbol --canonical
```

### 전체 파이프라인 (Step 1-3)

```bash
python run_pipeline.py --mode independent \
    --target vcf/disease_1.vcf vcf/disease_2.vcf vcf/disease_3.vcf \
    --reference vcf/ctrl_1.vcf vcf/ctrl_2.vcf vcf/ctrl_3.vcf \
    --de-genes deseq2_results.csv \
    --vep-results results/vep_output.txt \
    --output results/
```

### Config 파일 사용

```bash
python run_pipeline.py --config config.yaml
```

---

## 파라미터

| 파라미터 | 기본값 | 설명 |
|---------|--------|------|
| `--min-coverage` | 10 | 최소 read depth |
| `--min-qual` | 30 | 최소 QUAL score |
| `--min-alt-ratio` | 0.3 | Target에서의 최소 alt allele 비율 |
| `--max-ref-alt-ratio` | 0.1 | Reference에서의 최대 alt allele 비율 (이 이상이면 shared로 판정) |
| `--min-recurrence` | 2 | 최소 몇 개 target sample에서 발견되어야 하는지 |
| `--padj-cutoff` | 0.05 | DE gene 판정 adjusted p-value |
| `--log2fc-cutoff` | 1.0 | DE gene 판정 \|log2FC\| |

### 시나리오별 권장 파라미터

```
Disease model vs WT:
  --min-alt-ratio 0.3    (heterozygous variant 포함)
  --max-ref-alt-ratio 0.1
  --min-recurrence 2     (3개 중 2개 이상)

Tumor vs Normal:
  --min-alt-ratio 0.1    (subclonal mutation 허용)
  --max-ref-alt-ratio 0.05
  --min-recurrence 1     (paired이므로 1부터 유의미)
```

---

## 출력 파일

```
results/
├── differential_snps.csv        # Condition-specific SNPs
├── differential_indels.csv      # Condition-specific INDELs
├── differential_all_variants.csv
├── for_vep.vcf                  # VEP 입력용 VCF
├── candidate_genes.csv          # 전체 candidate gene list
├── candidates_tier1.csv         # High-impact variant + DE
├── candidates_tier2.csv         # Any variant + DE
├── candidates_tier3.csv         # High-impact variant, no DE
├── candidates_tier4.csv         # Low-impact variant only
└── candidates_tier5.csv         # DE only (variant 미발견)
```

### Candidate Gene Tier 체계

```
Tier 1: High-impact variant + DE gene          ← 최우선 후보
        (stop_gained, frameshift, splice 등 + padj<0.05 + |log2FC|>1)

Tier 2: Any variant + DE gene
        (missense, inframe 등 + DE)

Tier 3: High-impact variant, no DE
        (gain-of-function 가능성, protein-level 검증 필요)

Tier 4: Low-impact variant only
        (추가 기능 분석 필요)

Tier 5: DE only
        (variant가 regulatory region/intron에 있거나 RNA-seq coverage 부족 가능)
```

---

## bcftools vs GATK VCF 포맷 차이

이 파이프라인이 자동 처리하는 두 포맷의 핵심 차이:

```
┌──────────────┬─────────────────────────┬─────────────────────────┐
│              │ bcftools                │ GATK HaplotypeCaller    │
├──────────────┼─────────────────────────┼─────────────────────────┤
│ Allele depth │ INFO: DP4=F,R,F,R      │ FORMAT: AD=ref,alt      │
│ Total depth  │ INFO: DP               │ FORMAT: DP              │
│ Genotype     │ 없거나 GT in FORMAT    │ GT in FORMAT            │
│ Strand info  │ 있음 (Forward/Reverse) │ 없음                    │
│ INDEL 표시   │ INFO에 "INDEL" flag    │ REF/ALT 길이 차이       │
│ Quality      │ QUAL column            │ QUAL + FORMAT:GQ        │
└──────────────┴─────────────────────────┴─────────────────────────┘
```

---

## INDEL 처리 방식 (RNAmapper와의 차이)

| 항목 | 기존 RNAmapper | 이 파이프라인 |
|------|---------------|-------------|
| INDEL 자체 | 제거 | **보존, 별도 트랙 분석** |
| INDEL 근처 10bp SNP | 제거 | **flag만 추가 (near_indel=True), 유지** |
| VEP 입력 | SNP만 | **SNP + INDEL 모두** |
| 이유 | BSA marker 정확도 | Frameshift 등 고영향 variant가 candidate |

---

## 테스트 실행

```bash
# 테스트 데이터 생성 (bcftools + GATK 혼합)
python generate_test_data.py

# 파이프라인 실행
python run_pipeline.py --mode independent \
    --target test_data/vcf/disease_1.vcf test_data/vcf/disease_2.vcf test_data/vcf/disease_3.vcf \
    --reference test_data/vcf/ctrl_1.vcf test_data/vcf/ctrl_2.vcf test_data/vcf/ctrl_3.vcf \
    --de-genes test_data/deseq2_results.csv \
    --vep-results test_data/vep_output.txt \
    --output test_data/results/

# 기대 결과: 6개 disease-specific variant 검출 (SNP 3 + INDEL 3)
#            Tier 1 candidate: Apc, Brca1, Pten, Kras
```

---

## 향후 확장

- **Pathway enrichment**: candidate gene list → clusterProfiler / g:Profiler
- **Allele-Specific Expression (ASE)**: 같은 유전자에서 두 allele의 발현 비율 비교
- **Splicing analysis**: rMATS / SUPPA2로 alternative splicing 차이 분석
- **Mutation signature**: tumor somatic variant의 mutational signature 분석
- **Multi-sample genotyping**: GATK GenotypeGVCFs로 joint calling 후 분석
