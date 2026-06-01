# Sample Mapping Evidence: RNA-seq (SRR) ↔ WES (GSM)

## 검증 일시
2026-04-24

## 데이터 출처
- **RNA-seq**: GEO Series GSE297548 (SRA: SRP586338)
- **WES**: GEO Series GSE297544
- **논문**: Yan P et al., "Aging-associated differences in mammary tumor-initiating populations and immune evasion pathways in breast cancer" (Dana-Farber Cancer Institute)

## 매핑 방법

### 1단계: SRR → RNA-seq GSM 확인 (NCBI SRA API)
각 SRR accession에 대해 `https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc=SRR...`에서:
- `EXTERNAL_ID namespace="GEO"` → RNA-seq GSM accession
- `<TITLE>` → tumor name (e.g., "M20_4")
- `filename=` → 원본 FASTQ 파일명 (e.g., "20M-4_S1_L001_R1_001.fastq.gz")

### 2단계: WES GSM → tumor name 확인 (GEO)
GSE297544에서 각 WES sample의 GSM ID와 VCF 파일명으로 tumor name 확인:
- VCF 파일명 패턴: `GSM{id}_{age}m_{tumor_no}_tumor.pass.vcf`
- 예: GSM8994506_20m_4_tumor.pass.vcf → M20 tumor 4

### 3단계: tumor name으로 매칭
RNA-seq와 WES 모두 같은 tumor에서 유래:
- RNA-seq "M20_4" = WES "20m_4" = 같은 종양

### 검증 근거 (3중 확인)
| 증거 | 출처 | 예시 (SRR33625148) |
|------|------|-------------------|
| RNA GSM | SRA metadata | GSM8994567 |
| Tumor name | GEO sample title | M20_4 |
| 원본 FASTQ 파일명 | SRA original file | `20M-4_S1_L001_R1_001.fastq.gz` |
| WES VCF 파일 존재 | 로컬 파일 확인 | GSM8994506_20m_4_tumor.pass.vcf ✓ |

**3가지 독립적 증거가 모두 일치함.**

## 이전 매핑 오류
기존 `sample_mapping.tsv`는 GEO 등록 순서로 추정했으나, SRR 번호가 GSM 역순이었기 때문에 전부 틀렸음:
- 잘못된 매핑: SRR33625148 → M20_1
- 올바른 매핑: SRR33625148 → **M20_4**

## 주의사항
- M1_4와 M6_1은 RNA-seq 라이브러리에서 제외됨 (사후 시간 경과로 RNA 품질 부적합)
  - 논문: "Two tumors (1m-4 and 6m-1) were excluded from library preparation due to extended postmortem intervals"
- M1_8 (GSM8994489): WES VCF 존재, RNA-seq SRR33625164로 매칭됨

## 최종 매핑 테이블

| SRR | RNA GSM | Tumor | WES GSM | WES VCF | Verified |
|-----|---------|-------|---------|---------|----------|
| SRR33625148 | GSM8994567 | M20_4 | GSM8994506 | GSM8994506_20m_4_tumor.pass.vcf | ✓ |
| SRR33625149 | GSM8994566 | M20_3 | GSM8994505 | GSM8994505_20m_3_tumor.pass.vcf | ✓ |
| SRR33625150 | GSM8994565 | M20_2 | GSM8994504 | GSM8994504_20m_2_tumor.pass.vcf | ✓ |
| SRR33625151 | GSM8994564 | M20_1 | GSM8994503 | GSM8994503_20m_1_tumor.pass.vcf | ✓ |
| SRR33625152 | GSM8994563 | M12_3 | GSM8994502 | GSM8994502_12m_3_tumor.pass.vcf | ✓ |
| SRR33625153 | GSM8994562 | M12_2 | GSM8994501 | GSM8994501_12m_2_tumor.pass.vcf | ✓ |
| SRR33625154 | GSM8994561 | M12_1 | GSM8994500 | GSM8994500_12m_1_tumor.pass.vcf | ✓ |
| SRR33625155 | GSM8994560 | M6_2 | GSM8994499 | GSM8994499_6m_2_tumor.pass.vcf | ✓ |
| SRR33625156 | GSM8994559 | M3_8 | GSM8994497 | GSM8994497_3m_8_tumor.pass.vcf | ✓ |
| SRR33625157 | GSM8994558 | M3_7 | GSM8994496 | GSM8994496_3m_7_tumor.pass.vcf | ✓ |
| SRR33625158 | GSM8994557 | M3_6 | GSM8994495 | GSM8994495_3m_6_tumor.pass.vcf | ✓ |
| SRR33625159 | GSM8994556 | M3_5 | GSM8994494 | GSM8994494_3m_5_tumor.pass.vcf | ✓ |
| SRR33625160 | GSM8994555 | M3_4 | GSM8994493 | GSM8994493_3m_4_tumor.pass.vcf | ✓ |
| SRR33625161 | GSM8994554 | M3_3 | GSM8994492 | GSM8994492_3m_3_tumor.pass.vcf | ✓ |
| SRR33625162 | GSM8994553 | M3_2 | GSM8994491 | GSM8994491_3m_2_tumor.pass.vcf | ✓ |
| SRR33625163 | GSM8994552 | M3_1 | GSM8994490 | GSM8994490_3m_1_tumor.pass.vcf | ✓ |
| SRR33625164 | GSM8994551 | M1_8 | GSM8994489 | GSM8994489_1m_8_tumor.pass.vcf | ✓ |
| SRR33625165 | GSM8994550 | M1_7 | GSM8994488 | GSM8994488_1m_7_tumor.pass.vcf | ✓ |
| SRR33625166 | GSM8994549 | M1_6 | GSM8994487 | GSM8994487_1m_6_tumor.pass.vcf | ✓ |
| SRR33625167 | GSM8994548 | M1_5 | GSM8994486 | GSM8994486_1m_5_tumor.pass.vcf | ✓ |
| SRR33625168 | GSM8994547 | M1_3 | GSM8994484 | GSM8994484_1m_3_tumor.pass.vcf | ✓ |
| SRR33625169 | GSM8994546 | M1_2 | GSM8994483 | GSM8994483_1m_2_tumor.pass.vcf | ✓ |
| SRR33625170 | GSM8994545 | M1_1 | GSM8994482 | GSM8994482_1m_1_tumor.pass.vcf | ✓ |

**23/23 매핑 검증 완료. 모든 WES VCF 파일 존재 확인.**
