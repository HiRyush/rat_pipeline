# PHMG_IH: PHMG Inhalation Model — RNA-seq DNA Mutation Prediction

## 프로젝트 위치 및 목적

**Pipeline 적용 단계** (3단계 연구 설계 중 마지막)

```
PHMG_IT (DNA+RNA, Intratracheal) → 파이프라인 구축 & 최적화  ✅ 완료
mammary_cancer (DNA+RNA)         → 2차 검증 (범용성 확인)     ⬜ 진행 중
PHMG_IH (RNA only, Inhalation)   → DNA mutation prediction    ⬜ 이 프로젝트
```

PHMG_IT에서 최적화하고 mammary_cancer에서 검증한 RNA-seq variant calling 파이프라인을
**DNA ground truth 없이** RNA만 존재하는 Inhalation 데이터셋에 적용하여 DNA mutation을 예측한다.

## 데이터 정보

- **종:** Rat (Rattus norvegicus)
- **레퍼런스:** rn7 (mRatBN7.2) — PHMG_IT와 동일
- **투여 방법:** Inhalation (흡입)
- **RNA-seq:** Paired-end FASTQ (~145GB)
- **DNA WGS:** 없음 (RNA only dataset) — ground truth validation 불가, 의도된 설계
- **분석 목표:** Condition-specific DNA mutation prediction (PHMG-treated vs Control)

## 데이터 위치

| 데이터 | 경로 |
|--------|------|
| RNA FASTQ | `/home/yusanghyeon/RAT_project/PHMG_IH/fastq/` |
| Reference genome (rn7) | `/home/yusanghyeon/RAT_project/PHMG_IT/reference/rn7.fa` (공유) |
| STAR index (rn7) | `/home/yusanghyeon/RAT_project/PHMG_IT/reference/star_index/` (공유) |
| GTF annotation | `/home/yusanghyeon/RAT_project/PHMG_IT/reference/rn7.ncbiRefSeq.gtf` (공유) |

## 소프트웨어 환경

- PHMG_IT와 동일: conda env `rnaseq` (`/home/yusanghyeon/RAT_project/miniforge3/envs/rnaseq/`)
- STAR 2.7.11b, GATK 4.6.1.0, bcftools 1.21, samtools 1.21

## 파이프라인 (PHMG_IT에서 검증된 파라미터 적용)

1. **STAR alignment** — Arm C (aggressive, SplitNCigar 생략)
2. **MarkDuplicates**
3. **BQSR** (known_sites 필요 — PHMG_IT의 것 재사용 가능)
4. **HaplotypeCaller** or **bcftools mpileup+call**
5. **VariantFiltration / Filter optimization** — PHMG_IT Phase 5 결과 적용
6. **Differential variant analysis** — PHMG_IH 내부 Control vs Treated

## 주의사항

- DNA ground truth 없음 → Sensitivity/Precision 수치 계산 불가
- 검증은 PHMG_IT + mammary_cancer에서 이미 완료된 상태로 적용
- Reference는 PHMG_IT와 공유 (별도 복사 불필요)
