# rat_pipeline — RNA-only somatic-like variant detection (정리된 코드 저장소)

RNA-seq만으로 (matched DNA 없이) somatic-like / treatment-enriched non-germline variant를
검출하는 **method/pipeline**의 코드를 raw fastq → 최종 결과 순서대로 step별 정리한 곳.

> **현재 scope: PHMG_IT(Rat rn7) 전용.** mammary_cancer cross-validation은 아직 진행 중이라
> 이 repo에는 포함하지 않는다. mammary가 마무리되면 같은 step 레이아웃에 `mammary__*` prefix로 추가한다.

- **Deliverable은 method 자체** — PHMG_IT는 method를 개발·검증한 벤치마크(Validation 1), PHMG_IH는 RNA-only 최종 적용 대상.
- 데이터(BAM/VCF/FASTQ, ~540GB)는 추적하지 않는다. 스크립트는 내부에 **절대경로**를 박고 있어
  어느 폴더에서 실행해도 `/home/yusanghyeon/RAT_project/...`의 원본 데이터를 그대로 읽는다.
  경로 한눈에 보기: [`config/paths.env`](config/paths.env).
- 이 폴더는 원본 코드의 **비파괴 복사본**이다. 신규경로↔원본경로 대응은 [`MANIFEST.tsv`](MANIFEST.tsv).
  복사를 재현하는 스크립트는 [`build_repo.sh`](build_repo.sh).

## 확정 파이프라인 (Observation-first + Imputation-as-Filter)

```
raw fastq
  └─ 01_alignment        STAR (MAPQ=255 보정)                 → BAM
  └─ 02_variant_calling  MuTect2 force-call(공통 좌표) + Filter → PASS VCF
  └─ 03_differential     Fisher exact + recurrence + AF (treatment vs control)
  └─ 04_imputation       BEAGLE + HRDP panel → germline 추정
  └─ 05_intersect_filter (Step3 ∩ Imputed differential) + RNA editing 제거 → 최종 candidate
  └─ 06_validation       DNA truth(WGS/WES)로 PPV / recall 채점
```

> 핵심 교훈: **Imputation은 detector가 아니라 filter로 써야 한다** (detector=0개, filter=PPV 67%).
> Ablation(`analysis/ablation/reproduce_ablation.py`, 재현 가능)에서 imputation filter step이
> PPV를 ~2% → 68.9%로 끌어올린 것이 main claim. 최종 candidate=961, TRUE_SOMATIC=657, PPV 68.4%.

## 디렉토리

| 경로 | 내용 |
|---|---|
| `pipeline/00_reference` | HRDP panel 다운로드, genome/index 준비 |
| `pipeline/01_alignment` | STAR alignment (fastq → BAM) |
| `pipeline/02_variant_calling` | MuTect2 force-call + FilterMutectCalls |
| `pipeline/03_differential` | Fisher differential (`run_pipeline.py` = 핵심) |
| `pipeline/04_imputation` | BEAGLE + HRDP imputation |
| `pipeline/05_intersect_filter` | 교집합 + RNA editing/artifact filter (`observation_first_v2.py` = 최종) |
| `pipeline/06_validation` | DNA truth PPV/recall 채점 |
| `analysis/ablation` | component ablation |
| `analysis/coverage_gene_mapping` | 최종 candidate → gene mapping, coverage (numerator 961/1,356로 재계산 필요) |
| `analysis/ablation` | `reproduce_ablation.py` — 재현 가능한 ablation 캐스케이드 (A→B→C→D) |
| `analysis/recall_decomposition` | recall 분해 분석 (B1~B4, salvage) |
| `experimental/` | superseded 탐색 arm (freebayes, joint_calling, arms B/C/D, discordance detector) |
| `lib/coworker` | Fisher differential / integration Python 패키지 (modules/ 통째, step02·03이 사용) |
| `docs/` | method 설계·현황 문서 |

> 파일명 규칙: step 폴더 안에서 `phmg_it__*` prefix로 dataset 구분 (현재 PHMG_IT 전용).
> step이 1차축이라 mammary 추가 시 같은 단계 구현을 나란히 비교할 수 있다 (cross-dataset 일관성이 method paper의 향후 main claim).

## 데이터셋 (이 repo: PHMG_IT만)

| Dataset | 종/ref | 샘플 | DNA truth | 역할 | repo 포함 |
|---|---|---|---|---|---|
| PHMG_IT | Rat / rn7 | 15 (5C+10T) | WGS | 개발 + Validation 1 (PPV 68.4%) | ✅ |
| mammary_cancer | Rat / rn6 | 23 | WES | Validation 2 (cross-validation 진행 중) | ⏳ 미포함 |
| PHMG_IH | Rat / rn7 | TBD | 없음 | 최종 적용 (RNA-only) | — |

## 실행 메모

```bash
source /home/yusanghyeon/RAT_project/miniforge3/etc/profile.d/conda.sh
conda activate rnaseq            # STAR, GATK 4.6.1.0, bcftools, Beagle 5.5
# 각 step 스크립트를 순서대로 실행. 데이터 입출력 경로는 스크립트 내부에 하드코딩되어 있음.
```
