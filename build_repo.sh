#!/usr/bin/env bash
# build_repo.sh — RAT_project 분산 코드를 step별로 비파괴 복사 정리.
# 원본은 절대 건드리지 않음 (cp only). 재실행 안전 (idempotent).
set -euo pipefail

SRC="/home/yusanghyeon/RAT_project"
DST="$SRC/rat_pipeline"
MAN="$DST/MANIFEST.tsv"

mkdir -p "$DST"/{config,lib,docs}
mkdir -p "$DST"/pipeline/{00_reference,01_alignment,02_variant_calling,03_differential,04_imputation,05_intersect_filter,06_validation}
mkdir -p "$DST"/analysis/{ablation,coverage_gene_mapping,recall_decomposition,cross_validation}
mkdir -p "$DST"/experimental/{freebayes,joint_calling,arms_BCD,discordance}

printf 'new_path\torigin_path\n' > "$MAN"

# cp_one <src-relative-to-SRC> <dst-relative-to-DST>
cp_one() {
  local s="$SRC/$1" d="$DST/$2"
  if [[ ! -e "$s" ]]; then echo "  !! MISSING: $1" >&2; return 0; fi
  mkdir -p "$(dirname "$d")"
  cp -p "$s" "$d"
  printf '%s\t%s\n' "$2" "$1" >> "$MAN"
}

# ---- lib: coworker Python 패키지 통째 ----
cp -rp "$SRC/PHMG_IT/coworker" "$DST/lib/coworker"
printf 'lib/coworker/\tPHMG_IT/coworker/ (whole package, modules intact)\n' >> "$MAN"

# ---- 00_reference ----
cp_one PHMG_IT/imputation/scripts/01_download_hrdp.sh        pipeline/00_reference/phmg_it__01_download_hrdp.sh

# ---- 01_alignment ----
cp_one PHMG_IT/coworker/star_alignment_unsorted.sh          pipeline/01_alignment/phmg_it__star_alignment_unsorted.sh
cp_one PHMG_IT/scripts/star_alignment_aggressive.sh         pipeline/01_alignment/phmg_it__star_alignment_aggressive.sh
cp_one mammary_cancer/scripts/run_rnaseq_variant_calling.sh pipeline/01_alignment/mammary__run_rnaseq_variant_calling.sh
cp_one mammary_cancer/scripts/run_all_samples.sh            pipeline/01_alignment/mammary__run_all_samples.sh

# ---- 02_variant_calling ----
cp_one PHMG_IT/scripts/mutect2_scatter.sh                   pipeline/02_variant_calling/phmg_it__mutect2_scatter.sh
cp_one PHMG_IT/scripts/mutect2_batch.sh                     pipeline/02_variant_calling/phmg_it__mutect2_batch.sh
cp_one PHMG_IT/scripts/mutect2_per_sample.sh                pipeline/02_variant_calling/phmg_it__mutect2_per_sample.sh
cp_one PHMG_IT/scripts/run_remaining_9_samples.sh           pipeline/02_variant_calling/phmg_it__run_remaining_9_samples.sh
cp_one PHMG_IT/coworker/gatk_rnaseq_variant_calling_broad.sh pipeline/02_variant_calling/phmg_it__gatk_rnaseq_variant_calling_broad.sh
cp_one mammary_cancer/scripts/02_mutect2_per_sample.sh      pipeline/02_variant_calling/mammary__02_mutect2_per_sample.sh
cp_one mammary_cancer/scripts/03_run_parallel.sh            pipeline/02_variant_calling/mammary__03_run_parallel.sh
cp_one mammary_cancer/scripts/04_filter_mutect.sh           pipeline/02_variant_calling/mammary__04_filter_mutect.sh
cp_one mammary_cancer/scripts/11_force_call_scatter.sh      pipeline/02_variant_calling/mammary__11_force_call_scatter.sh
cp_one mammary_cancer/scripts/12_filter_force_call.sh       pipeline/02_variant_calling/mammary__12_filter_force_call.sh

# ---- 03_differential ----
cp_one PHMG_IT/scripts/differential_variant_analysis.sh     pipeline/03_differential/phmg_it__differential_variant_analysis.sh
cp_one PHMG_IT/coworker/run_pipeline.py                     pipeline/03_differential/phmg_it__run_pipeline.py
cp_one pipeline_candidates/scripts/step4_differential_filter.py pipeline/03_differential/phmg_it__step4_differential_filter.py
cp_one mammary_cancer/scripts/05_differential.sh            pipeline/03_differential/mammary__05_differential.sh

# ---- 04_imputation ----
cp_one PHMG_IT/imputation/scripts/02_run_imputation.sh      pipeline/04_imputation/phmg_it__02_run_imputation.sh
cp_one PHMG_IT/imputation/scripts/02_run_imputation_v2.sh   pipeline/04_imputation/phmg_it__02_run_imputation_v2.sh
cp_one PHMG_IT/imputation/scripts/03_phase_and_impute.sh    pipeline/04_imputation/phmg_it__03_phase_and_impute.sh
cp_one mammary_cancer/scripts/run_imputation.sh             pipeline/04_imputation/mammary__run_imputation.sh
cp_one mammary_cancer/scripts/run_imputation_with_gl.sh     pipeline/04_imputation/mammary__run_imputation_with_gl.sh

# ---- 05_intersect_filter ----
cp_one PHMG_IT/scripts/imputed_differential_analysis.sh     pipeline/05_intersect_filter/phmg_it__imputed_differential_analysis.sh
cp_one PHMG_IT/scripts/somatic_capture_imputed_baseline.sh  pipeline/05_intersect_filter/phmg_it__somatic_capture_imputed_baseline.sh
cp_one pipeline_candidates/scripts/observation_first_pipeline.py pipeline/05_intersect_filter/phmg_it__observation_first_pipeline.py
cp_one pipeline_candidates/scripts/observation_first_v2.py  pipeline/05_intersect_filter/phmg_it__observation_first_v2.py
cp_one pipeline_candidates/scripts/step5_deseq2_filter.R    pipeline/05_intersect_filter/phmg_it__step5_deseq2_filter.R
cp_one pipeline_candidates/scripts/step6_artifact_filter.py pipeline/05_intersect_filter/phmg_it__step6_artifact_filter.py
cp_one pipeline_candidates/scripts/step7_pop_gp_filter.py   pipeline/05_intersect_filter/phmg_it__step7_pop_gp_filter.py
cp_one mammary_cancer/scripts/06_intersection_rna_edit.sh   pipeline/05_intersect_filter/mammary__06_intersection_rna_edit.sh

# ---- 06_validation ----
cp_one pipeline_candidates/scripts/validation_dna_truth.py  pipeline/06_validation/phmg_it__validation_dna_truth.py
cp_one PHMG_IT/coworker/validation/benchmark.py             pipeline/06_validation/phmg_it__benchmark.py
cp_one PHMG_IT/coworker/validation/run_bcftools_isec.sh     pipeline/06_validation/phmg_it__run_bcftools_isec.sh
cp_one PHMG_IT/scripts/build_dna_truth_treatment_specific.sh pipeline/06_validation/phmg_it__build_dna_truth_treatment_specific.sh
cp_one PHMG_IT/scripts/evaluate_coverage.sh                 pipeline/06_validation/phmg_it__evaluate_coverage.sh
cp_one mammary_cancer/scripts/01_extract_exonic_bed.sh      pipeline/06_validation/mammary__01_extract_exonic_bed.sh
cp_one mammary_cancer/scripts/07_wes_validation.sh          pipeline/06_validation/mammary__07_wes_validation.sh
cp_one mammary_cancer/scripts/08_compare_to_phmg.sh         pipeline/06_validation/mammary__08_compare_to_phmg.sh
cp_one mammary_cancer/scripts/09_sample_identity_check.sh   pipeline/06_validation/mammary__09_sample_identity_check.sh

# ---- analysis: recall_decomposition (스크립트 묶음) ----
for f in 01_build_inputs.sh 02_classify_b1b2b3.sh 02b_reparse.sh 03_b4_kmer_scan.sh \
         03b_recover_p3p9_unmapped.sh 04_compute_recall.py 05_b4_salvage_to_numerator.py \
         06a_verify_salvage_quick.py 06b_verify_imputation.sh 06c_salvage_with_imputation_filter.py \
         07_augmented_realign.sh 08_final_salvage_recall.py 09_full_imputation_filter.sh; do
  cp_one "PHMG_IT/results/recall_decomposition/scripts/$f" "analysis/recall_decomposition/$f"
done
cp_one PHMG_IT/results/recall_decomposition/data/_kmer_scan.py analysis/recall_decomposition/_kmer_scan.py

# ---- analysis: 기타 특성화 ----
cp_one PHMG_IT/scripts/filter_optimization.sh              analysis/coverage_gene_mapping/phmg_it__filter_optimization.sh
cp_one PHMG_IT/scripts/multi_sample_callable.sh            analysis/coverage_gene_mapping/phmg_it__multi_sample_callable.sh
cp_one PHMG_IT/scripts/pipeline_multi_sample_callable.sh   analysis/coverage_gene_mapping/phmg_it__pipeline_multi_sample_callable.sh

# ---- experimental: superseded arms ----
cp_one PHMG_IT/scripts/freebayes_calling.sh               experimental/freebayes/freebayes_calling.sh
cp_one PHMG_IT/scripts/freebayes_eval_no_chr1.sh          experimental/freebayes/freebayes_eval_no_chr1.sh
cp_one PHMG_IT/scripts/joint_calling.sh                   experimental/joint_calling/joint_calling.sh
cp_one PHMG_IT/scripts/run_arm_B.sh                       experimental/arms_BCD/run_arm_B.sh
cp_one PHMG_IT/scripts/run_arm_C.sh                       experimental/arms_BCD/run_arm_C.sh
cp_one PHMG_IT/scripts/run_arm_C_batch.sh                 experimental/arms_BCD/run_arm_C_batch.sh
cp_one PHMG_IT/scripts/arm_D_alignment.sh                 experimental/arms_BCD/arm_D_alignment.sh
cp_one pipeline_candidates/scripts/step3_discordance_detection.sh experimental/discordance/step3_discordance_detection.sh

# ---- docs (method 문서 복사) ----
cp_one CLAUDE.md                                          docs/00_project_overview__CLAUDE.md
cp_one PHMG_IT/CLAUDE.md                                  docs/phmg_it__CLAUDE.md
cp_one PHMG_IH/CLAUDE.md                                  docs/phmg_ih__CLAUDE.md
cp_one pipeline_candidates/secondary_pipeline.md          docs/method_design__secondary_pipeline.md
cp_one pipeline_candidates/CLAUDE_20260506.md             docs/dev_log_20260506.md
cp_one mammary_cancer/CROSS_VALIDATION_20260528.md        docs/cross_validation_20260528.md
cp_one mammary_cancer/sample_mapping_evidence.md          docs/mammary_sample_mapping_evidence.md
cp_one PHMG_IT/results/recall_decomposition/README.md     docs/recall_decomposition__README.md

echo
echo "복사 완료. MANIFEST 항목 수: $(($(wc -l < "$MAN") - 1))"
