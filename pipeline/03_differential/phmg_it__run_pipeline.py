#!/usr/bin/env python3
"""
run_pipeline.py — Mouse Variant-Expression 통합 분석 파이프라인

사용법:
  python run_pipeline.py --config config.yaml
  python run_pipeline.py --mode independent \
      --target dis1.vcf dis2.vcf dis3.vcf \
      --reference ctrl1.vcf ctrl2.vcf ctrl3.vcf \
      --de-genes de_results.csv \
      --vep-results vep_output.txt \
      --output results/

VCF 파일은 bcftools mpileup/call 또는 GATK HaplotypeCaller 출력 모두 지원합니다.
"""

import argparse
import json
import sys
import os
from pathlib import Path
from datetime import datetime

# Pipeline modules
from modules.vcf_parser import (
    parse_vcf, filter_variants, flag_near_indel,
    split_snp_indel, write_vcf_for_vep, detect_vcf_source,
    records_to_table
)
from modules.variant_analyzer import (
    analyze_independent, analyze_paired,
    to_variant_records, write_results_csv
)
from modules.integration import (
    load_de_genes, load_vep_results, integrate,
    write_candidate_genes
)


def load_yaml_config(filepath: str) -> dict:
    """간단한 YAML 파서 (외부 의존성 없이)"""
    config = {}
    current_section = None

    with open(filepath, 'r') as f:
        for line in f:
            line = line.rstrip()
            if not line or line.startswith('#'):
                continue

            # Top-level key
            if not line.startswith(' ') and ':' in line:
                key, _, value = line.partition(':')
                key = key.strip()
                value = value.strip()
                if value:
                    # Simple key: value
                    config[key] = value
                else:
                    # Section header
                    current_section = key
                    config[current_section] = {}
            elif current_section and line.startswith('  '):
                # Nested key: value
                stripped = line.strip()
                if stripped.startswith('- '):
                    # List item
                    if not isinstance(config[current_section], list):
                        config[current_section] = []
                    config[current_section].append(stripped[2:].strip())
                elif ':' in stripped:
                    key, _, value = stripped.partition(':')
                    config[current_section][key.strip()] = value.strip()

    return config


def run_independent_analysis(args):
    """Strategy A: Independent samples (Disease vs WT)"""
    print("=" * 60)
    print("  Mouse Variant-Expression Pipeline")
    print(f"  Mode: Independent (Disease model vs WT)")
    print(f"  Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("-" * 60)
    print(f"  Parameters:")
    print(f"    min_coverage:        {args.min_coverage}")
    print(f"    min_qual:            {args.min_qual}")
    print(f"    min_alt_ratio:       {args.min_alt_ratio}")
    print(f"    max_ref_alt_ratio:   {args.max_ref_alt_ratio}")
    print(f"    min_delta_alt_ratio: {args.min_delta_alt_ratio}")
    print(f"    min_recurrence:      {args.min_recurrence}")
    print("=" * 60)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # === Step 1: Parse VCF files ===
    print(f"\n{'='*60}")
    print("STEP 1: Parsing VCF files")
    print(f"{'='*60}")

    target_samples = {}
    for vcf_path in args.target:
        name = Path(vcf_path).stem
        records = parse_vcf(vcf_path, min_qual=args.min_qual)
        records = flag_near_indel(records, distance=10)
        target_samples[name] = records

    reference_samples = {}
    for vcf_path in args.reference:
        name = Path(vcf_path).stem
        records = parse_vcf(vcf_path, min_qual=args.min_qual)
        reference_samples[name] = records

    # === Step 2: Differential Variant Analysis ===
    print(f"\n{'='*60}")
    print("STEP 2: Differential Variant Analysis")
    print(f"{'='*60}")

    diff_variants = analyze_independent(
        target_samples=target_samples,
        reference_samples=reference_samples,
        min_coverage=args.min_coverage,
        min_target_alt_ratio=args.min_alt_ratio,
        min_recurrence=args.min_recurrence,
        max_ref_alt_ratio=args.max_ref_alt_ratio,
        min_delta_alt_ratio=args.min_delta_alt_ratio,
    )

    # SNP/INDEL 별도 결과 파일
    snp_variants = [dv for dv in diff_variants if dv.variant_type.value == 'SNP']
    indel_variants = [dv for dv in diff_variants if dv.variant_type.value != 'SNP']

    write_results_csv(snp_variants, str(output_dir / "differential_snps.csv"))
    write_results_csv(indel_variants, str(output_dir / "differential_indels.csv"))
    write_results_csv(diff_variants, str(output_dir / "differential_all_variants.csv"))

    # VEP 입력용 VCF
    vep_records = to_variant_records(diff_variants)
    write_vcf_for_vep(vep_records, str(output_dir / "for_vep.vcf"))

    print(f"\n  >> VEP input written to: {output_dir / 'for_vep.vcf'}")
    print(f"  >> Run VEP:")
    print(f"     vep -i {output_dir / 'for_vep.vcf'} "
          f"-o {output_dir / 'vep_output.txt'} "
          f"--species mus_musculus --cache --sift b --polyphen b --symbol --canonical")

    # === Step 3: Integration (if DE genes and VEP results provided) ===
    if args.de_genes and args.vep_results:
        print(f"\n{'='*60}")
        print("STEP 3: Variant-Expression Integration")
        print(f"{'='*60}")

        de_genes = load_de_genes(args.de_genes,
                                 padj_cutoff=args.padj_cutoff,
                                 log2fc_cutoff=args.log2fc_cutoff)
        vep_annotations = load_vep_results(args.vep_results)

        candidates = integrate(
            diff_variants=diff_variants,
            vep_annotations=vep_annotations,
            de_genes=de_genes,
            padj_cutoff=args.padj_cutoff,
            log2fc_cutoff=args.log2fc_cutoff,
        )

        write_candidate_genes(candidates, str(output_dir / "candidate_genes.csv"))
        # Tier별 분리 출력
        for tier in range(1, 6):
            tier_cands = [c for c in candidates if c.tier == tier]
            if tier_cands:
                write_candidate_genes(
                    tier_cands,
                    str(output_dir / f"candidates_tier{tier}.csv"),
                    max_tier=tier
                )

    elif args.de_genes and not args.vep_results:
        print(f"\n  [NOTE] DE genes provided but VEP results not yet available.")
        print(f"  Run VEP first, then re-run with --vep-results to complete integration.")

    # === Summary ===
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"  Target samples:     {len(target_samples)}")
    print(f"  Reference samples:  {len(reference_samples)}")
    print(f"  Differential SNPs:  {len(snp_variants)}")
    print(f"  Differential INDELs: {len(indel_variants)}")
    print(f"  Output directory:   {output_dir}")
    print(f"{'='*60}\n")


def run_paired_analysis(args):
    """Strategy B: Paired samples (Tumor vs Normal)"""
    print("=" * 60)
    print("  Mouse Variant-Expression Pipeline")
    print(f"  Mode: Paired (Tumor vs Normal)")
    print(f"  Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("-" * 60)
    print(f"  Parameters:")
    print(f"    min_coverage:        {args.min_coverage}")
    print(f"    min_qual:            {args.min_qual}")
    print(f"    min_alt_ratio:       {args.min_alt_ratio}")
    print(f"    max_ref_alt_ratio:   {args.max_ref_alt_ratio}")
    print(f"    min_delta_alt_ratio: {args.min_delta_alt_ratio}")
    print(f"    min_recurrence:      {args.min_recurrence}")
    print("=" * 60)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Parse pairs
    if len(args.target) != len(args.reference):
        print("[ERROR] Paired mode requires equal number of target and reference files.",
              file=sys.stderr)
        sys.exit(1)

    # === Step 1: Parse VCF files ===
    print(f"\n{'='*60}")
    print("STEP 1: Parsing VCF files")
    print(f"{'='*60}")

    tumor_samples = {}
    normal_samples = {}
    pairs = []

    for tumor_path, normal_path in zip(args.target, args.reference):
        tumor_name = Path(tumor_path).stem
        normal_name = Path(normal_path).stem

        tumor_records = parse_vcf(tumor_path, min_qual=args.min_qual)
        tumor_records = flag_near_indel(tumor_records, distance=10)
        tumor_samples[tumor_name] = tumor_records

        normal_records = parse_vcf(normal_path, min_qual=args.min_qual)
        normal_samples[normal_name] = normal_records

        pairs.append((tumor_name, normal_name))

    # === Step 2: Somatic Variant Detection ===
    print(f"\n{'='*60}")
    print("STEP 2: Somatic Variant Detection")
    print(f"{'='*60}")

    diff_variants = analyze_paired(
        tumor_samples=tumor_samples,
        normal_samples=normal_samples,
        pairs=pairs,
        min_coverage=args.min_coverage,
        min_tumor_alt_ratio=args.min_alt_ratio,
        max_normal_alt_ratio=args.max_ref_alt_ratio,
        min_recurrence=args.min_recurrence,
        min_delta_alt_ratio=args.min_delta_alt_ratio,
    )

    # 결과 출력 (이하 independent와 동일)
    snp_variants = [dv for dv in diff_variants if dv.variant_type.value == 'SNP']
    indel_variants = [dv for dv in diff_variants if dv.variant_type.value != 'SNP']

    write_results_csv(snp_variants, str(output_dir / "somatic_snps.csv"))
    write_results_csv(indel_variants, str(output_dir / "somatic_indels.csv"))
    write_results_csv(diff_variants, str(output_dir / "somatic_all_variants.csv"))

    vep_records = to_variant_records(diff_variants)
    write_vcf_for_vep(vep_records, str(output_dir / "for_vep.vcf"))

    # Integration (if available)
    if args.de_genes and args.vep_results:
        print(f"\n{'='*60}")
        print("STEP 3: Variant-Expression Integration")
        print(f"{'='*60}")

        de_genes = load_de_genes(args.de_genes,
                                 padj_cutoff=args.padj_cutoff,
                                 log2fc_cutoff=args.log2fc_cutoff)
        vep_annotations = load_vep_results(args.vep_results)

        candidates = integrate(
            diff_variants=diff_variants,
            vep_annotations=vep_annotations,
            de_genes=de_genes,
            padj_cutoff=args.padj_cutoff,
            log2fc_cutoff=args.log2fc_cutoff,
        )

        write_candidate_genes(candidates, str(output_dir / "candidate_genes.csv"))

    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"  Pairs analyzed:     {len(pairs)}")
    print(f"  Somatic SNPs:       {len(snp_variants)}")
    print(f"  Somatic INDELs:     {len(indel_variants)}")
    print(f"  Output directory:   {output_dir}")
    print(f"{'='*60}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Mouse Variant-Expression Integration Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:

  # Disease model vs WT (independent)
  python run_pipeline.py --mode independent \\
      --target disease_1.vcf disease_2.vcf disease_3.vcf \\
      --reference control_1.vcf control_2.vcf control_3.vcf \\
      --output results/disease_vs_wt/

  # Tumor vs Normal (paired)
  python run_pipeline.py --mode paired \\
      --target tumor_1.vcf tumor_2.vcf tumor_3.vcf \\
      --reference normal_1.vcf normal_2.vcf normal_3.vcf \\
      --output results/tumor_vs_normal/

  # Full pipeline with DE genes + VEP
  python run_pipeline.py --mode independent \\
      --target dis1.vcf dis2.vcf dis3.vcf \\
      --reference ctrl1.vcf ctrl2.vcf ctrl3.vcf \\
      --de-genes deseq2_results.csv \\
      --vep-results vep_output.txt \\
      --output results/

  # With YAML config
  python run_pipeline.py --config config.yaml

VCF files from both bcftools and GATK HaplotypeCaller are supported.
The pipeline automatically detects the VCF source format.
        """
    )

    parser.add_argument('--config', type=str,
                        help='YAML config file (overrides other arguments)')
    parser.add_argument('--mode', choices=['independent', 'paired'],
                        default='independent',
                        help='Analysis mode: independent (disease vs WT) '
                             'or paired (tumor vs normal)')

    # Input files
    parser.add_argument('--target', nargs='+',
                        help='Target VCF files (disease or tumor)')
    parser.add_argument('--reference', nargs='+',
                        help='Reference VCF files (control or normal)')
    parser.add_argument('--de-genes', type=str,
                        help='DESeq2 results CSV (optional for step 3)')
    parser.add_argument('--vep-results', type=str,
                        help='VEP output file (optional for step 3)')

    # Output
    parser.add_argument('--output', type=str, default='results/',
                        help='Output directory')

    # Filtering parameters
    parser.add_argument('--min-coverage', type=int, default=10,
                        help='Minimum read depth (default: 10)')
    parser.add_argument('--min-qual', type=float, default=30.0,
                        help='Minimum QUAL score (default: 30)')
    parser.add_argument('--min-alt-ratio', type=float, default=0.3,
                        help='Min alt allele ratio in target (default: 0.3)')
    parser.add_argument('--max-ref-alt-ratio', type=float, default=0.1,
                        help='Max alt allele ratio in reference (default: 0.1)')
    parser.add_argument('--min-recurrence', type=int, default=2,
                        help='Min samples with variant (default: 2)')
    parser.add_argument('--min-delta-alt-ratio', type=float, default=0.2,
                        help='Min alt ratio difference between target and reference (default: 0.2)')

    # DE cutoffs
    parser.add_argument('--padj-cutoff', type=float, default=0.05,
                        help='Adjusted p-value cutoff for DE (default: 0.05)')
    parser.add_argument('--log2fc-cutoff', type=float, default=1.0,
                        help='|log2FC| cutoff for DE (default: 1.0)')

    args = parser.parse_args()

    # Config file override
    if args.config:
        config = load_yaml_config(args.config)
        # Map config to args (simplified)
        if 'mode' in config:
            args.mode = config['mode']
        if 'output' in config:
            args.output = config['output']
        if isinstance(config.get('target'), list):
            args.target = config['target']
        if isinstance(config.get('reference'), list):
            args.reference = config['reference']
        if isinstance(config.get('params'), dict):
            params = config['params']
            for key in ['min_coverage', 'min_qual', 'min_recurrence']:
                if key in params:
                    setattr(args, key.replace('-', '_'), int(params[key]))
            for key in ['min_alt_ratio', 'max_ref_alt_ratio',
                        'min_delta_alt_ratio',
                        'padj_cutoff', 'log2fc_cutoff']:
                if key in params:
                    setattr(args, key.replace('-', '_'), float(params[key]))
        if 'de_genes' in config:
            args.de_genes = config['de_genes']
        if 'vep_results' in config:
            args.vep_results = config['vep_results']

    # Validation
    if not args.target or not args.reference:
        parser.error("--target and --reference VCF files are required")

    # Run
    if args.mode == 'independent':
        run_independent_analysis(args)
    elif args.mode == 'paired':
        run_paired_analysis(args)


if __name__ == "__main__":
    main()
