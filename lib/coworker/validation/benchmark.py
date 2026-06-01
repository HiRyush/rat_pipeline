#!/usr/bin/env python3
"""
benchmark.py — variant_analyzer vs bcftools isec 정량 비교

Workflow:
  1) 동일 WGS VCF에 bcftools isec 실행 → gold standard (condition-specific variants)
  2) variant_analyzer 결과 (differential_all_variants.csv) 로드
  3) Position-level 비교 → TP, FP, FN 산출
  4) Sensitivity, Precision, F1, Jaccard 계산
  5) Chromosome별/Variant type별 성능 breakdown
  6) 결과 리포트 + CSV 출력

사용법:
  # Step 1: bcftools isec 실행 (bash)
  bash run_bcftools_isec.sh

  # Step 2: 벤치마킹
  python benchmark.py \
      --analyzer-csv differential_all_variants.csv \
      --isec-dir bcftools_isec_output/ \
      --output benchmark_results/

  # 또는 isec 대신 disease-only VCF 직접 지정
  python benchmark.py \
      --analyzer-csv differential_all_variants.csv \
      --disease-only-vcf bcftools_isec_output/0000.vcf \
      --shared-vcf bcftools_isec_output/0002.vcf \
      --output benchmark_results/
"""

import argparse
import csv
import gzip
import sys
import os
import re
from collections import defaultdict
from dataclasses import dataclass, field
from typing import List, Dict, Set, Tuple, Optional
from pathlib import Path


# ============================================================
# Data structures
# ============================================================

@dataclass
class BenchmarkVariant:
    chrom: str
    pos: int
    ref: str = ""
    alt: str = ""

    @property
    def position_key(self) -> str:
        """Position + allele 기반 매칭"""
        return f"{self.chrom}_{self.pos}_{self.alt}"

    @property
    def locus_key(self) -> str:
        """Position만으로 매칭 (allele 무관)"""
        return f"{self.chrom}_{self.pos}"


@dataclass
class BenchmarkMetrics:
    """성능 지표"""
    tp: int = 0          # True Positive: 둘 다 검출
    fp: int = 0          # False Positive: analyzer만 검출 (isec에 없음)
    fn: int = 0          # False Negative: isec만 검출 (analyzer에 없음)
    tn: int = 0          # True Negative: 둘 다 미검출 (해당 없음)

    @property
    def sensitivity(self) -> float:
        """Recall = TP / (TP + FN)"""
        return self.tp / (self.tp + self.fn) if (self.tp + self.fn) > 0 else 0.0

    @property
    def precision(self) -> float:
        """Precision = TP / (TP + FP)"""
        return self.tp / (self.tp + self.fp) if (self.tp + self.fp) > 0 else 0.0

    @property
    def f1_score(self) -> float:
        """F1 = 2 * (P * R) / (P + R)"""
        p, r = self.precision, self.sensitivity
        return 2 * p * r / (p + r) if (p + r) > 0 else 0.0

    @property
    def jaccard(self) -> float:
        """Jaccard Index = TP / (TP + FP + FN)"""
        total = self.tp + self.fp + self.fn
        return self.tp / total if total > 0 else 0.0

    def summary_str(self) -> str:
        return (f"TP={self.tp}  FP={self.fp}  FN={self.fn}  |  "
                f"Sensitivity={self.sensitivity:.4f}  "
                f"Precision={self.precision:.4f}  "
                f"F1={self.f1_score:.4f}  "
                f"Jaccard={self.jaccard:.4f}")


# ============================================================
# Parsing functions
# ============================================================

def parse_analyzer_csv(filepath: str) -> List[BenchmarkVariant]:
    """variant_analyzer의 differential_all_variants.csv 파싱"""
    variants = []
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            chrom = row['chrom'].strip()
            # chr prefix 통일
            if not chrom.startswith('chr'):
                chrom = f"chr{chrom}"
            try:
                pos = int(row['pos'])
            except ValueError:
                continue
            ref = row.get('ref', '').strip()
            alt = row.get('alt', '').strip()
            variants.append(BenchmarkVariant(chrom=chrom, pos=pos, ref=ref, alt=alt))

    print(f"[Analyzer] Loaded {len(variants)} variants from {filepath}",
          file=sys.stderr)
    return variants


def _open_vcf(filepath: str):
    """VCF 또는 VCF.GZ 파일을 자동 감지하여 열기"""
    if filepath.endswith('.gz') or filepath.endswith('.bgz'):
        return gzip.open(filepath, 'rt', encoding='utf-8')
    else:
        return open(filepath, 'r', encoding='utf-8')


def parse_vcf_variants(filepath: str) -> List[BenchmarkVariant]:
    """VCF / VCF.GZ 파일에서 variant 위치 추출 (bcftools isec 출력 등)"""
    variants = []
    with _open_vcf(filepath) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue

            chrom = fields[0].strip()
            if not chrom.startswith('chr'):
                chrom = f"chr{chrom}"

            try:
                pos = int(fields[1])
            except ValueError:
                continue

            ref = fields[3].strip()
            alt_field = fields[4].strip()

            for alt in alt_field.split(','):
                alt = alt.strip()
                if alt == '.' or alt == '*':
                    continue
                variants.append(BenchmarkVariant(
                    chrom=chrom, pos=pos, ref=ref, alt=alt))

    print(f"[VCF] Loaded {len(variants)} variants from {filepath}",
          file=sys.stderr)
    return variants


def parse_isec_directory(isec_dir: str) -> Dict[str, List[BenchmarkVariant]]:
    """
    bcftools isec 출력 디렉토리 파싱

    bcftools isec -p output_dir disease.vcf.gz control.vcf.gz
    출력:
      0000.vcf — disease에만 있는 variant
      0001.vcf — control에만 있는 variant
      0002.vcf — 공통 variant (disease 기준)
      0003.vcf — 공통 variant (control 기준)

    bcftools isec -p output_dir control.vcf.gz disease.vcf.gz (순서 반대)
    출력:
      0000.vcf — control에만 있는 variant
      0001.vcf — disease에만 있는 variant
      0002.vcf — 공통 variant (control 기준)
      0003.vcf — 공통 variant (disease 기준)
    """
    result = {}

    for fname in ['0000.vcf', '0001.vcf', '0002.vcf', '0003.vcf',
                   '0000.vcf.gz', '0001.vcf.gz', '0002.vcf.gz', '0003.vcf.gz']:
        fpath = os.path.join(isec_dir, fname)
        if os.path.exists(fpath):
            key = fname.replace('.vcf.gz', '').replace('.vcf', '')
            result[key] = parse_vcf_variants(fpath)

    # README.txt에서 파일 매핑 정보 추출
    readme_path = os.path.join(isec_dir, 'README.txt')
    if os.path.exists(readme_path):
        print(f"\n[isec README]", file=sys.stderr)
        with open(readme_path, 'r') as f:
            for line in f:
                print(f"  {line.rstrip()}", file=sys.stderr)
        print(file=sys.stderr)

    return result


# ============================================================
# Comparison logic
# ============================================================

def compare_variant_sets(analyzer_variants: List[BenchmarkVariant],
                         gold_standard: List[BenchmarkVariant],
                         match_mode: str = "position_allele"
                         ) -> Tuple[BenchmarkMetrics, Dict]:
    """
    두 variant set 비교

    match_mode:
      "position_allele" — chrom + pos + alt가 모두 일치해야 match (엄격)
      "position_only"   — chrom + pos만 일치하면 match (관대)
    """
    if match_mode == "position_allele":
        analyzer_set = {v.position_key for v in analyzer_variants}
        gold_set = {v.position_key for v in gold_standard}
    elif match_mode == "position_only":
        analyzer_set = {v.locus_key for v in analyzer_variants}
        gold_set = {v.locus_key for v in gold_standard}
    else:
        raise ValueError(f"Unknown match_mode: {match_mode}")

    tp_set = analyzer_set & gold_set
    fp_set = analyzer_set - gold_set
    fn_set = gold_set - analyzer_set

    metrics = BenchmarkMetrics(
        tp=len(tp_set),
        fp=len(fp_set),
        fn=len(fn_set),
    )

    details = {
        'tp_variants': sorted(tp_set),
        'fp_variants': sorted(fp_set),
        'fn_variants': sorted(fn_set),
    }

    return metrics, details


def compare_by_chromosome(analyzer_variants: List[BenchmarkVariant],
                          gold_standard: List[BenchmarkVariant],
                          match_mode: str = "position_allele"
                          ) -> Dict[str, BenchmarkMetrics]:
    """Chromosome별 성능 분석"""
    # Chromosome별 분류
    analyzer_by_chr = defaultdict(list)
    gold_by_chr = defaultdict(list)

    for v in analyzer_variants:
        analyzer_by_chr[v.chrom].append(v)
    for v in gold_standard:
        gold_by_chr[v.chrom].append(v)

    all_chroms = sorted(set(list(analyzer_by_chr.keys()) +
                            list(gold_by_chr.keys())),
                        key=lambda x: (len(x), x))

    results = {}
    for chrom in all_chroms:
        metrics, _ = compare_variant_sets(
            analyzer_by_chr.get(chrom, []),
            gold_by_chr.get(chrom, []),
            match_mode=match_mode
        )
        results[chrom] = metrics

    return results


def compare_by_variant_type(analyzer_variants: List[BenchmarkVariant],
                            gold_standard: List[BenchmarkVariant],
                            analyzer_csv_path: str,
                            match_mode: str = "position_allele"
                            ) -> Dict[str, BenchmarkMetrics]:
    """SNP vs INDEL 별 성능 분석"""
    # variant_analyzer CSV에서 type 정보 가져오기
    type_map = {}
    with open(analyzer_csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            chrom = row['chrom'].strip()
            if not chrom.startswith('chr'):
                chrom = f"chr{chrom}"
            key = f"{chrom}_{row['pos']}_{row['alt']}"
            type_map[key] = row.get('variant_type', 'SNP')

    # Gold standard에서 type 판별 (REF/ALT 길이)
    def get_type(v):
        if v.position_key in type_map:
            return type_map[v.position_key]
        if len(v.ref) == len(v.alt):
            return 'SNP'
        return 'INDEL'

    # Type별 분류
    analyzer_snp = [v for v in analyzer_variants if get_type(v) == 'SNP']
    analyzer_indel = [v for v in analyzer_variants if get_type(v) != 'SNP']
    gold_snp = [v for v in gold_standard if len(v.ref) == len(v.alt) == 1]
    gold_indel = [v for v in gold_standard if len(v.ref) != 1 or len(v.alt) != 1]

    results = {}
    snp_metrics, _ = compare_variant_sets(analyzer_snp, gold_snp, match_mode)
    indel_metrics, _ = compare_variant_sets(analyzer_indel, gold_indel, match_mode)
    results['SNP'] = snp_metrics
    results['INDEL'] = indel_metrics

    return results


# ============================================================
# Output generation
# ============================================================

def write_report(metrics: BenchmarkMetrics,
                 details: Dict,
                 chr_metrics: Dict[str, BenchmarkMetrics],
                 type_metrics: Dict[str, BenchmarkMetrics],
                 output_dir: str,
                 match_mode: str,
                 analyzer_total: int,
                 gold_total: int):
    """종합 리포트 생성"""
    os.makedirs(output_dir, exist_ok=True)

    report_path = os.path.join(output_dir, "benchmark_report.txt")
    with open(report_path, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("  BENCHMARK REPORT: variant_analyzer vs bcftools isec\n")
        f.write("=" * 70 + "\n\n")

        f.write(f"Match mode: {match_mode}\n")
        f.write(f"Analyzer variants:     {analyzer_total}\n")
        f.write(f"Gold standard variants: {gold_total}\n\n")

        f.write("-" * 70 + "\n")
        f.write("OVERALL METRICS\n")
        f.write("-" * 70 + "\n")
        f.write(f"  True Positives  (TP): {metrics.tp:>6}  "
                f"(both detected)\n")
        f.write(f"  False Positives (FP): {metrics.fp:>6}  "
                f"(analyzer only — not in isec)\n")
        f.write(f"  False Negatives (FN): {metrics.fn:>6}  "
                f"(isec only — missed by analyzer)\n\n")
        f.write(f"  Sensitivity (Recall):  {metrics.sensitivity:.4f}  "
                f"(= TP / (TP+FN) = {metrics.tp}/{metrics.tp+metrics.fn})\n")
        f.write(f"  Precision:             {metrics.precision:.4f}  "
                f"(= TP / (TP+FP) = {metrics.tp}/{metrics.tp+metrics.fp})\n")
        f.write(f"  F1 Score:              {metrics.f1_score:.4f}\n")
        f.write(f"  Jaccard Index:         {metrics.jaccard:.4f}\n\n")

        # Interpretation
        f.write("-" * 70 + "\n")
        f.write("INTERPRETATION\n")
        f.write("-" * 70 + "\n")
        if metrics.sensitivity >= 0.8:
            f.write("  Sensitivity >= 0.80: GOOD — 대부분의 gold standard variant를 검출함\n")
        elif metrics.sensitivity >= 0.5:
            f.write("  Sensitivity 0.50-0.79: MODERATE — 일부 variant를 놓침\n")
        else:
            f.write("  Sensitivity < 0.50: LOW — 상당수 variant를 놓침 (filtering 검토 필요)\n")

        if metrics.precision >= 0.8:
            f.write("  Precision >= 0.80: GOOD — false positive가 적음\n")
        elif metrics.precision >= 0.5:
            f.write("  Precision 0.50-0.79: MODERATE — 일부 false positive 존재\n")
        else:
            f.write("  Precision < 0.50: LOW — false positive가 많음 (filtering 강화 필요)\n")
        f.write("\n")

        # By variant type
        f.write("-" * 70 + "\n")
        f.write("BY VARIANT TYPE\n")
        f.write("-" * 70 + "\n")
        f.write(f"  {'Type':<10} {'TP':>6} {'FP':>6} {'FN':>6}  "
                f"{'Sens':>8} {'Prec':>8} {'F1':>8}\n")
        f.write(f"  {'-'*10} {'-'*6} {'-'*6} {'-'*6}  "
                f"{'-'*8} {'-'*8} {'-'*8}\n")
        for vtype, m in type_metrics.items():
            f.write(f"  {vtype:<10} {m.tp:>6} {m.fp:>6} {m.fn:>6}  "
                    f"{m.sensitivity:>8.4f} {m.precision:>8.4f} {m.f1_score:>8.4f}\n")
        f.write("\n")

        # By chromosome
        f.write("-" * 70 + "\n")
        f.write("BY CHROMOSOME\n")
        f.write("-" * 70 + "\n")
        f.write(f"  {'Chr':<8} {'TP':>6} {'FP':>6} {'FN':>6}  "
                f"{'Sens':>8} {'Prec':>8} {'F1':>8}\n")
        f.write(f"  {'-'*8} {'-'*6} {'-'*6} {'-'*6}  "
                f"{'-'*8} {'-'*8} {'-'*8}\n")
        for chrom, m in chr_metrics.items():
            if m.tp + m.fp + m.fn == 0:
                continue
            f.write(f"  {chrom:<8} {m.tp:>6} {m.fp:>6} {m.fn:>6}  "
                    f"{m.sensitivity:>8.4f} {m.precision:>8.4f} {m.f1_score:>8.4f}\n")
        f.write("\n")

        f.write("=" * 70 + "\n")

    print(f"\n[Report] Written to {report_path}", file=sys.stderr)

    # Detailed CSV outputs
    # TP variants
    tp_path = os.path.join(output_dir, "true_positives.csv")
    with open(tp_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['variant_key', 'status'])
        for v in details['tp_variants']:
            writer.writerow([v, 'TP'])

    # FP variants (analyzer only)
    fp_path = os.path.join(output_dir, "false_positives.csv")
    with open(fp_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['variant_key', 'status', 'note'])
        for v in details['fp_variants']:
            writer.writerow([v, 'FP', 'detected by analyzer but not by isec'])

    # FN variants (isec only)
    fn_path = os.path.join(output_dir, "false_negatives.csv")
    with open(fn_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['variant_key', 'status', 'note'])
        for v in details['fn_variants']:
            writer.writerow([v, 'FN', 'detected by isec but not by analyzer'])

    print(f"[Output] TP: {tp_path}", file=sys.stderr)
    print(f"[Output] FP: {fp_path}", file=sys.stderr)
    print(f"[Output] FN: {fn_path}", file=sys.stderr)

    return report_path


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Benchmark variant_analyzer vs bcftools isec",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:

  # bcftools isec 디렉토리 사용
  python benchmark.py \\
      --analyzer-csv differential_all_variants.csv \\
      --isec-dir isec_output/ \\
      --isec-disease-index 0 \\
      --output benchmark_results/

  # 개별 VCF 파일 지정
  python benchmark.py \\
      --analyzer-csv differential_all_variants.csv \\
      --disease-only-vcf isec_output/0000.vcf \\
      --output benchmark_results/

  # Position-only matching (allele 무관)
  python benchmark.py \\
      --analyzer-csv differential_all_variants.csv \\
      --disease-only-vcf isec_output/0000.vcf \\
      --match-mode position_only \\
      --output benchmark_results/

Match modes:
  position_allele — chrom + pos + alt 모두 일치 (엄격, 기본값)
  position_only   — chrom + pos만 일치 (관대, multi-allelic site 허용)
        """
    )

    parser.add_argument('--analyzer-csv', required=True,
                        help='variant_analyzer output CSV '
                             '(differential_all_variants.csv)')

    # Gold standard input (택 1)
    gold_group = parser.add_mutually_exclusive_group()
    gold_group.add_argument('--isec-dir',
                            help='bcftools isec output directory')
    gold_group.add_argument('--disease-only-vcf',
                            help='Disease-only VCF '
                                 '(bcftools isec의 0000.vcf 또는 0001.vcf)')

    parser.add_argument('--isec-disease-index', type=int, default=0,
                        choices=[0, 1],
                        help='isec 출력에서 disease-only VCF 인덱스 '
                             '(0: disease가 첫 번째 입력, 1: 두 번째 입력)')
    parser.add_argument('--shared-vcf',
                        help='Shared variants VCF (optional, for context)')

    parser.add_argument('--match-mode', default='position_allele',
                        choices=['position_allele', 'position_only'],
                        help='Matching 기준 (default: position_allele)')
    parser.add_argument('--output', default='benchmark_results/',
                        help='Output directory')

    args = parser.parse_args()

    # Load analyzer results
    analyzer_variants = parse_analyzer_csv(args.analyzer_csv)

    # Load gold standard
    gold_standard = None

    if args.isec_dir:
        isec_data = parse_isec_directory(args.isec_dir)
        # disease-only variants
        disease_key = f"000{args.isec_disease_index}"
        if disease_key in isec_data:
            gold_standard = isec_data[disease_key]
            print(f"[Gold Standard] Using {disease_key}.vcf "
                  f"as disease-only variants", file=sys.stderr)
        else:
            print(f"[ERROR] {disease_key}.vcf not found in {args.isec_dir}",
                  file=sys.stderr)
            sys.exit(1)

    elif args.disease_only_vcf:
        gold_standard = parse_vcf_variants(args.disease_only_vcf)

    else:
        print("[ERROR] Provide either --isec-dir or --disease-only-vcf",
              file=sys.stderr)
        sys.exit(1)

    if not gold_standard:
        print("[ERROR] No gold standard variants loaded", file=sys.stderr)
        sys.exit(1)

    # Run comparison
    print(f"\n{'='*60}", file=sys.stderr)
    print(f"  Benchmarking: {len(analyzer_variants)} analyzer variants "
          f"vs {len(gold_standard)} gold standard variants", file=sys.stderr)
    print(f"  Match mode: {args.match_mode}", file=sys.stderr)
    print(f"{'='*60}\n", file=sys.stderr)

    # Overall metrics
    metrics, details = compare_variant_sets(
        analyzer_variants, gold_standard, args.match_mode)
    print(f"  Overall: {metrics.summary_str()}", file=sys.stderr)

    # By chromosome
    chr_metrics = compare_by_chromosome(
        analyzer_variants, gold_standard, args.match_mode)

    # By variant type
    type_metrics = compare_by_variant_type(
        analyzer_variants, gold_standard,
        args.analyzer_csv, args.match_mode)

    for vtype, m in type_metrics.items():
        print(f"  {vtype}: {m.summary_str()}", file=sys.stderr)

    # Generate report
    report_path = write_report(
        metrics=metrics,
        details=details,
        chr_metrics=chr_metrics,
        type_metrics=type_metrics,
        output_dir=args.output,
        match_mode=args.match_mode,
        analyzer_total=len(analyzer_variants),
        gold_total=len(gold_standard),
    )

    # Print summary to stdout
    print(f"\nSensitivity: {metrics.sensitivity:.4f}")
    print(f"Precision:   {metrics.precision:.4f}")
    print(f"F1 Score:    {metrics.f1_score:.4f}")
    print(f"Jaccard:     {metrics.jaccard:.4f}")


if __name__ == "__main__":
    main()
