"""
variant_analyzer.py — Condition-specific variant detection

두 조건 (Control vs Disease / Normal vs Tumor) 간 variant 차이를 분석합니다.

분석 전략:
  A) Disease model vs WT → Independent samples → Recurrence 기반
  B) Tumor vs Normal → Paired samples → 개체 내 비교

SNP와 INDEL을 별도 트랙으로 분석하되, 동일한 로직을 적용합니다.
"""

import csv
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple
from pathlib import Path
from scipy import stats as scipy_stats

from modules.vcf_parser import (
    VariantRecord, VariantType, VCFSource,
    parse_vcf, filter_variants, flag_near_indel,
    split_snp_indel, write_vcf_for_vep
)


@dataclass
class DifferentialVariant:
    """조건 간 비교 결과를 담는 구조"""
    chrom: str
    pos: int
    ref: str
    alt: str
    variant_type: VariantType
    indel_size: int = 0
    is_frameshift: bool = False
    near_indel: bool = False
    # Target (disease/tumor) 쪽 정보
    target_alt_ratios: List[float] = field(default_factory=list)
    target_depths: List[int] = field(default_factory=list)
    target_samples_found: int = 0
    target_mean_alt_ratio: float = 0.0
    # Reference (control/normal) 쪽 정보
    ref_alt_ratios: List[float] = field(default_factory=list)
    ref_depths: List[int] = field(default_factory=list)
    ref_samples_found: int = 0
    ref_mean_alt_ratio: float = 0.0
    # 통계
    recurrence: int = 0
    fisher_pvalue: float = 1.0
    delta_alt_ratio: float = 0.0  # target - reference

    @property
    def position_key(self) -> str:
        return f"{self.chrom}_{self.pos}_{self.alt}"


def _build_position_index(records: List[VariantRecord]) -> Dict[str, VariantRecord]:
    """position_key → VariantRecord 매핑"""
    index = {}
    for r in records:
        index[r.position_key] = r
    return index


def _build_locus_index(records: List[VariantRecord]) -> Dict[str, List[VariantRecord]]:
    """locus_key → [VariantRecord, ...] 매핑 (같은 위치의 여러 variant)"""
    index = defaultdict(list)
    for r in records:
        index[r.locus_key].append(r)
    return index


def analyze_independent(target_samples: Dict[str, List[VariantRecord]],
                        reference_samples: Dict[str, List[VariantRecord]],
                        min_coverage: int = 10,
                        min_target_alt_ratio: float = 0.3,
                        min_recurrence: int = 2,
                        max_ref_alt_ratio: float = 0.1,
                        min_delta_alt_ratio: float = 0.2
                        ) -> List[DifferentialVariant]:
    """
    Strategy A: Independent samples (Disease model vs WT)

    각 target sample에서 발견되고, reference에서는 없거나 낮은 빈도인 variant를 찾습니다.
    여러 target sample에서 반복 출현하는 variant에 높은 신뢰도를 부여합니다.

    Parameters:
        target_samples: {"sample_name": [VariantRecord, ...]} — disease samples
        reference_samples: {"sample_name": [VariantRecord, ...]} — control samples
        min_coverage: 최소 read depth
        min_target_alt_ratio: target에서의 최소 alt allele 비율
        min_recurrence: 최소 몇 개 target sample에서 발견되어야 하는지
        max_ref_alt_ratio: reference에서 이 비율 이하여야 condition-specific으로 판정
        min_delta_alt_ratio: target - reference 간 alt ratio 최소 차이
    """
    print(f"\n[Analyze Independent] {len(target_samples)} target vs "
          f"{len(reference_samples)} reference samples", file=sys.stderr)
    print(f"  Params: coverage>={min_coverage}, target_alt>={min_target_alt_ratio}, "
          f"ref_alt<={max_ref_alt_ratio}, delta>={min_delta_alt_ratio}, "
          f"recurrence>={min_recurrence}", file=sys.stderr)

    # Step 1: Target samples에서 quality filter 통과한 variant 수집
    target_variants_by_pos = defaultdict(list)  # position_key → [(sample, record)]

    for sample_name, records in target_samples.items():
        filtered = filter_variants(records,
                                   min_coverage=min_coverage,
                                   min_alt_ratio=min_target_alt_ratio)
        for r in filtered:
            target_variants_by_pos[r.position_key].append((sample_name, r))

    print(f"  Unique variant positions in target: "
          f"{len(target_variants_by_pos)}", file=sys.stderr)

    # Step 2: Recurrence filter
    recurring = {k: v for k, v in target_variants_by_pos.items()
                 if len(v) >= min_recurrence}

    print(f"  After recurrence filter (>={min_recurrence}): "
          f"{len(recurring)}", file=sys.stderr)

    # Step 3: Reference index 구축
    ref_indices = {}
    for sample_name, records in reference_samples.items():
        ref_indices[sample_name] = _build_position_index(records)

    # Step 4: Reference에서 확인하고 condition-specific variant 선별
    results = []

    for pos_key, target_hits in recurring.items():
        sample_names, target_records = zip(*target_hits)
        representative = target_records[0]

        # Reference에서 같은 위치 확인
        ref_ratios = []
        ref_depths = []
        ref_found = 0

        for ref_name, ref_idx in ref_indices.items():
            if pos_key in ref_idx:
                ref_rec = ref_idx[pos_key]
                if ref_rec.total_depth >= min_coverage:
                    ref_ratios.append(ref_rec.alt_ratio)
                    ref_depths.append(ref_rec.total_depth)
                    ref_found += 1

        # Reference에서의 평균 alt ratio
        mean_ref_ratio = sum(ref_ratios) / len(ref_ratios) if ref_ratios else 0.0

        # Condition-specific 판정
        if mean_ref_ratio > max_ref_alt_ratio:
            continue  # Reference에서도 흔한 variant → skip

        # Target 정보 집계
        target_ratios = [r.alt_ratio for r in target_records]
        target_depths_list = [r.total_depth for r in target_records]
        mean_target_ratio = sum(target_ratios) / len(target_ratios)

        # Delta alt ratio filter (target - reference)
        delta = mean_target_ratio - mean_ref_ratio
        if delta < min_delta_alt_ratio:
            continue  # target-reference difference insufficient

        # Fisher's exact test (optional, coverage가 충분할 때)
        pvalue = 1.0
        if ref_ratios and target_ratios:
            # 대표값으로 2x2 table 구성
            t_alt = sum(r.alt_count for r in target_records)
            t_ref = sum(r.ref_count for r in target_records)
            r_alt = sum(ref_indices[rn][pos_key].alt_count
                        for rn in ref_indices
                        if pos_key in ref_indices[rn]
                        and ref_indices[rn][pos_key].total_depth >= min_coverage)
            r_ref = sum(ref_indices[rn][pos_key].ref_count
                        for rn in ref_indices
                        if pos_key in ref_indices[rn]
                        and ref_indices[rn][pos_key].total_depth >= min_coverage)
            if t_alt + t_ref > 0 and r_alt + r_ref > 0:
                try:
                    _, pvalue = scipy_stats.fisher_exact(
                        [[t_alt, t_ref], [r_alt, r_ref]])
                except Exception:
                    pvalue = 1.0

        dv = DifferentialVariant(
            chrom=representative.chrom,
            pos=representative.pos,
            ref=representative.ref,
            alt=representative.alt,
            variant_type=representative.variant_type,
            indel_size=representative.indel_size,
            is_frameshift=representative.is_frameshift,
            near_indel=representative.near_indel,
            target_alt_ratios=target_ratios,
            target_depths=target_depths_list,
            target_samples_found=len(target_hits),
            target_mean_alt_ratio=mean_target_ratio,
            ref_alt_ratios=ref_ratios,
            ref_depths=ref_depths,
            ref_samples_found=ref_found,
            ref_mean_alt_ratio=mean_ref_ratio,
            recurrence=len(target_hits),
            fisher_pvalue=pvalue,
            delta_alt_ratio=mean_target_ratio - mean_ref_ratio,
        )
        results.append(dv)

    results.sort(key=lambda x: (x.fisher_pvalue, -x.recurrence))

    snp_count = sum(1 for r in results if r.variant_type == VariantType.SNP)
    indel_count = sum(1 for r in results if r.variant_type in
                      (VariantType.INSERTION, VariantType.DELETION))
    print(f"  Condition-specific variants: {len(results)} "
          f"(SNPs: {snp_count}, INDELs: {indel_count})", file=sys.stderr)

    return results


def analyze_paired(tumor_samples: Dict[str, List[VariantRecord]],
                   normal_samples: Dict[str, List[VariantRecord]],
                   pairs: List[Tuple[str, str]],
                   min_coverage: int = 10,
                   min_tumor_alt_ratio: float = 0.1,
                   max_normal_alt_ratio: float = 0.05,
                   min_recurrence: int = 1,
                   min_delta_alt_ratio: float = 0.1
                   ) -> List[DifferentialVariant]:
    """
    Strategy B: Paired samples (Tumor vs Normal, same individual)

    각 pair 내에서 tumor에만 있는 somatic variant를 찾고,
    여러 pair에서 반복되는 것을 우선합니다.

    Parameters:
        tumor_samples: {"sample_name": [VariantRecord, ...]}
        normal_samples: {"sample_name": [VariantRecord, ...]}
        pairs: [(tumor_name, normal_name), ...] — 매칭 정보
        min_tumor_alt_ratio: tumor에서의 최소 alt ratio (subclonal 허용 위해 낮게)
        max_normal_alt_ratio: normal에서 이 이하여야 somatic
    """
    print(f"\n[Analyze Paired] {len(pairs)} pairs", file=sys.stderr)

    somatic_by_pos = defaultdict(list)  # position_key → [(pair_idx, tumor_record)]

    for pair_idx, (tumor_name, normal_name) in enumerate(pairs):
        tumor_records = tumor_samples.get(tumor_name, [])
        normal_records = normal_samples.get(normal_name, [])

        tumor_filtered = filter_variants(tumor_records,
                                         min_coverage=min_coverage,
                                         min_alt_ratio=min_tumor_alt_ratio)
        normal_index = _build_position_index(normal_records)

        for t_rec in tumor_filtered:
            is_somatic = False

            if t_rec.position_key not in normal_index:
                # Normal에 아예 없음 → somatic candidate
                is_somatic = True
            else:
                n_rec = normal_index[t_rec.position_key]
                if n_rec.total_depth >= min_coverage:
                    if n_rec.alt_ratio <= max_normal_alt_ratio:
                        # Normal에서 거의 없음 → somatic candidate
                        is_somatic = True
                else:
                    # Normal coverage 부족 → 판단 불가, 보수적으로 포함
                    is_somatic = True

            if is_somatic:
                somatic_by_pos[t_rec.position_key].append((pair_idx, t_rec))

    print(f"  Unique somatic positions: {len(somatic_by_pos)}", file=sys.stderr)

    # Recurrence filter & 결과 구성
    results = []
    for pos_key, hits in somatic_by_pos.items():
        if len(hits) < min_recurrence:
            continue

        pair_indices, tumor_records = zip(*hits)
        representative = tumor_records[0]

        # Normal 쪽 정보 수집
        ref_ratios = []
        ref_depths = []
        for _, (tumor_name, normal_name) in enumerate(pairs):
            normal_records = normal_samples.get(normal_name, [])
            normal_index = _build_position_index(normal_records)
            if pos_key in normal_index:
                n_rec = normal_index[pos_key]
                ref_ratios.append(n_rec.alt_ratio)
                ref_depths.append(n_rec.total_depth)

        target_ratios = [r.alt_ratio for r in tumor_records]
        target_depths_list = [r.total_depth for r in tumor_records]

        dv = DifferentialVariant(
            chrom=representative.chrom,
            pos=representative.pos,
            ref=representative.ref,
            alt=representative.alt,
            variant_type=representative.variant_type,
            indel_size=representative.indel_size,
            is_frameshift=representative.is_frameshift,
            near_indel=representative.near_indel,
            target_alt_ratios=target_ratios,
            target_depths=target_depths_list,
            target_samples_found=len(hits),
            target_mean_alt_ratio=(sum(target_ratios) / len(target_ratios)),
            ref_alt_ratios=ref_ratios,
            ref_depths=ref_depths,
            ref_samples_found=len(ref_ratios),
            ref_mean_alt_ratio=(sum(ref_ratios) / len(ref_ratios)
                                if ref_ratios else 0.0),
            recurrence=len(hits),
            fisher_pvalue=1.0,
            delta_alt_ratio=0.0,
        )
        dv.delta_alt_ratio = dv.target_mean_alt_ratio - dv.ref_mean_alt_ratio

        # Delta alt ratio filter
        if dv.delta_alt_ratio < min_delta_alt_ratio:
            continue

        results.append(dv)

    results.sort(key=lambda x: (-x.recurrence, -x.delta_alt_ratio))

    snp_count = sum(1 for r in results if r.variant_type == VariantType.SNP)
    indel_count = sum(1 for r in results if r.variant_type in
                      (VariantType.INSERTION, VariantType.DELETION))
    print(f"  Somatic variants: {len(results)} "
          f"(SNPs: {snp_count}, INDELs: {indel_count})", file=sys.stderr)

    return results


def to_variant_records(diff_variants: List[DifferentialVariant]) -> List[VariantRecord]:
    """DifferentialVariant → VariantRecord 변환 (VEP 출력용)"""
    records = []
    for dv in diff_variants:
        r = VariantRecord(
            chrom=dv.chrom, pos=dv.pos, ref=dv.ref, alt=dv.alt,
            qual=0.0, filter_status="PASS",
            alt_ratio=dv.target_mean_alt_ratio,
            total_depth=max(dv.target_depths) if dv.target_depths else 0,
            variant_type=dv.variant_type,
            indel_size=dv.indel_size,
            is_frameshift=dv.is_frameshift,
            near_indel=dv.near_indel,
        )
        records.append(r)
    return records


def write_results_csv(diff_variants: List[DifferentialVariant],
                      output_path: str):
    """결과를 CSV로 출력"""
    if not diff_variants:
        print(f"[WARNING] No variants to write to {output_path}", file=sys.stderr)
        return

    fieldnames = [
        'chrom', 'pos', 'ref', 'alt', 'variant_type',
        'indel_size', 'is_frameshift', 'near_indel',
        'recurrence', 'target_mean_alt_ratio', 'ref_mean_alt_ratio',
        'delta_alt_ratio', 'fisher_pvalue',
        'target_samples_found', 'ref_samples_found',
        'target_alt_ratios', 'target_depths',
        'ref_alt_ratios', 'ref_depths',
    ]

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for dv in diff_variants:
            writer.writerow({
                'chrom': dv.chrom,
                'pos': dv.pos,
                'ref': dv.ref,
                'alt': dv.alt,
                'variant_type': dv.variant_type.value,
                'indel_size': dv.indel_size,
                'is_frameshift': dv.is_frameshift,
                'near_indel': dv.near_indel,
                'recurrence': dv.recurrence,
                'target_mean_alt_ratio': round(dv.target_mean_alt_ratio, 4),
                'ref_mean_alt_ratio': round(dv.ref_mean_alt_ratio, 4),
                'delta_alt_ratio': round(dv.delta_alt_ratio, 4),
                'fisher_pvalue': f"{dv.fisher_pvalue:.6e}",
                'target_samples_found': dv.target_samples_found,
                'ref_samples_found': dv.ref_samples_found,
                'target_alt_ratios': ';'.join(f"{r:.3f}" for r in dv.target_alt_ratios),
                'target_depths': ';'.join(str(d) for d in dv.target_depths),
                'ref_alt_ratios': ';'.join(f"{r:.3f}" for r in dv.ref_alt_ratios),
                'ref_depths': ';'.join(str(d) for d in dv.ref_depths),
            })

    print(f"  Written {len(diff_variants)} variants to {output_path}",
          file=sys.stderr)
