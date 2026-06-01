"""
vcf_parser.py — Unified VCF parser for bcftools and GATK HaplotypeCaller (macOS/GZ supported)
"""

import re
import csv
import sys
import gzip
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple
from pathlib import Path
from enum import Enum


class VCFSource(Enum):
    BCFTOOLS = "bcftools"
    GATK = "gatk"
    UNKNOWN = "unknown"


class VariantType(Enum):
    SNP = "SNP"
    INSERTION = "INSERTION"
    DELETION = "DELETION"
    MNP = "MNP"
    COMPLEX = "COMPLEX"


@dataclass
class VariantRecord:
    """통합된 variant 데이터 구조 — bcftools/GATK 공통"""
    chrom: str
    pos: int
    ref: str
    alt: str
    qual: float
    filter_status: str
    ref_count: int = 0
    alt_count: int = 0
    total_depth: int = 0
    alt_ratio: float = 0.0
    fwd_ref: int = -1
    rev_ref: int = -1
    fwd_alt: int = -1
    rev_alt: int = -1
    genotype: str = ""
    genotype_quality: int = -1
    variant_type: VariantType = VariantType.SNP
    indel_size: int = 0
    is_frameshift: bool = False
    source: VCFSource = VCFSource.UNKNOWN
    near_indel: bool = False
    raw_info: str = ""

    @property
    def is_indel(self) -> bool:
        return self.variant_type in (VariantType.INSERTION, VariantType.DELETION)

    @property
    def is_snp(self) -> bool:
        return self.variant_type == VariantType.SNP

    @property
    def position_key(self) -> str:
        """Fixes the AttributeError in variant_analyzer.py"""
        return f"{self.chrom}_{self.pos}_{self.alt}"

    @property
    def locus_key(self) -> str:
        return f"{self.chrom}_{self.pos}"


def detect_vcf_source(filepath: str) -> VCFSource:
    """Detects VCF source with GZIP support"""
    is_gz = str(filepath).endswith('.gz')
    opener = gzip.open if is_gz else open
    mode = 'rt' if is_gz else 'r'
    
    try:
        with opener(filepath, mode) as f:
            for line in f:
                if not line.startswith('#'):
                    break
                line_lower = line.lower()
                if 'bcftools' in line_lower or 'samtools' in line_lower:
                    return VCFSource.BCFTOOLS
                if 'haplotypecaller' in line_lower or 'gatk' in line_lower:
                    return VCFSource.GATK
                if line.startswith('##FORMAT=<ID=AD'):
                    return VCFSource.GATK
                if '##INFO=<ID=DP4' in line or '##INFO=<ID=I16' in line:
                    return VCFSource.BCFTOOLS
    except Exception:
        pass
    return VCFSource.UNKNOWN


def classify_variant(ref: str, alt: str) -> Tuple[VariantType, int, bool]:
    ref_len, alt_len = len(ref), len(alt)
    if ref_len == 1 and alt_len == 1:
        return VariantType.SNP, 0, False
    elif alt_len > ref_len:
        size = alt_len - ref_len
        return VariantType.INSERTION, size, (size % 3 != 0)
    elif ref_len > alt_len:
        size = ref_len - alt_len
        return VariantType.DELETION, size, (size % 3 != 0)
    return VariantType.MNP, 0, False


def _parse_bcftools_allele_counts(info: str) -> Dict:
    result = {'fwd_ref': 0, 'rev_ref': 0, 'fwd_alt': 0, 'rev_alt': 0, 'dp': 0}
    dp4_match = re.search(r'DP4=(\d+),(\d+),(\d+),(\d+)', info)
    if dp4_match:
        result['fwd_ref'] = int(dp4_match.group(1))
        result['rev_ref'] = int(dp4_match.group(2))
        result['fwd_alt'] = int(dp4_match.group(3))
        result['rev_alt'] = int(dp4_match.group(4))
        result['dp'] = sum(result[k] for k in ['fwd_ref', 'rev_ref', 'fwd_alt', 'rev_alt'])
    return result


def _parse_gatk_allele_counts(format_str: str, sample_str: str) -> Dict:
    result = {'ref_count': 0, 'alt_count': 0, 'dp': 0, 'gt': './.', 'gq': 0}
    fmt_fields = format_str.split(':')
    sample_fields = sample_str.split(':')
    field_map = {k: v for k, v in zip(fmt_fields, sample_fields)}
    if 'GT' in field_map: result['gt'] = field_map['GT']
    if 'AD' in field_map and ',' in field_map['AD']:
        parts = field_map['AD'].split(',')
        result['ref_count'], result['alt_count'] = int(parts[0]), int(parts[1])
    if 'DP' in field_map:
        try: result['dp'] = int(field_map['DP'])
        except: pass
    if result['dp'] == 0: result['dp'] = result['ref_count'] + result['alt_count']
    return result


def parse_vcf(filepath: str, chromosomes: Optional[List[str]] = None, min_qual: float = 30.0) -> List[VariantRecord]:
    source = detect_vcf_source(filepath)
    if source == VCFSource.UNKNOWN:
        source = VCFSource.BCFTOOLS

    is_gz = str(filepath).endswith('.gz')
    opener = gzip.open if is_gz else open
    mode = 'rt' if is_gz else 'r'

    records = []
    with opener(filepath, mode) as f:
        for line in f:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len(fields) < 8: continue
            
            chrom, pos, _, ref, alt_field, qual_val, filter_status, info = fields[:8]
            try:
                qual = float(qual_val) if qual_val != '.' else 0.0
            except ValueError: qual = 0.0
            
            if qual < min_qual and qual_val != '.': continue

            for alt in alt_field.split(','):
                vtype, indel_size, is_fs = classify_variant(ref, alt)
                record = VariantRecord(
                    chrom=chrom, pos=int(pos), ref=ref, alt=alt,
                    qual=qual, filter_status=filter_status,
                    variant_type=vtype, indel_size=indel_size,
                    is_frameshift=is_fs, source=source, raw_info=info
                )

                if source == VCFSource.BCFTOOLS:
                    counts = _parse_bcftools_allele_counts(info)
                    record.fwd_ref, record.rev_ref = counts['fwd_ref'], counts['rev_ref']
                    record.fwd_alt, record.rev_alt = counts['fwd_alt'], counts['rev_alt']
                    record.ref_count = record.fwd_ref + record.rev_ref
                    record.alt_count = record.fwd_alt + record.rev_alt
                    record.total_depth = counts['dp']
                elif source == VCFSource.GATK and len(fields) >= 10:
                    counts = _parse_gatk_allele_counts(fields[8], fields[9])
                    record.ref_count, record.alt_count = counts['ref_count'], counts['alt_count']
                    record.total_depth, record.genotype = counts['dp'], counts['gt']

                if record.total_depth > 0:
                    record.alt_ratio = record.alt_count / record.total_depth
                records.append(record)
    return records


def flag_near_indel(records: List[VariantRecord], distance: int = 10) -> List[VariantRecord]:
    by_chrom = {}
    for r in records:
        by_chrom.setdefault(r.chrom, []).append(r)
    for chrom, chrom_records in by_chrom.items():
        indel_positions = [r.pos for r in chrom_records if r.is_indel]
        for r in chrom_records:
            if r.is_snp:
                for ip in indel_positions:
                    if abs(r.pos - ip) <= distance:
                        r.near_indel = True
                        break
    return records


def split_snp_indel(records: List[VariantRecord]) -> Tuple[List[VariantRecord], List[VariantRecord]]:
    return [r for r in records if r.is_snp], [r for r in records if r.is_indel]


def filter_variants(records: List[VariantRecord], min_coverage: int = 10, min_alt_ratio: float = 0.0, max_alt_ratio: float = 1.0, require_both_strands: bool = False) -> List[VariantRecord]:
    filtered = []
    for r in records:
        if r.total_depth < min_coverage: continue
        if r.alt_ratio < min_alt_ratio or r.alt_ratio > max_alt_ratio: continue
        filtered.append(r)
    return filtered


def records_to_table(records: List[VariantRecord]) -> List[Dict]:
    return [{'chrom': r.chrom, 'pos': r.pos, 'ref': r.ref, 'alt': r.alt, 'variant_type': r.variant_type.value, 'alt_ratio': round(r.alt_ratio, 4), 'total_depth': r.total_depth} for r in records]


def write_vcf_for_vep(records: List[VariantRecord], output_path: str):
    with open(output_path, 'w') as f:
        f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for r in records:
            f.write(f"{r.chrom}\t{r.pos}\t.\t{r.ref}\t{r.alt}\t{r.qual:.1f}\tPASS\tDP={r.total_depth}\n")