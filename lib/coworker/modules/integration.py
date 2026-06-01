"""
integration.py — Variant + Differential Expression 통합 분석

Variant 결과와 DESeq2 결과를 통합하여 Tier별 candidate gene list를 생성합니다.

Tier 구조:
  Tier 1: High-impact variant + DE gene (최우선 후보)
  Tier 2: Any variant + DE gene
  Tier 3: DE gene only (variant 미발견 — regulatory 또는 coverage 부족 가능)
  Tier 4: Variant only (발현 변화 없음 — gain-of-function 가능)
"""

import csv
import math
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Set
from pathlib import Path

from modules.variant_analyzer import DifferentialVariant
from modules.vcf_parser import VariantType


@dataclass
class DEGene:
    """DESeq2 결과에서 읽어온 DE gene 정보"""
    gene_id: str
    gene_symbol: str = ""
    log2fc: float = 0.0
    padj: float = 1.0
    base_mean: float = 0.0


@dataclass
class VEPAnnotation:
    """VEP 결과에서 읽어온 variant annotation"""
    chrom: str = ""
    pos: int = 0
    allele: str = ""
    gene_id: str = ""
    gene_symbol: str = ""
    consequence: str = ""
    impact: str = ""  # HIGH, MODERATE, LOW, MODIFIER
    amino_acids: str = ""
    protein_position: str = ""
    sift: str = ""
    polyphen: str = ""
    extra: str = ""


@dataclass
class CandidateGene:
    """최종 candidate gene — variant + expression 통합"""
    gene_symbol: str
    gene_id: str = ""
    tier: int = 4
    # Variant 정보
    variants: List[Dict] = field(default_factory=list)
    n_variants: int = 0
    has_high_impact: bool = False
    has_frameshift: bool = False
    max_recurrence: int = 0
    # Expression 정보
    log2fc: float = 0.0
    padj: float = 1.0
    is_de: bool = False
    # Annotation
    consequences: List[str] = field(default_factory=list)
    sift_scores: List[str] = field(default_factory=list)
    polyphen_scores: List[str] = field(default_factory=list)

    @property
    def priority_score(self) -> float:
        """높을수록 우선순위 높음"""
        score = 0.0
        # Tier 기반
        score += (5 - self.tier) * 100
        # High impact bonus
        if self.has_high_impact:
            score += 50
        if self.has_frameshift:
            score += 30
        # Recurrence bonus
        score += self.max_recurrence * 10
        # DE significance bonus
        if self.is_de and self.padj > 0:
            score += min(50, -10 * math.log10(max(self.padj, 1e-300)))
        # Fold change bonus
        score += min(20, abs(self.log2fc) * 5)
        return score


# === High-impact consequence 정의 ===
HIGH_IMPACT_CONSEQUENCES = {
    'transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant',
    'stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost',
}

MODERATE_IMPACT_CONSEQUENCES = {
    'missense_variant', 'inframe_insertion', 'inframe_deletion',
    'protein_altering_variant', 'coding_sequence_variant',
}


def load_de_genes(filepath: str,
                  padj_cutoff: float = 0.05,
                  log2fc_cutoff: float = 1.0) -> Dict[str, DEGene]:
    """
    DESeq2 결과 CSV 읽기

    예상 컬럼: gene_id (or rownames), baseMean, log2FoldChange, padj
    """
    de_genes = {}

    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames

        # 컬럼명 자동 감지
        gene_col = None
        for candidate in ['gene_id', 'gene', 'Gene', 'symbol', 'gene_symbol', '']:
            if candidate in fieldnames:
                gene_col = candidate
                break
        # 첫 번째 unnamed 컬럼이 gene ID인 경우 (DESeq2 기본 출력)
        if gene_col == '' or gene_col is None:
            gene_col = fieldnames[0]

        l2fc_col = None
        for candidate in ['log2FoldChange', 'log2fc', 'logFC', 'lfc']:
            if candidate in fieldnames:
                l2fc_col = candidate
                break

        padj_col = None
        for candidate in ['padj', 'FDR', 'adj.P.Val', 'q_value']:
            if candidate in fieldnames:
                padj_col = candidate
                break

        basemean_col = None
        for candidate in ['baseMean', 'AveExpr', 'base_mean']:
            if candidate in fieldnames:
                basemean_col = candidate
                break

        if not l2fc_col or not padj_col:
            print(f"[ERROR] Cannot find log2FC or padj columns in {filepath}",
                  file=sys.stderr)
            print(f"  Available columns: {fieldnames}", file=sys.stderr)
            return de_genes

        for row in reader:
            gene_id = row[gene_col].strip()
            if not gene_id or gene_id == 'NA':
                continue

            try:
                l2fc = float(row[l2fc_col]) if row[l2fc_col] != 'NA' else 0.0
                padj = float(row[padj_col]) if row[padj_col] != 'NA' else 1.0
                basemean = (float(row[basemean_col])
                            if basemean_col and row.get(basemean_col, 'NA') != 'NA'
                            else 0.0)
            except (ValueError, KeyError):
                continue

            de_genes[gene_id] = DEGene(
                gene_id=gene_id,
                log2fc=l2fc,
                padj=padj,
                base_mean=basemean,
            )

    # padj + log2fc cutoff 적용
    sig_count = sum(1 for g in de_genes.values()
                    if g.padj < padj_cutoff and abs(g.log2fc) >= log2fc_cutoff)
    print(f"[INFO] Loaded {len(de_genes)} genes from DE results, "
          f"{sig_count} significant (padj<{padj_cutoff}, |log2FC|>={log2fc_cutoff})",
          file=sys.stderr)

    return de_genes


def load_vep_results(filepath: str) -> List[VEPAnnotation]:
    """VEP 결과 파일 파싱"""
    annotations = []

    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'):
                # Header line에서 컬럼명 추출
                if line.startswith('#Uploaded_variation') or line.startswith('#Up'):
                    header = line.strip('#').strip().split('\t')
                continue

            fields = line.strip().split('\t')
            if len(fields) < 7:
                continue

            ann = VEPAnnotation()

            # VEP 기본 출력 포맷
            # #Uploaded_variation  Location  Allele  Gene  Feature  Feature_type
            # Consequence  cDNA_position  CDS_position  Protein_position
            # Amino_acids  Codons  Existing_variation  Extra

            if len(fields) >= 1:
                # Location에서 chrom:pos 추출
                loc = fields[1] if len(fields) > 1 else ""
                if ':' in loc:
                    parts = loc.split(':')
                    ann.chrom = parts[0]
                    try:
                        pos_str = parts[1].split('-')[0]
                        ann.pos = int(pos_str)
                    except (ValueError, IndexError):
                        pass

            if len(fields) >= 3:
                ann.allele = fields[2]
            if len(fields) >= 4:
                ann.gene_id = fields[3]
            if len(fields) >= 7:
                ann.consequence = fields[6]
            if len(fields) >= 10:
                ann.protein_position = fields[9]
            if len(fields) >= 11:
                ann.amino_acids = fields[10]
            if len(fields) >= 14:
                ann.extra = fields[13]

            # Extra에서 SYMBOL 추출
            if ann.extra:
                for part in ann.extra.split(';'):
                    if part.startswith('SYMBOL='):
                        ann.gene_symbol = part.split('=')[1]
                    elif part.startswith('IMPACT='):
                        ann.impact = part.split('=')[1]
                    elif part.startswith('SIFT='):
                        ann.sift = part.split('=')[1]
                    elif part.startswith('PolyPhen='):
                        ann.polyphen = part.split('=')[1]

            annotations.append(ann)

    print(f"[INFO] Loaded {len(annotations)} VEP annotations from {filepath}",
          file=sys.stderr)
    return annotations


def integrate(diff_variants: List[DifferentialVariant],
              vep_annotations: List[VEPAnnotation],
              de_genes: Dict[str, DEGene],
              padj_cutoff: float = 0.05,
              log2fc_cutoff: float = 1.0
              ) -> List[CandidateGene]:
    """
    Variant + VEP annotation + DE gene 결과를 통합하여 candidate gene list 생성

    Returns: Tier별로 정렬된 CandidateGene 리스트
    """
    print(f"\n[Integration] {len(diff_variants)} variants, "
          f"{len(vep_annotations)} annotations, "
          f"{len(de_genes)} DE genes", file=sys.stderr)

    # Step 1: VEP annotation을 position으로 인덱싱
    vep_by_pos = defaultdict(list)
    for ann in vep_annotations:
        key = f"{ann.chrom}_{ann.pos}"
        vep_by_pos[key].append(ann)

    # Step 2: Variant → Gene 매핑
    gene_variants = defaultdict(list)  # gene_symbol → [variant_info]

    for dv in diff_variants:
        pos_key = f"{dv.chrom}_{dv.pos}"
        matched_anns = vep_by_pos.get(pos_key, [])

        if matched_anns:
            for ann in matched_anns:
                gene = ann.gene_symbol if ann.gene_symbol else ann.gene_id
                if not gene:
                    continue
                gene_variants[gene].append({
                    'chrom': dv.chrom,
                    'pos': dv.pos,
                    'ref': dv.ref,
                    'alt': dv.alt,
                    'variant_type': dv.variant_type.value,
                    'indel_size': dv.indel_size,
                    'is_frameshift': dv.is_frameshift,
                    'consequence': ann.consequence,
                    'impact': ann.impact,
                    'amino_acids': ann.amino_acids,
                    'sift': ann.sift,
                    'polyphen': ann.polyphen,
                    'recurrence': dv.recurrence,
                    'target_mean_alt_ratio': dv.target_mean_alt_ratio,
                    'fisher_pvalue': dv.fisher_pvalue,
                    'gene_id': ann.gene_id,
                })
        else:
            # VEP annotation 없음 (intergenic 등)
            gene_variants["_unannotated_"].append({
                'chrom': dv.chrom,
                'pos': dv.pos,
                'ref': dv.ref,
                'alt': dv.alt,
                'variant_type': dv.variant_type.value,
                'indel_size': dv.indel_size,
                'is_frameshift': dv.is_frameshift,
                'consequence': 'unknown',
                'impact': 'unknown',
                'recurrence': dv.recurrence,
                'target_mean_alt_ratio': dv.target_mean_alt_ratio,
            })

    # Step 3: Significant DE genes 세트
    sig_de = {gene_id for gene_id, g in de_genes.items()
              if g.padj < padj_cutoff and abs(g.log2fc) >= log2fc_cutoff}

    # DE gene을 gene_symbol로도 매칭 가능하게 확장
    all_de_symbols = set()
    de_by_symbol = {}
    for gene_id, g in de_genes.items():
        all_de_symbols.add(gene_id)
        de_by_symbol[gene_id] = g
        if g.gene_symbol:
            all_de_symbols.add(g.gene_symbol)
            de_by_symbol[g.gene_symbol] = g

    # Step 4: Candidate gene 생성 및 Tier 할당
    candidates = []
    genes_with_variants = set(gene_variants.keys()) - {"_unannotated_"}

    for gene, var_list in gene_variants.items():
        if gene == "_unannotated_":
            continue

        cand = CandidateGene(gene_symbol=gene)

        # Variant 정보
        cand.variants = var_list
        cand.n_variants = len(var_list)
        cand.consequences = list(set(v.get('consequence', '') for v in var_list))
        cand.sift_scores = [v.get('sift', '') for v in var_list if v.get('sift')]
        cand.polyphen_scores = [v.get('polyphen', '') for v in var_list
                                if v.get('polyphen')]
        cand.max_recurrence = max(v.get('recurrence', 0) for v in var_list)

        # Gene ID (VEP에서)
        gene_ids = [v.get('gene_id', '') for v in var_list if v.get('gene_id')]
        if gene_ids:
            cand.gene_id = gene_ids[0]

        # High impact 확인
        all_consequences = set()
        for v in var_list:
            for c in v.get('consequence', '').split(','):
                all_consequences.add(c.strip())

        cand.has_high_impact = bool(all_consequences & HIGH_IMPACT_CONSEQUENCES)
        cand.has_frameshift = any(v.get('is_frameshift', False) for v in var_list)

        # DE 정보 매칭
        de_match = de_by_symbol.get(gene) or de_by_symbol.get(cand.gene_id)
        if de_match:
            cand.log2fc = de_match.log2fc
            cand.padj = de_match.padj
            cand.is_de = (de_match.padj < padj_cutoff
                          and abs(de_match.log2fc) >= log2fc_cutoff)

        # Tier 할당
        if cand.has_high_impact and cand.is_de:
            cand.tier = 1  # High-impact variant + DE
        elif cand.is_de:  # Any variant + DE
            cand.tier = 2
        elif cand.has_high_impact:  # High-impact variant, no DE
            cand.tier = 3
        else:
            cand.tier = 4  # Low-impact variant only

        candidates.append(cand)

    # Step 5: DE only genes (variant 없음)
    for gene_id, g in de_genes.items():
        if (g.padj < padj_cutoff and abs(g.log2fc) >= log2fc_cutoff
                and gene_id not in genes_with_variants):
            cand = CandidateGene(
                gene_symbol=gene_id,
                gene_id=gene_id,
                tier=5,  # DE only
                log2fc=g.log2fc,
                padj=g.padj,
                is_de=True,
            )
            candidates.append(cand)

    # Sort by tier, then priority score
    candidates.sort(key=lambda c: (c.tier, -c.priority_score))

    # Summary
    tier_counts = defaultdict(int)
    for c in candidates:
        tier_counts[c.tier] += 1

    print(f"\n  === Candidate Gene Summary ===", file=sys.stderr)
    print(f"  Tier 1 (High-impact variant + DE): {tier_counts[1]}", file=sys.stderr)
    print(f"  Tier 2 (Any variant + DE):         {tier_counts[2]}", file=sys.stderr)
    print(f"  Tier 3 (High-impact, no DE):       {tier_counts[3]}", file=sys.stderr)
    print(f"  Tier 4 (Low-impact variant only):  {tier_counts[4]}", file=sys.stderr)
    print(f"  Tier 5 (DE only, no variant):      {tier_counts[5]}", file=sys.stderr)
    print(f"  Total candidates: {len(candidates)}", file=sys.stderr)

    return candidates


def write_candidate_genes(candidates: List[CandidateGene],
                          output_path: str,
                          max_tier: int = 5):
    """Candidate gene list를 CSV로 출력"""
    filtered = [c for c in candidates if c.tier <= max_tier]

    fieldnames = [
        'tier', 'gene_symbol', 'gene_id',
        'n_variants', 'consequences', 'has_high_impact', 'has_frameshift',
        'max_recurrence', 'sift', 'polyphen',
        'log2fc', 'padj', 'is_de',
        'priority_score',
        'variant_details',
    ]

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for c in filtered:
            # Variant details를 compact string으로
            var_details = []
            for v in c.variants:
                detail = (f"{v.get('chrom','')}:{v.get('pos','')} "
                          f"{v.get('ref','')}->{v.get('alt','')} "
                          f"({v.get('consequence','')}) "
                          f"[recurrence={v.get('recurrence',0)}]")
                var_details.append(detail)

            writer.writerow({
                'tier': c.tier,
                'gene_symbol': c.gene_symbol,
                'gene_id': c.gene_id,
                'n_variants': c.n_variants,
                'consequences': '|'.join(c.consequences),
                'has_high_impact': c.has_high_impact,
                'has_frameshift': c.has_frameshift,
                'max_recurrence': c.max_recurrence,
                'sift': '|'.join(c.sift_scores) if c.sift_scores else '',
                'polyphen': '|'.join(c.polyphen_scores) if c.polyphen_scores else '',
                'log2fc': round(c.log2fc, 4) if c.log2fc else '',
                'padj': f"{c.padj:.6e}" if c.padj < 1.0 else '',
                'is_de': c.is_de,
                'priority_score': round(c.priority_score, 1),
                'variant_details': ' | '.join(var_details),
            })

    print(f"\n[Output] Written {len(filtered)} candidate genes to {output_path}",
          file=sys.stderr)
