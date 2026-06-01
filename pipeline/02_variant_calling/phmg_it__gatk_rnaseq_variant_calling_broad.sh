#!/bin/bash
# ==============================================================================
# RNA-seq Variant Calling Pipeline — Broadened Detection
# ==============================================================================
# Modified from GATK Best Practices for RNA-seq
# Key changes: relaxed HaplotypeCaller & VariantFiltration for maximum sensitivity
# ==============================================================================

# Reference genome
REF="reference/mouse_genome_GRCm39.fna"

# Sample list
SAMPLES=("disease_5" "disease_6")

# Output directory
OUTDIR="gatk_output_broad"
mkdir -p $OUTDIR

# Number of threads
THREADS=8

# Loop through all samples
for SAMPLE in "${SAMPLES[@]}"; do
    echo "=========================================="
    echo "Processing sample: $SAMPLE"
    echo "=========================================="
    
    # 1. Add Read Groups
    echo "[1/6] Adding read groups for $SAMPLE..."
    gatk AddOrReplaceReadGroups \
         -I bam_dz/${SAMPLE}Aligned.sorted.out.bam \
         -O ${OUTDIR}/${SAMPLE}_rg.bam \
         -RGID ${SAMPLE} \
         -RGLB lib1 \
         -RGPL ILLUMINA \
         -RGPU unit1 \
         -RGSM ${SAMPLE} \
         --CREATE_INDEX true
    
    # 2. Mark Duplicates
    echo "[2/6] Marking duplicates for $SAMPLE..."
    gatk MarkDuplicates \
         -I ${OUTDIR}/${SAMPLE}_rg.bam \
         -O ${OUTDIR}/${SAMPLE}_marked.bam \
         -M ${OUTDIR}/${SAMPLE}_metrics.txt \
         --CREATE_INDEX true
    
    # =========================================================================
    # 3. Split N Cigar Reads
    # =========================================================================
    # --skip-mapping-quality-transform: STAR assigns MAPQ=255 to unique reads,
    #   which SplitNCigarReads normally converts to 60. Skipping this preserves
    #   original MAPQ values, allowing more reads through downstream filters.
    # =========================================================================
    echo "[3/6] Splitting N cigar reads for $SAMPLE..."
    gatk SplitNCigarReads \
         -R $REF \
         -I ${OUTDIR}/${SAMPLE}_marked.bam \
         -O ${OUTDIR}/${SAMPLE}_split.bam
    
    # =========================================================================
    # 4. Variant Calling — HaplotypeCaller (★ most impact)
    # =========================================================================
    #
    # Parameter changes and rationale:
    #
    # --standard-min-confidence-threshold-for-calling 10.0
    #   [WAS: 20.0]  Lower confidence threshold → low-expression gene variants
    #   included in raw VCF. VariantFiltration handles quality control downstream.
    #
    # --minimum-mapping-quality 10
    #   [WAS: 20 (default)]  Include reads with lower MAPQ → variants in
    #   paralog regions (Speer4, Gm-series) and multi-copy gene families.
    #
    # --min-base-quality-score 6
    #   [WAS: 10 (default)]  Lower base quality bases contribute as evidence →
    #   helps in low-coverage regions near splice junctions.
    #
    # --active-probability-threshold 0.001
    #   [WAS: 0.002 (default)]  More genomic regions scanned for variants →
    #   low-frequency variant regions not skipped.
    #
    # --max-reads-per-alignment-start 0
    #   [WAS: 50 (default)]  Disable downsampling → high-expression genes
    #   (Ig locus, Myc) fully utilized without read cap.
    #
    # --dont-use-soft-clipped-bases
    #   [KEPT]  Still recommended for RNA-seq to avoid splice junction artifacts.
    #
    # =========================================================================
    echo "[4/6] Calling variants for $SAMPLE..."
    gatk HaplotypeCaller \
         -R $REF \
         -I ${OUTDIR}/${SAMPLE}_split.bam \
         -O ${OUTDIR}/${SAMPLE}_raw.vcf \
         --dont-use-soft-clipped-bases \
         --standard-min-confidence-threshold-for-calling 10.0 \
         --minimum-mapping-quality 10 \
         --min-base-quality-score 6 \
         --active-probability-threshold 0.001 \
         --max-reads-per-alignment-start 0 \
         --native-pair-hmm-threads ${THREADS}
    
    # =========================================================================
    # 5. Variant Filtration
    # =========================================================================
    #
    # QD < 1.0    [WAS: 2.0]  RNA-seq has inherently lower QD in low-coverage
    #             genes. Relaxing prevents over-filtering expressed variants.
    #
    # FS > 60.0   [WAS: 30.0]  RNA-seq naturally exhibits strand bias because
    #             transcripts are single-stranded. WGS threshold (30) is too
    #             aggressive for RNA data.
    #
    # MQ < 30.0   [WAS: 40.0]  Allows variants in regions with moderate mapping
    #             quality, recovering paralog/repetitive region variants.
    #
    # =========================================================================
    echo "[5/6] Filtering variants for $SAMPLE..."
    gatk VariantFiltration \
         -R $REF \
         -V ${OUTDIR}/${SAMPLE}_raw.vcf \
         -O ${OUTDIR}/${SAMPLE}_filtered.vcf \
         --filter-expression "QD < 1.0"  --filter-name "QD1" \
         --filter-expression "FS > 60.0" --filter-name "FS60" \
         --filter-expression "MQ < 30.0" --filter-name "MQ30"
    
    # 6. Extract PASS only
    echo "[6/6] Extracting PASS variants for $SAMPLE..."
    gatk SelectVariants \
         -R $REF \
         -V ${OUTDIR}/${SAMPLE}_filtered.vcf \
         -O ${OUTDIR}/${SAMPLE}_final.vcf \
         --exclude-filtered
    
    # =========================================================================
    # Stats
    # =========================================================================
    RAW_COUNT=$(grep -c -v "^#" ${OUTDIR}/${SAMPLE}_raw.vcf)
    FILT_COUNT=$(grep -c -v "^#" ${OUTDIR}/${SAMPLE}_filtered.vcf)
    PASS_COUNT=$(grep -c -v "^#" ${OUTDIR}/${SAMPLE}_final.vcf)
    
    echo ""
    echo "  Results for $SAMPLE:"
    echo "    Raw variants:      $RAW_COUNT"
    echo "    After filtration:  $FILT_COUNT"
    echo "    PASS variants:     $PASS_COUNT"
    echo ""
    echo "Completed: $SAMPLE"
    echo ""
done

echo "=========================================="
echo "All samples processed successfully!"
echo "=========================================="
echo ""
echo "Parameter summary (vs original):"
echo "  HaplotypeCaller:"
echo "    confidence threshold:  20 → 10"
echo "    min mapping quality:   20 → 10"
echo "    min base quality:      10 → 6"
echo "    active probability:    0.002 → 0.001"
echo "    max reads/start:       50 → unlimited"
echo "  SplitNCigarReads:"
echo "    skip MAPQ transform:   no → yes"
echo "  VariantFiltration:"
echo "    QD:                    <2.0 → <1.0"
echo "    FS:                    >30 → >60"
echo "    MQ:                    <40 → <30"
echo "=========================================="
