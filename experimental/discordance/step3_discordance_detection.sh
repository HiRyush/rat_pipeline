#!/bin/bash
# ==============================================================================
# Step 3: Discordance Detection — Observed RNA-seq vs Imputed Germline
# ==============================================================================
# Per-sample: MuTect2 PASS variant (observed ALT) + Imputed 0|0 (expected REF)
# = Discordance → somatic candidate
# ==============================================================================
set -euo pipefail
source ~/RAT_project/miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

MUTECT_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/results/mutect2/force_filtered"
IMPUTED_DIR="/home/yusanghyeon/RAT_project/PHMG_IT/imputation/imputed"
OUT_DIR="/home/yusanghyeon/RAT_project/pipeline_candidates/results/discordance"
mkdir -p "$OUT_DIR"

SAMPLES="C1 C2 C3 C4 C5 P1 P2 P3 P4 P5 P6 P7 P8 P9 P10"

echo "================================================================"
echo "Step 3: Discordance Detection"
echo "  Observed (MuTect2 PASS ALT) vs Imputed (0|0 REF)"
echo "================================================================"

for S in $SAMPLES; do
    echo ""
    echo "--- Processing $S ---"

    # 1. Extract MuTect2 PASS variants with DP>=5 and AF>=0.05
    #    (observed ALT alleles in RNA-seq)
    echo "  Extracting PASS variants from MuTect2..."
    bcftools view -f PASS "$MUTECT_DIR/${S}.force_filtered.vcf.gz" | \
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%DP\t%AF]\n' | \
        awk -F'\t' '$5 >= 5 && $6 >= 0.05' | \
        sort -k1,1 -k2,2n > "$OUT_DIR/${S}_observed.tsv"

    OBS_COUNT=$(wc -l < "$OUT_DIR/${S}_observed.tsv")
    echo "  Observed PASS variants (DP>=5, AF>=0.05): $OBS_COUNT"

    # 2. Merge per-chromosome imputed VCFs and extract 0|0 (REF) sites
    #    These are positions where imputation predicts germline = REF
    echo "  Extracting imputed REF (0|0) sites..."

    # Build list of imputed positions with their genotype
    > "$OUT_DIR/${S}_imputed_ref.tsv"
    for CHR in $(seq 1 20); do
        IMP_VCF="$IMPUTED_DIR/${S}_chr${CHR}.vcf.gz"
        if [ -f "$IMP_VCF" ]; then
            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%DS]\n' "$IMP_VCF" | \
                awk -F'\t' '$5 == "0|0"' >> "$OUT_DIR/${S}_imputed_ref.tsv"
        fi
    done
    sort -k1,1 -k2,2n "$OUT_DIR/${S}_imputed_ref.tsv" -o "$OUT_DIR/${S}_imputed_ref.tsv"

    IMP_REF_COUNT=$(wc -l < "$OUT_DIR/${S}_imputed_ref.tsv")
    echo "  Imputed REF (0|0) sites: $IMP_REF_COUNT"

    # 3. Find discordances: positions where observed=ALT but imputed=REF(0|0)
    #    Join on chrom+pos+ref+alt
    echo "  Finding discordances (observed ALT + imputed REF)..."

    # Match by chrom:pos — observed ALT at site where imputation says 0|0
    awk -F'\t' '{print $1":"$2}' "$OUT_DIR/${S}_observed.tsv" | sort > "$OUT_DIR/${S}_obs_keys.tmp"
    awk -F'\t' '{print $1":"$2}' "$OUT_DIR/${S}_imputed_ref.tsv" | sort > "$OUT_DIR/${S}_imp_keys.tmp"

    # Intersection: sites that are BOTH observed-ALT AND imputed-REF
    comm -12 "$OUT_DIR/${S}_obs_keys.tmp" "$OUT_DIR/${S}_imp_keys.tmp" > "$OUT_DIR/${S}_discord_keys.tmp"

    DISCORD_COUNT=$(wc -l < "$OUT_DIR/${S}_discord_keys.tmp")
    echo "  Discordances found: $DISCORD_COUNT"

    # 4. Build discordance output with full info
    #    chrom, pos, ref, alt, obs_DP, obs_AF, imp_GT, imp_DS
    if [ "$DISCORD_COUNT" -gt 0 ]; then
        # Create indexed lookup for observed data
        awk -F'\t' '{print $1":"$2"\t"$0}' "$OUT_DIR/${S}_observed.tsv" | sort -k1,1 > "$OUT_DIR/${S}_obs_indexed.tmp"
        awk -F'\t' '{print $1":"$2"\t"$5"\t"$6}' "$OUT_DIR/${S}_imputed_ref.tsv" | sort -k1,1 > "$OUT_DIR/${S}_imp_indexed.tmp"

        # Join on key
        join -t$'\t' -1 1 -2 1 "$OUT_DIR/${S}_obs_indexed.tmp" "$OUT_DIR/${S}_imp_indexed.tmp" | \
            awk -F'\t' 'BEGIN{OFS="\t"} {
                # key, chrom, pos, ref, alt, obs_dp, obs_af, imp_gt, imp_ds
                split($1, a, ":");
                print a[1], a[2], $4, $5, $6, $7, $8, $9
            }' | grep -v "^$" | sort -k1,1 -k2,2n > "$OUT_DIR/${S}_discordances.tsv"

        FINAL_COUNT=$(wc -l < "$OUT_DIR/${S}_discordances.tsv")
        echo "  Final discordance output: $FINAL_COUNT"
    else
        > "$OUT_DIR/${S}_discordances.tsv"
        echo "  Final discordance output: 0"
    fi

    # Cleanup temp files
    rm -f "$OUT_DIR/${S}_obs_keys.tmp" "$OUT_DIR/${S}_imp_keys.tmp" \
          "$OUT_DIR/${S}_discord_keys.tmp" "$OUT_DIR/${S}_obs_indexed.tmp" \
          "$OUT_DIR/${S}_imp_indexed.tmp"

    echo "  Done: $S"
done

# 5. Summary
echo ""
echo "================================================================"
echo "Summary: Discordance counts per sample"
echo "================================================================"
echo -e "Sample\tObserved\tImputed_REF\tDiscordances" > "$OUT_DIR/discordance_summary.tsv"
for S in $SAMPLES; do
    OBS=$(wc -l < "$OUT_DIR/${S}_observed.tsv")
    IMP=$(wc -l < "$OUT_DIR/${S}_imputed_ref.tsv")
    DIS=$(wc -l < "$OUT_DIR/${S}_discordances.tsv")
    echo -e "${S}\t${OBS}\t${IMP}\t${DIS}" >> "$OUT_DIR/discordance_summary.tsv"
    echo "  $S: Observed=$OBS, Imputed_REF=$IMP, Discordances=$DIS"
done

echo ""
echo "Results saved to: $OUT_DIR"
echo "Done."
