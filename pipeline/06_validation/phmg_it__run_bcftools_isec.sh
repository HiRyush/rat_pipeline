#!/bin/bash
# ============================================================
# run_bcftools_isec.sh
#
# bcftools isec을 이용한 condition-specific variant 추출
# variant_analyzer 결과와 비교하기 위한 gold standard 생성
#
# Compatible with bcftools >= 1.3
#
# 사용법:
#   bash run_bcftools_isec.sh \
#       --disease disease_1.vcf.gz disease_2.vcf.gz \
#       --control control_1.vcf.gz control_2.vcf.gz \
#       --output isec_results
#
# 이미 merge된 VCF가 있으면:
#   bash run_bcftools_isec.sh \
#       --disease-merged disease_merged.vcf.gz \
#       --control-merged control_merged.vcf.gz \
#       --output isec_results
# ============================================================

set -euo pipefail

# Default values
OUTPUT_DIR="isec_results"
DISEASE_FILES=()
CONTROL_FILES=()
DISEASE_MERGED=""
CONTROL_MERGED=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --disease)
            shift
            while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
                DISEASE_FILES+=("$1")
                shift
            done
            ;;
        --control)
            shift
            while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
                CONTROL_FILES+=("$1")
                shift
            done
            ;;
        --disease-merged)
            DISEASE_MERGED="$2"
            shift 2
            ;;
        --control-merged)
            CONTROL_MERGED="$2"
            shift 2
            ;;
        --output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

mkdir -p "${OUTPUT_DIR}"

echo "============================================================"
echo "  bcftools isec — Gold Standard Generation"
echo "============================================================"

# Check dependencies
if ! command -v bcftools &> /dev/null; then
    echo "ERROR: bcftools not found."
    exit 1
fi
if ! command -v tabix &> /dev/null; then
    echo "ERROR: tabix not found (part of htslib)."
    exit 1
fi

BCFTOOLS_VER=$(bcftools 2>&1 | head -1 || true)
echo "  ${BCFTOOLS_VER}"
echo ""


# ============================================================
# Helper: ensure VCF is bgzipped + tabix indexed
# ============================================================
ensure_gz_indexed() {
    local input_vcf="$1"
    local output_gz="$2"   # used only if input is plain VCF

    if [[ "${input_vcf}" == *.vcf.gz || "${input_vcf}" == *.bgz ]]; then
        # Already gzipped — ensure tabix index exists
        if [[ ! -f "${input_vcf}.tbi" ]]; then
            echo "    tabix: ${input_vcf}"
            tabix -p vcf "${input_vcf}"
        fi
        echo "${input_vcf}"
    else
        # Plain VCF → bgzip copy → tabix
        echo "    bgzip+tabix: ${input_vcf}"
        bcftools view "${input_vcf}" -Oz -o "${output_gz}"
        tabix -p vcf "${output_gz}"
        echo "${output_gz}"
    fi
}


# ============================================================
# Helper: rename sample in VCF (avoids merge duplicate error)
# bcftools reheader -s works in bcftools >= 1.3
# ============================================================
rename_and_prepare() {
    local input_vcf="$1"
    local new_sample_name="$2"
    local output_gz="$3"

    local sample_file="${output_gz%.vcf.gz}.sname.txt"
    echo "${new_sample_name}" > "${sample_file}"

    # reheader -o preserves input format (gz in → gz out)
    bcftools reheader -s "${sample_file}" -o "${output_gz}" "${input_vcf}"
    tabix -p vcf "${output_gz}"

    rm -f "${sample_file}"
}


# ============================================================
# Helper: merge multiple VCFs
# Uses bcftools merge — options BEFORE input files
# ============================================================
merge_vcfs() {
    local output_path="$1"
    shift
    local vcf_files=("$@")

    if [[ ${#vcf_files[@]} -eq 1 ]]; then
        cp "${vcf_files[0]}" "${output_path}"
        cp "${vcf_files[0]}.tbi" "${output_path}.tbi" 2>/dev/null || \
            tabix -p vcf "${output_path}"
        return
    fi

    echo "  Merging ${#vcf_files[@]} VCF files..."

    # bcftools merge: all options before positional args
    bcftools merge \
        -Oz -o "${output_path}" \
        "${vcf_files[@]}"

    tabix -p vcf "${output_path}"
}


# ============================================================
# Main: Prepare and merge VCFs
# ============================================================

if [[ -n "${DISEASE_MERGED}" && -n "${CONTROL_MERGED}" ]]; then
    echo "Using pre-merged VCFs"
    DISEASE_GZ=$(ensure_gz_indexed "${DISEASE_MERGED}" "${OUTPUT_DIR}/disease_merged.vcf.gz")
    CONTROL_GZ=$(ensure_gz_indexed "${CONTROL_MERGED}" "${OUTPUT_DIR}/control_merged.vcf.gz")

elif [[ ${#DISEASE_FILES[@]} -gt 0 && ${#CONTROL_FILES[@]} -gt 0 ]]; then
    echo "Step 1: Preparing VCF files (rename samples to avoid conflicts)"
    echo "  Disease: ${#DISEASE_FILES[@]} files"
    echo "  Control: ${#CONTROL_FILES[@]} files"
    echo ""

    DISEASE_PREPARED=()
    for i in "${!DISEASE_FILES[@]}"; do
        vcf="${DISEASE_FILES[$i]}"
        prepared="${OUTPUT_DIR}/disease_s${i}.vcf.gz"
        echo "  [disease_${i}] ${vcf}"
        rename_and_prepare "${vcf}" "disease_${i}" "${prepared}"
        DISEASE_PREPARED+=("${prepared}")
    done

    CONTROL_PREPARED=()
    for i in "${!CONTROL_FILES[@]}"; do
        vcf="${CONTROL_FILES[$i]}"
        prepared="${OUTPUT_DIR}/control_s${i}.vcf.gz"
        echo "  [control_${i}] ${vcf}"
        rename_and_prepare "${vcf}" "control_${i}" "${prepared}"
        CONTROL_PREPARED+=("${prepared}")
    done

    echo ""
    echo "Step 2: Merging VCF files"

    DISEASE_GZ="${OUTPUT_DIR}/disease_merged.vcf.gz"
    CONTROL_GZ="${OUTPUT_DIR}/control_merged.vcf.gz"

    merge_vcfs "${DISEASE_GZ}" "${DISEASE_PREPARED[@]}"
    merge_vcfs "${CONTROL_GZ}" "${CONTROL_PREPARED[@]}"
else
    echo "ERROR: Provide either --disease/--control or --disease-merged/--control-merged"
    exit 1
fi

echo ""
echo "  Disease VCF: ${DISEASE_GZ}"
echo "  Control VCF: ${CONTROL_GZ}"


# ============================================================
# Run bcftools isec
# ============================================================
echo ""
echo "Step 3: Running bcftools isec"

ISEC_DIR="${OUTPUT_DIR}/isec"
bcftools isec \
    -p "${ISEC_DIR}" \
    "${DISEASE_GZ}" \
    "${CONTROL_GZ}"


# ============================================================
# Results summary
# ============================================================
echo ""
echo "============================================================"
echo "  bcftools isec Results"
echo "============================================================"

# Auto-detect output format
if [[ -f "${ISEC_DIR}/0000.vcf" ]]; then
    ISEC_EXT="vcf"
elif [[ -f "${ISEC_DIR}/0000.vcf.gz" ]]; then
    ISEC_EXT="vcf.gz"
else
    echo "ERROR: No isec output files found in ${ISEC_DIR}"
    exit 1
fi

count_variants() {
    local f="$1"
    if [[ "${f}" == *.gz ]]; then
        zgrep -c -v "^#" "${f}" 2>/dev/null || echo "0"
    else
        grep -c -v "^#" "${f}" 2>/dev/null || echo "0"
    fi
}

DISEASE_ONLY="${ISEC_DIR}/0000.${ISEC_EXT}"
CONTROL_ONLY="${ISEC_DIR}/0001.${ISEC_EXT}"
SHARED="${ISEC_DIR}/0002.${ISEC_EXT}"

N_DISEASE_ONLY=$(count_variants "${DISEASE_ONLY}")
N_CONTROL_ONLY=$(count_variants "${CONTROL_ONLY}")
N_SHARED=$(count_variants "${SHARED}")

echo ""
echo "  Disease-only (0000.${ISEC_EXT}):  ${N_DISEASE_ONLY}"
echo "  Control-only (0001.${ISEC_EXT}):  ${N_CONTROL_ONLY}"
echo "  Shared       (0002.${ISEC_EXT}):  ${N_SHARED}"
echo ""
echo "  Gold standard: ${DISEASE_ONLY}"
echo ""
echo "============================================================"
echo "  Next: Run benchmark"
echo "============================================================"
echo ""
echo "  python benchmark.py \\"
echo "      --analyzer-csv differential_all_variants.csv \\"
echo "      --disease-only-vcf ${DISEASE_ONLY} \\"
echo "      --output benchmark_results/"
echo ""
