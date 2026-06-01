#!/bin/bash
set -euo pipefail

# HRDP Reference Panel 다운로드 (rn7/mRatBN7.2)
# 고유 strain만 선별 (중복 연도/NG 제거, 최신 우선)

BASE_URL="https://download.rgd.mcw.edu/strain_specific_variants/Dwinell_MCW_HybridRatDiversityProgram/Feb2024/mRatBN7/JointAnalysis_split_by_SingleSample/SNPs"
OUT="/home/yusanghyeon/RAT_project/PHMG_IT/imputation/ref_panel/hrdp_vcf"
mkdir -p "${OUT}"

# 선별 strain 목록 (92개 중 BN 계열 제거 → 다양한 ~80 strain)
# BN은 reference strain이라 variant가 거의 없어 imputation에 도움 안 됨
STRAINS=(
    ACI_N_2020
    BDIX_NemOdaMcwi_2022
    BUF_MnaMcwi_2022
    BXH2_CubMcwi_2023
    BXH3_CubMcwi_2020
    BXH6_CubMcwi_2021
    DA_OlaHsd_2019
    F344_DuCrl_2021
    F344_NCrl_2019
    F344_NHsd_2021
    F344_StmMcwi_2019
    FHH_EurMcwi_2019
    FXLE12_StmMcwi_2023
    FXLE13_StmMcwi_2023
    FXLE14_StmMcwi_2023
    FXLE15_StmMcwi_2023
    FXLE16_Stm_2020
    FXLE17_StmMcwi_2023
    FXLE18_Stm_2020
    FXLE19_StmMcwi_2022
    FXLE20_StmMcwi_2023
    FXLE21_StmMcwi_2023
    FXLE22_StmMcwi_2023
    FXLE23_StmMcwi_2023
    FXLE25_StmMcwi_2023
    FXLE26_StmMcwi_2023
    GK_FarMcwi_2019
    HXB2_IpcvMcwi_2019
    HXB4_IpcvMcwi_2020
    HXB10_IpcvMcwi_2019
    HXB17_IpcvMcwi_2021
    HXB18_IpcvMcwi_2021
    HXB20_IpcvMcwi_2020
    HXB23_IpcvMcwi_2021
    HXB31_IpcvMcwi_2019
    LEW_Crl_2019
    LEXF1A_Stm_2019
    LEXF1C_StmMcwi_2021
    LEXF2A_StmMcwi_2023
    LEXF2B_Stm_2019
    LEXF2C_StmMcwi_2022
    LEXF3_Stm_2020
    LEXF4_Stm_2020
    LEXF5_StmMcwi_2023
    LEXF6B_StmMcwi_2023
    LEXF7A_Stm_2019NG
    LEXF7B_StmMcwi_2022
    LEXF7C_StmMcwi_2022
    LEXF8A_StmMcwi_2021
    LEXF8D_StmMcwi_2023
    LEXF9_StmMcwi_2022
    LEXF10A_StmMcwi_2020
    LEXF10B_StmMcwi_2022
    LEXF10C_StmMcwi_2023
    LEXF11_Stm_2020
    LE_Stm_2019
    LH_MavRrrcAek_2020
    LL_MavRrrcAek_2020
    LN_MavRrrcAek_2020
    M520_N_2020
    MR_N_2020
    MWF_SimwMcwi_2019
    PVG_Seac_2019
    RCS_LavRrrcMcwi_2021
    SHRSP_A3NCrl_2019
    SHR_NCrl_2021
    SHR_NHsd_2021
    SHR_OlalpcvMcwi_2019
    SR_JrHsd_2020
    SS_JrHsdMcwi_2019
    SS_HsdMcwiCrl_2021
    WAG_RijCrl_2020
    WKY_NCrl_2019
    WKY_NHsd_2022
    WN_N_2020
)

echo "============================================"
echo "HRDP VCF Download — ${#STRAINS[@]} strains"
echo "============================================"

N_PARALLEL=5
RUNNING=0
DONE=0
FAIL=0

for STRAIN in "${STRAINS[@]}"; do
    FILE="${STRAIN}_rnBN7_SNPs_PASS.vcf.gz"
    if [ -f "${OUT}/${FILE}" ] && [ -s "${OUT}/${FILE}" ]; then
        ((DONE++))
        continue
    fi

    while [ ${RUNNING} -ge ${N_PARALLEL} ]; do
        wait -n 2>/dev/null || true
        RUNNING=$((RUNNING - 1))
    done

    (
        wget -q -O "${OUT}/${FILE}" "${BASE_URL}/${FILE}" && \
        wget -q -O "${OUT}/${FILE}.tbi" "${BASE_URL}/${FILE}.tbi" && \
        echo "[OK] ${STRAIN}" || \
        { echo "[FAIL] ${STRAIN}"; rm -f "${OUT}/${FILE}" "${OUT}/${FILE}.tbi"; }
    ) &
    RUNNING=$((RUNNING + 1))
done

wait
DONE=$(ls "${OUT}"/*.vcf.gz 2>/dev/null | wc -l)
echo ""
echo "============================================"
echo "완료: ${DONE}/${#STRAINS[@]} strains"
echo "Disk: $(df -h /home/yusanghyeon | tail -1 | awk '{print $4}') free"
echo "============================================"
