#!/usr/bin/env bash
set -euo pipefail

analysis="/eos/home-l/leyao/26JJ/JPsiJPsi/Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer"
source_tag="newdata167_accmix23_effmix19_sps0p85_dps0p15_lumi36p684_pairsym_cseed50_root640_strategy2_nominal_v4_20260902"
campaign_tag="newdata167_accmix23_effmix19_sps0p85_dps0p15_lumi36p684_pairsym_cseed50_root640_strategy2_diff35_v1_20260902"
weight_name="WeightData_newdata167_accmix23_effmix19_sps0p85_dps0p15_pairsym_cseed50_v1_20260902.root"
source_dir="${analysis}/fit_results/${source_tag}"
result_dir="${analysis}/fit_results/${campaign_tag}"
total_dir="${result_dir}/total_models"

if [[ "$(root-config --version)" != "6.40.02" ]]; then
    echo "This workflow requires ROOT 6.40.02" >&2
    exit 3
fi
required=(
    "${source_dir}/inputs/${weight_name}"
    "${source_dir}/total/Model_4D_tot.root"
    "${source_dir}/total/Fit_4D_tot_native_asymptotic_unseeded.root"
    "${analysis}/Fit_4D_tot.cpp"
    "${analysis}/Fit_valid.cpp"
    "${analysis}/Plot_4D.hpp"
)
for path in "${required[@]}"; do
    test -s "${path}"
done
if [[ -e "${result_dir}" ]]; then
    echo "Refusing to overwrite existing campaign: ${result_dir}" >&2
    exit 17
fi

mkdir -p "${result_dir}/inputs" "${total_dir}/logs" "${result_dir}/jobs" "${result_dir}/rows" "${result_dir}/logs"
cp "${source_dir}/inputs/${weight_name}" "${result_dir}/inputs/"
cp "${source_dir}/total/Model_4D_tot.root" "${source_dir}/total/Fit_4D_tot_native_asymptotic_unseeded.root" "${total_dir}/"
cp "${analysis}/Fit_4D_tot.cpp" "${analysis}/Plot_4D.hpp" "${total_dir}/"
ln -s "../inputs/${weight_name}" "${total_dir}/WeightData.root"

cd "${total_dir}"
reference_strategy=2
if ! root -l -b -q 'Fit_4D_tot.cpp(true,false,false,true,false,2)' > logs/reference_strategy2.log 2>&1 ||
   [[ ! -s Model_4D_tot_ref.root ||
      ! -s Fit_4D_tot_ref_native_asymptotic_unseeded.root ]]; then
    reference_strategy=1
    root -l -b -q 'Fit_4D_tot.cpp(true,false,false,true,false,1)' > logs/reference_strategy1.log 2>&1
fi
test -s Model_4D_tot_ref.root
test -s Fit_4D_tot_ref_native_asymptotic_unseeded.root

cd "${analysis}"
root -l -b -q "Fit_valid.cpp+(\"total\",0,0,\"${campaign_tag}\",\"${total_dir}\")" > "${result_dir}/logs/total_validation.log" 2>&1
tail -n 1 "${result_dir}/logs/total_validation.log" > "${result_dir}/total_row.csv"
awk -F, 'NF != 27 || $1 != "total" || $2 != "total" ||
    $22 != 0 || $23 != 3 || $25 != 0 || $26 != 3 {exit 1}' "${result_dir}/total_row.csv"

{
    echo "status=prepared"
    echo "campaign_tag=${campaign_tag}"
    echo "source_total_tag=${source_tag}"
    echo "weightdata=${result_dir}/inputs/${weight_name}"
    echo "nominal_model=${total_dir}/Model_4D_tot.root"
    echo "reference_model=${total_dir}/Model_4D_tot_ref.root"
    echo "reference_total_fit=ROOT native AsymptoticError(true), Strategy(${reference_strategy}), no external seed"
    echo "reference_strategy_policy=try Strategy(2); if it fails before producing accepted artifacts, retry Strategy(1)"
    echo "differential_fit=all total shape parameters fixed; five yields float; failed fits may retry with Strategy(2)+Offset(true); failed lower-boundary n_Sig_Comb and/or n_Comb_Comb may be fixed to zero"
    echo "systematic_status=base_fit_table_only; formal result is fit_results/_CURRENT/systematics"
    echo "base_correction_systematic=legacy 0.85:0.15 sys3; total relative 0.08805293089840492"
    echo "root_version=$(root-config --version)"
    echo "prepared_at=$(date --iso-8601=seconds)"
    sha256sum "${result_dir}/inputs/${weight_name}" "${total_dir}/Model_4D_tot.root" "${total_dir}/Fit_4D_tot_native_asymptotic_unseeded.root" "${total_dir}/Model_4D_tot_ref.root" "${total_dir}/Fit_4D_tot_ref_native_asymptotic_unseeded.root" "${result_dir}/total_row.csv" Fit_4D_diff.cpp Fit_valid.cpp Plot_4D.hpp collect_differential_newdata167.py run_differential_mix0p85_pairsym_condor.sh
} > "${result_dir}/prepare.metadata.txt"

echo "DIFFERENTIAL_CAMPAIGN_PREPARED result=${result_dir}"
