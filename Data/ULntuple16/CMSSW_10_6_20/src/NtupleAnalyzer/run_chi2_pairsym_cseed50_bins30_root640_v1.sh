#!/usr/bin/env bash
set -euo pipefail

analysis="/eos/home-l/leyao/26JJ/JPsiJPsi/Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer"
upstream_tag="newdata167_accmix23_effmix19_lumi36p684_pairsym_cseed50_root640_strategy2_nominalonly_v1_20260902"
result_tag="gof_root640_chi2_integral_newdata167_pairsym_cseed50_strategy2_nominalonly_bins30_v1_20260902"
result_dir="${analysis}/fit_results/${result_tag}"
weight_data="${analysis}/fit_results/${upstream_tag}/inputs/WeightData_newdata167_accmix23_effmix19_pairsym_cseed50_v1_20260902.root"
model_file="${analysis}/fit_results/${upstream_tag}/total/Model_4D_tot.root"
fit_result_file="${analysis}/fit_results/${upstream_tag}/total/Fit_4D_tot_native_asymptotic_unseeded.root"
projection_bins=30
base_seed=20260902

cd "${analysis}"
if [[ "$(root-config --version)" != "6.40.02" ]]; then
    echo "This workflow requires ROOT 6.40.02" >&2
    exit 2
fi
required_files=("${weight_data}" "${model_file}" "${fit_result_file}" Fit_Check.cpp Fit_Check_cpp.so)
for required in "${required_files[@]}"; do
    test -s "${required}"
done
if [[ Fit_Check.cpp -nt Fit_Check_cpp.so ]]; then
    echo "Fit_Check_cpp.so is stale; refusing implicit ACLiC rebuild" >&2
    exit 2
fi
if [[ -e "${result_dir}" ]]; then
    echo "Refusing to overwrite existing result: ${result_dir}" >&2
    exit 17
fi

root -l -b -q "Fit_Check.cpp+(0,0,\"${result_tag}\",${base_seed},${projection_bins},false,\"${weight_data}\",\"${model_file}\",\"${fit_result_file}\",false)"

required_outputs=(configuration.txt observed_chi2.csv projection_bins.csv toy_results.csv fit_attempts.csv)
for required in "${required_outputs[@]}"; do
    test -s "${result_dir}/${required}"
done
observed_rows="$(awk 'END {print NR - 1}' "${result_dir}/observed_chi2.csv")"
projection_rows="$(awk 'END {print NR - 1}' "${result_dir}/projection_bins.csv")"
test "${observed_rows}" -eq 4
test "${projection_rows}" -eq $((4 * projection_bins))
grep -Fxq "projection_pair_symmetrization=disabled: Jpsi1 and Jpsi2 projected separately" "${result_dir}/configuration.txt"

expected_total="$(awk -F= '$1 == "nominal_expected_events" {print $2}' "${result_dir}/configuration.txt")"
awk -F, -v expected_total="${expected_total}" -v bins="${projection_bins}" '
    NR > 1 {sums[$1] += $8; counts[$1]++}
    END {
        for(name in sums) {
            projections++
            difference=sums[name]-expected_total
            if(difference<0) difference=-difference
            if(counts[name]!=bins || difference/expected_total>1e-7) failed=1
        }
        if(projections!=4 || failed) exit 1
    }
' "${result_dir}/projection_bins.csv"

plot_count="$(find "${result_dir}/plots" -maxdepth 1 -type f \( -name 'projection_*.pdf' -o -name 'projection_*.png' \) | wc -l)"
test "${plot_count}" -eq 8
{
    echo "status=complete"
    echo "task=chi2_newdata167_pairsym_cseed50_strategy2_nominalonly_bins30_root640_v1"
    echo "result_tag=${result_tag}"
    echo "projection_bins=${projection_bins}"
    echo "tree_label_randomization=chensh-compatible C srand/rand seed 50, reseeded per input ROOT file"
    echo "projection_pair_symmetrization=disabled"
    echo "total_fit_minimizer_strategy=2"
    echo "fit_and_weights=from ${upstream_tag}"
    echo "projection_integration=adaptive Gauss-Kronrod interval integral normalized by the corresponding full-range integral"
    echo "integration_tolerances=absolute 1e-10; relative 1e-9; maximum 100000 subintervals"
    echo "root_version=$(root-config --version)"
    echo "completed_at=$(date --iso-8601=seconds)"
    sha256sum Fit_Check.cpp Fit_Check_cpp.so "${weight_data}" "${model_file}" "${fit_result_file}" "${result_dir}/configuration.txt" "${result_dir}/observed_chi2.csv" "${result_dir}/projection_bins.csv"
} > "${result_dir}/run.metadata.txt"

echo "CHI2_PAIRSYM_COMPLETE result=${result_dir}/observed_chi2.csv"
