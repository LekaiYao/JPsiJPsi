#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 || ! "$1" =~ ^[0-9]+$ ]]; then
    echo "Usage: $0 JOB_INDEX" >&2
    exit 2
fi

job_index="$1"
number_of_jobs=11
minimum_bins=10
bin_step=5
base_seed=20260902
campaign_tag="gof_root640_chi2_integral_newdata167_sps0p85_dps0p15_pairsym_cseed50_strategy2_v1_20260902"
analysis_dir="/eos/home-l/leyao/26JJ/JPsiJPsi/Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer"
upstream_tag="newdata167_accmix23_effmix19_sps0p85_dps0p15_lumi36p684_pairsym_cseed50_root640_strategy2_nominal_v4_20260902"
weight_name="WeightData_newdata167_accmix23_effmix19_sps0p85_dps0p15_pairsym_cseed50_v1_20260902.root"
weight_data="${analysis_dir}/fit_results/${upstream_tag}/inputs/${weight_name}"
model_file="${analysis_dir}/fit_results/${upstream_tag}/total/Model_4D_tot.root"
fit_result_file="${analysis_dir}/fit_results/${upstream_tag}/total/Fit_4D_tot_native_asymptotic_unseeded.root"

if (( job_index < 0 || job_index >= number_of_jobs )); then
    echo "JOB_INDEX must be in [0,$((number_of_jobs - 1))]" >&2
    exit 2
fi

projection_bins=$((minimum_bins + bin_step * job_index))
result_tag="${campaign_tag}_bins${projection_bins}"
result_dir="${analysis_dir}/fit_results/${result_tag}"

cd "${analysis_dir}"
worker_cache="${_CONDOR_SCRATCH_DIR:-/tmp}/xdg-cache-chi2-${job_index}"
mkdir -p "${worker_cache}"
export XDG_CACHE_HOME="${worker_cache}"
if [[ "$(root-config --version)" != "6.40.02" ]]; then
    echo "This workflow requires ROOT 6.40.02" >&2
    exit 2
fi
required_files=("${weight_data}" "${model_file}" "${fit_result_file}" Fit_Check.cpp Fit_Check_cpp.so)
for required in "${required_files[@]}"; do
    test -s "${required}"
done
if [[ Fit_Check.cpp -nt Fit_Check_cpp.so ]]; then
    echo "Fit_Check_cpp.so is stale; refusing worker-side ACLiC compilation" >&2
    exit 2
fi
if [[ -e "${result_dir}" ]]; then
    echo "Refusing to overwrite existing result: ${result_dir}" >&2
    exit 2
fi

root -l -b -q -e "gSystem->Load(\"Fit_Check_cpp.so\"); Fit_Check(0,0,\"${result_tag}\",${base_seed},${projection_bins},false,\"${weight_data}\",\"${model_file}\",\"${fit_result_file}\");"

required_outputs=(configuration.txt observed_chi2.csv projection_bins.csv toy_results.csv fit_attempts.csv)
for required in "${required_outputs[@]}"; do
    test -s "${result_dir}/${required}"
done
test "$(awk 'END {print NR - 1}' "${result_dir}/observed_chi2.csv")" -eq 4
test "$(awk 'END {print NR - 1}' "${result_dir}/projection_bins.csv")" -eq $((4 * projection_bins))
plot_count="$(find "${result_dir}/plots" -maxdepth 1 -type f \( -name 'projection_*.pdf' -o -name 'projection_*.png' \) | wc -l)"
test "${plot_count}" -eq 8

{
    echo "status=complete"
    echo "campaign_tag=${campaign_tag}"
    echo "job_index=${job_index}"
    echo "result_tag=${result_tag}"
    echo "projection_bins=${projection_bins}"
    echo "toys=0"
    echo "input_label_symmetrization=chensh-compatible C srand/rand seed 50, reseeded per input ROOT file"
    echo "projection_pair_symmetrization=disabled"
    echo "projection_sumw2=standard event-level"
    echo "projection_integration=adaptive Gauss-Kronrod interval integral normalized by the corresponding full-range integral"
    echo "integration_tolerances=absolute 1e-10; relative 1e-9; maximum 100000 subintervals"
    echo "worker_runtime=precompiled ROOT library loaded directly; per-job XDG cache"
    echo "root_version=$(root-config --version)"
    echo "completed_at=$(date --iso-8601=seconds)"
    sha256sum Fit_Check.cpp Fit_Check_cpp.so "${weight_data}" "${model_file}" "${fit_result_file}" "${result_dir}/configuration.txt" "${result_dir}/observed_chi2.csv" "${result_dir}/projection_bins.csv"
} > "${result_dir}/job.metadata.txt"
echo "CHI2_JOB_COMPLETE index=${job_index} bins=${projection_bins}"
