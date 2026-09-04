#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 || ! "$1" =~ ^[0-9]+$ ]]; then
    echo "Usage: $0 JOB_INDEX" >&2
    exit 2
fi

job_index="$1"
number_of_jobs=100
toys_per_job=50
global_base_seed=20260902
seed_stride=104729
projection_bins=20
campaign_tag="gof_root640_toys5000_yieldonly_ncomb0fallback_newdata167_accmix23_effmix19_v1_20260902"
analysis_dir="/eos/home-l/leyao/26JJ/JPsiJPsi/Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer"
upstream_tag="newdata167_accmix23_effmix19_lumi36p684_root640_v1_20260902"
weight_data="${analysis_dir}/fit_results/${upstream_tag}/inputs/WeightData_newdata167_accmix23_effmix19_v1_20260902.root"
model_file="${analysis_dir}/fit_results/${upstream_tag}/total/Model_4D_tot.root"
fit_result_file="${analysis_dir}/fit_results/${upstream_tag}/total/Fit_4D_tot_native_asymptotic_unseeded.root"

if (( job_index < 0 || job_index >= number_of_jobs )); then
    echo "JOB_INDEX must be in [0,$((number_of_jobs - 1))]" >&2
    exit 2
fi

global_start=$((job_index * toys_per_job))
global_end=$((global_start + toys_per_job - 1))
job_base_seed=$((global_base_seed + seed_stride * global_start))
printf -v job_suffix "%03d" "${job_index}"
job_tag="${campaign_tag}_job${job_suffix}"
result_dir="${analysis_dir}/fit_results/${job_tag}"

cd "${analysis_dir}"
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
    echo "Refusing to overwrite existing shard: ${result_dir}" >&2
    exit 2
fi

root -l -b -q "Fit_Check.cpp+(0,${toys_per_job},\"${job_tag}\",${job_base_seed},${projection_bins},false,\"${weight_data}\",\"${model_file}\",\"${fit_result_file}\")"

results_file="${result_dir}/toy_results.csv"
result_rows="$(awk 'END {print NR - 1}' "${results_file}")"
accepted_rows="$(awk -F, 'NR > 1 && $7 == 1 {count++} END {print count + 0}' "${results_file}")"
if [[ "${result_rows}" -ne "${toys_per_job}" ||
      "${accepted_rows}" -ne "${toys_per_job}" ]]; then
    echo "Shard validation failed: rows=${result_rows}, accepted=${accepted_rows}" >&2
    exit 1
fi
if [[ -e "${result_dir}/toy_0_fit_plots" ||
      -e "${result_dir}/toy_0_fit.root" ]]; then
    echo "Production shard unexpectedly created toy-0 diagnostics" >&2
    exit 1
fi

{
    echo "status=complete"
    echo "campaign_tag=${campaign_tag}"
    echo "job_index=${job_index}"
    echo "job_tag=${job_tag}"
    echo "global_toy_start=${global_start}"
    echo "global_toy_end=${global_end}"
    echo "local_toy_start=0"
    echo "local_toy_end=$((toys_per_job - 1))"
    echo "toys_requested=${toys_per_job}"
    echo "toys_accepted=${accepted_rows}"
    echo "global_base_seed=${global_base_seed}"
    echo "job_base_seed=${job_base_seed}"
    echo "seed_mapping=job_base_seed=global_base_seed+104729*global_toy_start"
    echo "toy0_diagnostics=disabled"
    echo "fit_parameters=yield-only; all nominal shape parameters fixed"
    echo "n_Comb_Comb_rule=lower boundary 0; after a failed boundary fit fix at 0 and refit"
    echo "fit_gate=status=0,covQual=3,EDM<0.01,n_P_P error finite and positive"
    echo "root_version=$(root-config --version)"
    echo "completed_at=$(date --iso-8601=seconds)"
    sha256sum Fit_Check.cpp Fit_Check_cpp.so "${weight_data}" "${model_file}" "${fit_result_file}" run_gof_toys5000_newdata167_condor.sh "${result_dir}/configuration.txt" "${results_file}" "${result_dir}/fit_attempts.csv"
} > "${result_dir}/job.metadata.txt"
echo "TOY_JOB_COMPLETE index=${job_index} range=${global_start}-${global_end}"
