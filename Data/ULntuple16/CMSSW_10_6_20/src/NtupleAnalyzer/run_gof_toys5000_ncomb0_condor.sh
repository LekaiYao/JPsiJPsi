#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 || ! "$1" =~ ^[0-9]+$ ]]; then
    echo "Usage: $0 JOB_INDEX" >&2
    exit 2
fi

job_index="$1"
number_of_jobs=100
toys_per_job=50
global_base_seed=20260901
seed_stride=104729
projection_bins=20
make_toy0_diagnostics=false
campaign_tag="gof_root640_toys5000_yieldonly_ncomb0fallback_condor_v1_20260901"
analysis_dir="/eos/home-l/leyao/26JJ/JPsiJPsi/Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer"

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
for required in WeightData.root Model_4D_tot.root \
    Fit_4D_tot_native_asymptotic_unseeded.root Fit_Check.cpp Fit_Check_cpp.so; do
    if [[ ! -e "${required}" ]]; then
        echo "Missing required input: ${required}" >&2
        exit 2
    fi
done
if [[ Fit_Check.cpp -nt Fit_Check_cpp.so ]]; then
    echo "Fit_Check_cpp.so is older than Fit_Check.cpp; refusing worker-side ACLiC compilation" >&2
    exit 2
fi
if [[ -e "${result_dir}" ]]; then
    echo "Refusing to overwrite existing shard: ${result_dir}" >&2
    exit 2
fi

echo "Starting job ${job_index}: global toys ${global_start}-${global_end}, job base seed ${job_base_seed}"
root -l -b -q \
    "Fit_Check.cpp+(0,${toys_per_job},\"${job_tag}\",${job_base_seed},${projection_bins},${make_toy0_diagnostics})"

results_file="${result_dir}/toy_results.csv"
if [[ ! -s "${results_file}" ]]; then
    echo "Missing toy results: ${results_file}" >&2
    exit 1
fi
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

metadata_tmp="${result_dir}/job.metadata.txt.tmp.$$"
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
    sha256sum Fit_Check.cpp Fit_Check_cpp.so WeightData.root Model_4D_tot.root \
        Fit_4D_tot_native_asymptotic_unseeded.root \
        run_gof_toys5000_ncomb0_condor.sh \
        "${result_dir}/configuration.txt" "${results_file}" \
        "${result_dir}/fit_attempts.csv"
} > "${metadata_tmp}"
mv "${metadata_tmp}" "${result_dir}/job.metadata.txt"
echo "Completed job ${job_index}: ${result_dir}"
