#!/usr/bin/env bash
set -euo pipefail

run_tag="${1:-gof_root640_v1}"
number_of_toys="${2:-5000}"
chunk_size="${3:-50}"
base_seed="${4:-20260901}"
result_dir="fit_results/${run_tag}"

if [[ ! "${run_tag}" =~ ^[A-Za-z0-9_-]+$ ]]; then
    echo "Invalid run tag: ${run_tag}" >&2
    exit 2
fi
if (( number_of_toys <= 0 || chunk_size <= 0 )); then
    echo "Toy and chunk counts must be positive" >&2
    exit 2
fi
if [[ -e "${result_dir}" ]]; then
    echo "Refusing to overwrite an existing run: ${result_dir}" >&2
    exit 2
fi

mkdir -p "${result_dir}/logs"
start=0
while (( start < number_of_toys )); do
    count="${chunk_size}"
    if (( start + count > number_of_toys )); then
        count=$((number_of_toys - start))
    fi
    end=$((start + count - 1))
    echo "Toys ${start}-${end} / $((number_of_toys - 1))"
    root -l -b -q \
        "Fit_Check.cpp+(${start},${count},\"${run_tag}\",${base_seed})" \
        > "${result_dir}/logs/toys_${start}_${end}.log" 2>&1
    start=$((start + count))
done

root -l -b -q "Fit_pull.cpp+(\"${run_tag}\",${number_of_toys})" \
    > "${result_dir}/logs/pull_summary.log" 2>&1

{
    echo "status=complete"
    echo "run_tag=${run_tag}"
    echo "number_of_toys=${number_of_toys}"
    echo "chunk_size=${chunk_size}"
    echo "base_seed=${base_seed}"
    echo "root_version=$(root-config --version)"
    echo "git_head=$(git rev-parse HEAD)"
    echo "generated_at=$(date --iso-8601=seconds)"
    echo "core_validation=historical toy-fit pull distribution redone with weighted Poisson-resampled toys and ROOT native AsymptoticError(true)"
    echo "auxiliary_chi2=four 20-bin nominal-PDF-quantile projections computed from each unbinned fit; no binned refit and empirical toy p-values"
    sha256sum WeightData.root Model_4D_tot.root \
        Fit_4D_tot_native_asymptotic_unseeded.root \
        Fit_Check.cpp Fit_pull.cpp run_gof_toys_root640.sh \
        "${result_dir}/configuration.txt" \
        "${result_dir}/projection_bins.csv" \
        "${result_dir}/observed_chi2.csv" \
        "${result_dir}/toy_results.csv" \
        "${result_dir}/gof_summary.txt" \
        "${result_dir}/pull_fit.root"
} > "${result_dir}/run.metadata.txt"

echo "Completed ${number_of_toys} toys: ${result_dir}"
