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
base_seed=20260901
campaign_tag="gof_root640_chi2_integral_binningscan_v1_20260901"
analysis_dir="/eos/home-l/leyao/26JJ/JPsiJPsi/Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer"

if (( job_index < 0 || job_index >= number_of_jobs )); then
    echo "JOB_INDEX must be in [0,$((number_of_jobs - 1))]" >&2
    exit 2
fi

projection_bins=$((minimum_bins + bin_step * job_index))
result_tag="${campaign_tag}_bins${projection_bins}"
result_dir="${analysis_dir}/fit_results/${result_tag}"

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
    echo "Refusing to overwrite existing result: ${result_dir}" >&2
    exit 2
fi

echo "Starting chi2 binning scan job ${job_index}: ${projection_bins} bins"
root -l -b -q \
    "Fit_Check.cpp+(0,0,\"${result_tag}\",${base_seed},${projection_bins},false)"

for required in configuration.txt observed_chi2.csv projection_bins.csv \
    toy_results.csv fit_attempts.csv; do
    if [[ ! -s "${result_dir}/${required}" ]]; then
        echo "Missing or empty result file: ${result_dir}/${required}" >&2
        exit 1
    fi
done

observed_rows="$(awk 'END {print NR - 1}' "${result_dir}/observed_chi2.csv")"
projection_rows="$(awk 'END {print NR - 1}' "${result_dir}/projection_bins.csv")"
if [[ "${observed_rows}" -ne 4 ||
      "${projection_rows}" -ne $((4 * projection_bins)) ]]; then
    echo "CSV validation failed: observed_rows=${observed_rows}, projection_rows=${projection_rows}" >&2
    exit 1
fi

expected_total="$(awk -F= '$1 == "nominal_expected_events" {print $2}' \
    "${result_dir}/configuration.txt")"
awk -F, -v expected_total="${expected_total}" -v bins="${projection_bins}" '
    NR > 1 {
        sums[$1] += $8
        counts[$1]++
    }
    END {
        projections = 0
        for(name in sums) {
            projections++
            difference = sums[name] - expected_total
            if(difference < 0) difference = -difference
            if(counts[name] != bins || difference / expected_total > 1e-7) {
                printf "Expected-yield closure failed for %s: bins=%d sum=%.17g expected=%.17g\n",
                    name, counts[name], sums[name], expected_total > "/dev/stderr"
                failed = 1
            }
        }
        if(projections != 4 || failed) exit 1
    }
' "${result_dir}/projection_bins.csv"

plot_count="$(find "${result_dir}/plots" -maxdepth 1 -type f \
    \( -name 'projection_*.pdf' -o -name 'projection_*.png' \) | wc -l)"
if [[ "${plot_count}" -ne 8 ]]; then
    echo "Plot validation failed: expected 8 PDF/PNG files, found ${plot_count}" >&2
    exit 1
fi

metadata_tmp="${result_dir}/job.metadata.txt.tmp.$$"
{
    echo "status=complete"
    echo "campaign_tag=${campaign_tag}"
    echo "job_index=${job_index}"
    echo "result_tag=${result_tag}"
    echo "projection_bins=${projection_bins}"
    echo "toys=0"
    echo "numeric_outputs=observed_chi2.csv,projection_bins.csv"
    echo "plot_outputs=4 projection PDFs and 4 projection PNGs"
    echo "projection_integration=adaptive Gauss-Kronrod interval integral normalized by the corresponding full-range integral"
    echo "integration_tolerances=absolute 1e-10; relative 1e-9; maximum 100000 subintervals"
    echo "root_version=$(root-config --version)"
    echo "completed_at=$(date --iso-8601=seconds)"
    sha256sum Fit_Check.cpp Fit_Check_cpp.so WeightData.root Model_4D_tot.root \
        Fit_4D_tot_native_asymptotic_unseeded.root \
        "${result_dir}/configuration.txt" \
        "${result_dir}/observed_chi2.csv" \
        "${result_dir}/projection_bins.csv"
} > "${metadata_tmp}"
mv "${metadata_tmp}" "${result_dir}/job.metadata.txt"
echo "Completed chi2 binning scan job ${job_index}: ${result_dir}"
