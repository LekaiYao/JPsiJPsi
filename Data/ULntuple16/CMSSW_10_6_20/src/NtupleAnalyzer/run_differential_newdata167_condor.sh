#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 || ! "$1" =~ ^[0-9]+$ ]]; then
    echo "Usage: $0 JOB_INDEX" >&2
    exit 2
fi

job_index="$1"
analysis_dir="/eos/home-l/leyao/26JJ/JPsiJPsi/Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer"
upstream_tag="newdata167_accmix23_effmix19_lumi36p684_root640_v1_20260902"
campaign_tag="newdata167_accmix23_effmix19_lumi36p684_root640_diff35_v1_20260902"
weight_data="${analysis_dir}/fit_results/${upstream_tag}/inputs/WeightData_newdata167_accmix23_effmix19_v1_20260902.root"
total_dir="${analysis_dir}/fit_results/${upstream_tag}/total"
result_dir="${analysis_dir}/fit_results/${campaign_tag}"

bins=(
    "delta_y 0 0.5"
    "delta_y 0.5 1"
    "delta_y 1 1.5"
    "delta_y 1.5 2"
    "delta_y 2 2.5"
    "delta_y 2.5 4"
    "delta_phi 0 0.3927"
    "delta_phi 0.3927 0.7854"
    "delta_phi 0.7854 1.1781"
    "delta_phi 1.1781 1.5708"
    "delta_phi 1.5708 1.9635"
    "delta_phi 1.9635 2.3562"
    "delta_phi 2.3562 2.7489"
    "delta_phi 2.7489 3.1416"
    "evt_mass 7.5 17.5"
    "evt_mass 17.5 27.5"
    "evt_mass 27.5 37.5"
    "evt_mass 37.5 47.5"
    "evt_mass 47.5 57.5"
    "evt_mass 57.5 67.5"
    "evt_mass 67.5 107.5"
    "evt_y 0 0.4"
    "evt_y 0.4 0.8"
    "evt_y 0.8 1.2"
    "evt_y 1.2 1.6"
    "evt_y 1.6 2.0"
    "evt_pt 0 5"
    "evt_pt 5 10"
    "evt_pt 10 15"
    "evt_pt 15 20"
    "evt_pt 20 25"
    "evt_pt 25 30"
    "evt_pt 30 35"
    "evt_pt 35 40"
    "evt_pt 40 80"
)

if (( job_index < 0 || job_index >= ${#bins[@]} )); then
    echo "JOB_INDEX must be in [0,$((${#bins[@]} - 1))]" >&2
    exit 2
fi

read -r variable minimum maximum <<< "${bins[${job_index}]}"
printf -v job_suffix "%02d" "${job_index}"
row_file="${result_dir}/rows/${job_suffix}.csv"
metadata_file="${result_dir}/jobs/${job_suffix}.metadata.txt"
if [[ -e "${row_file}" || -e "${metadata_file}" ]]; then
    echo "Refusing to overwrite differential job ${job_index}" >&2
    exit 2
fi

cd "${analysis_dir}"
if [[ "$(root-config --version)" != "6.40.02" ]]; then
    echo "This workflow requires ROOT 6.40.02" >&2
    exit 2
fi
required_files=("${weight_data}" "${total_dir}/Model_4D_tot.root" "${total_dir}/Model_4D_tot_ref.root" "${total_dir}/Fit_4D_tot_native_asymptotic_unseeded.root" "${total_dir}/Fit_4D_tot_ref_native_asymptotic_unseeded.root" Fit_4D_diff.cpp Fit_4D_diff_cpp.so Fit_valid.cpp Fit_valid_cpp.so Plot_4D.hpp)
for required in "${required_files[@]}"; do
    test -s "${required}"
done
if [[ Fit_4D_diff.cpp -nt Fit_4D_diff_cpp.so ||
      Fit_valid.cpp -nt Fit_valid_cpp.so ||
      Plot_4D.hpp -nt Fit_4D_diff_cpp.so ]]; then
    echo "A compiled macro is stale; refusing worker-side ACLiC compilation" >&2
    exit 2
fi

nominal_log="${result_dir}/logs/${job_suffix}_${variable}_nominal.log"
reference_log="${result_dir}/logs/${job_suffix}_${variable}_reference.log"
validation_log="${result_dir}/logs/${job_suffix}_${variable}_validation.log"
root -l -b -q "Fit_4D_diff.cpp+(\"${variable}\",${minimum},${maximum},false,true,\"${campaign_tag}\",\"${weight_data}\",\"${total_dir}\",false)" > "${nominal_log}" 2>&1
root -l -b -q "Fit_4D_diff.cpp+(\"${variable}\",${minimum},${maximum},true,true,\"${campaign_tag}\",\"${weight_data}\",\"${total_dir}\",false)" > "${reference_log}" 2>&1
root -l -b -q "Fit_valid.cpp+(\"${variable}\",${minimum},${maximum},\"${campaign_tag}\",\"${total_dir}\")" > "${validation_log}" 2>&1
tail -n 1 "${validation_log}" > "${row_file}.tmp"
awk -F, -v variable="${variable}" -v minimum="${minimum}" -v maximum="${maximum}" '
    NF != 27 || $1 != "differential" || $2 != variable ||
    $3 + 0 != minimum + 0 || $4 + 0 != maximum + 0 ||
    $22 != 0 || $23 != 3 || $25 != 0 || $26 != 3 {exit 1}
' "${row_file}.tmp"
mv "${row_file}.tmp" "${row_file}"

nominal_result="$(sed -n 's/^Accepted differential result: \([^,]*\),.*/\1/p' "${nominal_log}" | tail -n 1)"
reference_result="$(sed -n 's/^Accepted differential result: \([^,]*\),.*/\1/p' "${reference_log}" | tail -n 1)"
test -n "${nominal_result}" && test -n "${reference_result}"
test -s "${analysis_dir}/${nominal_result}"
test -s "${analysis_dir}/${reference_result}"

{
    echo "status=complete"
    echo "campaign_tag=${campaign_tag}"
    echo "job_index=${job_index}"
    echo "variable=${variable}"
    echo "bin_min=${minimum}"
    echo "bin_max=${maximum}"
    echo "nominal_result=${nominal_result}"
    echo "reference_result=${reference_result}"
    echo "error_convention=ROOT native AsymptoticError(true), no external seed"
    echo "fit_gate=status=0,covQual=3,EDM<0.01"
    echo "root_version=$(root-config --version)"
    echo "completed_at=$(date --iso-8601=seconds)"
    sha256sum "${weight_data}" "${total_dir}/Model_4D_tot.root" "${total_dir}/Model_4D_tot_ref.root" Fit_4D_diff.cpp Fit_4D_diff_cpp.so Fit_valid.cpp Fit_valid_cpp.so Plot_4D.hpp "${analysis_dir}/${nominal_result}" "${analysis_dir}/${reference_result}" "${row_file}"
} > "${metadata_file}.tmp"
mv "${metadata_file}.tmp" "${metadata_file}"
echo "DIFFERENTIAL_JOB_COMPLETE index=${job_index} variable=${variable} range=${minimum},${maximum}"
