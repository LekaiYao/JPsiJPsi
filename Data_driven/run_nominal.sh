#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <unique-run-tag>" >&2
  exit 64
fi

TAG="$1"
[[ "${TAG}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || {
  echo "Unsafe run tag: ${TAG}" >&2
  exit 64
}

REPO="/eos/home-l/leyao/26JJ/JPsiJPsi"
source "${REPO}/Data_driven/require_root640.sh"
RESULT="${REPO}/Data_driven/results/${TAG}"
MODEL="${REPO}/Data_driven/inputs/Model_4D_tot.root"
RAW_DATA_BASE="${JJ_DATA_NTUPLE_BASE:-/eos/user/c/chensh/JPsiJPsi/Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer}"
INPUT_DIR="${RESULT}/input"
SPLOT_DIR="${RESULT}/splot"
PP_DIR="${RESULT}/pp_data_2d"
TEMPLATE_DIR="${RESULT}/templates_2d"
WEIGHT_DATA="${INPUT_DIR}/WeightData_mJJ7p5.root"
CANDIDATES="${INPUT_DIR}/jpsi12_candidates.root"
MIXED="${RESULT}/mixed_dps_allpairs.root"

if [ -e "${RESULT}" ]; then
  echo "Refusing to overwrite existing tag: ${RESULT}" >&2
  exit 68
fi
mkdir -p "${INPUT_DIR}" "${SPLOT_DIR}" "${PP_DIR}" "${TEMPLATE_DIR}" "${RESULT}/logs"
cd "${REPO}"

cd "${REPO}/Data_driven/inputs"
root -l -b -q \
  "${REPO}/Data_driven/build_uniform_data_input.cpp(\"${WEIGHT_DATA}\")" \
  > "${RESULT}/logs/build_uniform_data.log" 2>&1
cd "${REPO}"
root -l -b -q \
  "Data_driven/audit_uniform_data_input.cpp(\"${WEIGHT_DATA}\",\"${INPUT_DIR}/selection_audit.txt\")" \
  > "${RESULT}/logs/audit_uniform_data.log" 2>&1
root -l -b -q \
  "Data_driven/build_jpsi12_candidates.cpp+(\"${CANDIDATES}\",\"data\",\"${RAW_DATA_BASE}\")" \
  > "${RESULT}/logs/build_candidates.log" 2>&1

root -l -b -q \
  "Data_driven/fit_pp_2d_adaptive.cpp+(\"${PP_DIR}\",\"${WEIGHT_DATA}\",true,false,true,1.8,1.5707963267948966,\"${MODEL}\")" \
  > "${RESULT}/logs/fit_pp.log" 2>&1
grep -q '^accepted_fits=12$' "${PP_DIR}/summary.txt"

for slot in 1 2; do
  root -l -b -q \
    "Data_driven/fit_single_jpsi_global_calibrated_splot.cpp+(${slot},\"${CANDIDATES}\",\"${MODEL}\",\"${SPLOT_DIR}/jpsi${slot}_sweights.root\")" \
    > "${RESULT}/logs/fit_jpsi${slot}.log" 2>&1
  grep -q '^fit_mode=single_global_calibrated_two_stage$' \
    "${SPLOT_DIR}/global_calibrated_splot_summary_jpsi${slot}.txt"
done

root -l -b -q \
  "Data_driven/build_crossslot_splot_mixed.cpp+(0,0,\"${SPLOT_DIR}/jpsi1_sweights.root\",\"${SPLOT_DIR}/jpsi2_sweights.root\",\"${MIXED}\")" \
  > "${RESULT}/logs/build_mixing.log" 2>&1
grep -q '^mixing_mode=all_cross_event_pairs$' "${RESULT}/logs/build_mixing.log"

root -l -b -q \
  "Data_driven/build_2d_templates_adaptive.cpp+(\"${TEMPLATE_DIR}\",\"${PP_DIR}/fit_results.csv\",\"${MIXED}\",1.8,1.5707963267948966,0.0,true)" \
  > "${RESULT}/logs/build_templates.log" 2>&1
F_DPS=$(awk -F= '$1=="f_dps"{print $2}' "${TEMPLATE_DIR}/summary.txt")
test -n "${F_DPS}"

MODEL_SHA=$(sha256sum "${MODEL}" | awk '{print $1}')
WEIGHT_SHA=$(sha256sum "${WEIGHT_DATA}" | awk '{print $1}')
CANDIDATE_SHA=$(sha256sum "${CANDIDATES}" | awk '{print $1}')
{
  echo "status=candidate_complete_pending_user_confirmation"
  echo "run_tag=${TAG}"
  echo "method=global_calibrated_allpairs_singleCR_12cell_mergeddy"
  echo "f_dps=${F_DPS}"
  echo "root_version=$(root-config --version)"
  echo "git_head=$(git rev-parse HEAD)"
  echo "model=${MODEL}"
  echo "model_sha256=${MODEL_SHA}"
  echo "weight_data=${WEIGHT_DATA}"
  echo "weight_data_sha256=${WEIGHT_SHA}"
  echo "candidates=${CANDIDATES}"
  echo "candidates_sha256=${CANDIDATE_SHA}"
  echo "data_ntuple_base=${RAW_DATA_BASE}"
  echo "pp_fit_error_convention=ROOT_native_AsymptoticError_true"
  echo "selection=mJJ_ge_7p5_and_nominal_rephrase_selection"
  echo "generated_at=$(date --iso-8601=seconds)"
  sha256sum Data_driven/build_uniform_data_input.cpp \
    Data_driven/build_jpsi12_candidates.cpp \
    Data_driven/fit_pp_2d_adaptive.cpp \
    Data_driven/fit_single_jpsi_global_calibrated_splot.cpp \
    Data_driven/build_crossslot_splot_mixed.cpp \
    Data_driven/build_2d_templates_adaptive.cpp \
    Data_driven/run_nominal.sh
} > "${RESULT}/run.metadata.txt"
printf 'candidate_complete_pending_user_confirmation\n' > "${RESULT}/complete.marker"
echo "DATA_DRIVEN_NOMINAL_CANDIDATE_COMPLETE tag=${TAG} f_DPS=${F_DPS}"
