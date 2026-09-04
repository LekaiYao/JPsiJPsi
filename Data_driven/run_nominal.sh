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
INPUT_MANIFEST="${REPO}/Data_driven/inputs/input_manifest.txt"
ACCEPTANCE="${REPO}/Data_driven/inputs/acceptance_sps_full10_v1.txt"
EFFICIENCY="${REPO}/Data_driven/inputs/efficiency_sps0p8_dps0p2_dedup60_v1.txt"
SELECTION_SOURCE="${REPO}/Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/rephrase.cpp"
RAW_DATA_BASE="${JJ_DATA_NTUPLE_BASE:-/eos/user/c/chensh/JPsiJPsi/Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer}"
INPUT_DIR="${RESULT}/input"
DATA_INVENTORY="${INPUT_DIR}/data_ntuple_inventory.tsv"
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
DATA_FILES=0
{
  printf 'era\tindex\tpath\tsize_bytes\tmtime_epoch\tsha256\n'
  for spec in B:20 C:9 D:14 E:3 F:8 G:29 H:36; do
    era="${spec%:*}"
    last="${spec#*:}"
    for i in $(seq 1 "${last}"); do
      path="${RAW_DATA_BASE}/${era}/Ntuple_2016_${era}_${i}.root"
      test -s "${path}"
      size=$(stat -c %s "${path}")
      mtime=$(stat -c %Y "${path}")
      checksum=$(sha256sum "${path}" | awk '{print $1}')
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${era}" "${i}" "${path}" "${size}" "${mtime}" "${checksum}"
      DATA_FILES=$((DATA_FILES + 1))
    done
  done
} > "${DATA_INVENTORY}"
test "${DATA_FILES}" -eq 119
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
MANIFEST_SHA=$(sha256sum "${INPUT_MANIFEST}" | awk '{print $1}')
ACCEPTANCE_SHA=$(sha256sum "${ACCEPTANCE}" | awk '{print $1}')
EFFICIENCY_SHA=$(sha256sum "${EFFICIENCY}" | awk '{print $1}')
SELECTION_SHA=$(sha256sum "${SELECTION_SOURCE}" | awk '{print $1}')
INVENTORY_SHA=$(sha256sum "${DATA_INVENTORY}" | awk '{print $1}')
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
  echo "data_ntuple_files=${DATA_FILES}"
  echo "data_ntuple_inventory=${DATA_INVENTORY}"
  echo "data_ntuple_inventory_sha256=${INVENTORY_SHA}"
  echo "input_manifest=${INPUT_MANIFEST}"
  echo "input_manifest_sha256=${MANIFEST_SHA}"
  echo "acceptance=${ACCEPTANCE}"
  echo "acceptance_sha256=${ACCEPTANCE_SHA}"
  echo "efficiency=${EFFICIENCY}"
  echo "efficiency_sha256=${EFFICIENCY_SHA}"
  echo "selection_source=${SELECTION_SOURCE}"
  echo "selection_source_sha256=${SELECTION_SHA}"
  echo "fit_seed_policy=deterministic_multistart_no_random_seed"
  echo "mixing_seed_policy=exhaustive_all_cross_event_pairs_requested_0_seed_0"
  echo "metadata_binning_source=actual_cells_size"
  echo "pp_fit_error_convention=ROOT_native_AsymptoticError_true"
  echo "pp_comb_comb_boundary_strategy=fix_zero_and_refit_on_root_range_exception"
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
