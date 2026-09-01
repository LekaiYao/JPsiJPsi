#!/bin/bash
set -euo pipefail

if [ "$#" -lt 3 ] || [ "$#" -gt 5 ]; then
  echo "Usage: $0 <replica_id> <base_seed> <output_dir> [legacy14|nominal12_mergeddy] [absy1p2|global_calibrated]" >&2
  exit 64
fi
REPO="/eos/home-l/leyao/26JJ/JPsiJPsi"
source "${REPO}/Data_driven/require_root640.sh"
MODEL="${JJ_DD_MODEL:-${REPO}/Data_driven/inputs/Model_4D_tot.root}"
SOURCE_CANDIDATES="${JJ_DD_SOURCE_CANDIDATES:-${REPO}/Data_driven/results/input/jpsi12_candidates.root}"
SOURCE_DATA="${JJ_DD_SOURCE_DATA:-${REPO}/Data_driven/results/input/WeightData_mJJ7p5.root}"

REPLICA_ID="$1"; BASE_SEED="$2"; OUTPUT_DIR="$3"
BINNING_MODE="${4:-legacy14}"
SINGLE_J_MODE="${5:-absy1p2}"
case "${BINNING_MODE}" in
  legacy14)
    PP_FIT_CALL="Data_driven/fit_pp_2d_adaptive.cpp+(\"${OUTPUT_DIR}/pp_data_2d\",\"${OUTPUT_DIR}/input/WeightData_mJJ7p5.root\",false,false,false,1.8,1.5707963267948966,\"${MODEL}\")"
    EXPECTED_FITS=14
    PP_REFIT_LABEL="adaptive_14_cells_legacy"
    ;;
  nominal12_mergeddy)
    PP_FIT_CALL="Data_driven/fit_pp_2d_adaptive.cpp+(\"${OUTPUT_DIR}/pp_data_2d\",\"${OUTPUT_DIR}/input/WeightData_mJJ7p5.root\",true,false,true,1.8,1.5707963267948966,\"${MODEL}\")"
    EXPECTED_FITS=12
    PP_REFIT_LABEL="adaptive_12_cells_single_cr_high_dphi_columns_merged_across_abs_delta_y"
    ;;
  *)
    echo "Unsupported binning mode: ${BINNING_MODE}" >&2
    exit 64
    ;;
esac
case "${SINGLE_J_MODE}" in
  absy1p2|global_calibrated) ;;
  *) echo "Unsupported single-J mode: ${SINGLE_J_MODE}" >&2; exit 64 ;;
esac
INPUT_DIR="${OUTPUT_DIR}/input"; SPLOT_DIR="${OUTPUT_DIR}/splot"
PP_DIR="${OUTPUT_DIR}/pp_data_2d"; TEMPLATE_DIR="${OUTPUT_DIR}/templates_2d"
MIXED="${OUTPUT_DIR}/mixed_dps_allpairs.root"; TIMING="${OUTPUT_DIR}/timing.csv"

mkdir -p "${OUTPUT_DIR}" "${SPLOT_DIR}"
cd "${REPO}"
echo "stage,start_epoch,end_epoch,wall_seconds" > "${TIMING}"
TOTAL_START=$(date +%s)

run_stage() {
  local name="$1"; shift; local start end
  start=$(date +%s); "$@"; end=$(date +%s)
  echo "${name},${start},${end},$((end-start))" >> "${TIMING}"
}

run_stage resample root -l -b -q \
  "Data_driven/build_full_bootstrap_replica_input.cpp+(${REPLICA_ID},${BASE_SEED},\"${INPUT_DIR}\",\"${SOURCE_CANDIDATES}\",\"${SOURCE_DATA}\")"
for slot in 1 2; do
  if [ "${SINGLE_J_MODE}" = "global_calibrated" ]; then
    run_stage "global_calibrated_splot_jpsi${slot}" root -l -b -q \
      "Data_driven/fit_single_jpsi_global_calibrated_splot.cpp+(${slot},\"${INPUT_DIR}/jpsi12_candidates.root\",\"${MODEL}\",\"${SPLOT_DIR}/jpsi${slot}_sweights.root\")" \
      > "${OUTPUT_DIR}/fit_jpsi${slot}.log" 2>&1
    grep -q '^fit_mode=single_global_calibrated_two_stage$' \
      "${SPLOT_DIR}/global_calibrated_splot_summary_jpsi${slot}.txt"
  else
    run_stage "category_splot_jpsi${slot}" root -l -b -q \
      "Data_driven/fit_single_jpsi_category_splot.cpp+(${slot},\"absy\",1.2,\"${INPUT_DIR}/jpsi12_candidates.root\",\"${MODEL}\",\"${SPLOT_DIR}/jpsi${slot}_sweights.root\")" \
      > "${OUTPUT_DIR}/fit_jpsi${slot}.log" 2>&1
    grep -q '^calibration_mode=independent_factorized_categories$' \
      "${SPLOT_DIR}/category_splot_summary_jpsi${slot}.txt"
  fi
done
run_stage all_cross_event_mixing root -l -b -q \
  "Data_driven/build_crossslot_splot_mixed.cpp+(0,0,\"${SPLOT_DIR}/jpsi1_sweights.root\",\"${SPLOT_DIR}/jpsi2_sweights.root\",\"${MIXED}\")" \
  > "${OUTPUT_DIR}/build_mixing.log" 2>&1
grep -q '^mixing_mode=all_cross_event_pairs$' "${OUTPUT_DIR}/build_mixing.log"
run_stage "pp_${EXPECTED_FITS}cell_fit" root -l -b -q \
  "${PP_FIT_CALL}" \
  > "${OUTPUT_DIR}/fit_pp.log" 2>&1
grep -q "^accepted_fits=${EXPECTED_FITS}$" "${PP_DIR}/summary.txt"
run_stage nominal_cr_normalization root -l -b -q \
  "Data_driven/build_2d_templates_adaptive.cpp+(\"${TEMPLATE_DIR}\",\"${PP_DIR}/fit_results.csv\",\"${MIXED}\",1.8,1.5707963267948966,0.0,true)" \
  > "${OUTPUT_DIR}/build_templates.log" 2>&1
grep -q '^f_dps=' "${TEMPLATE_DIR}/summary.txt"

TOTAL_END=$(date +%s); echo "total,${TOTAL_START},${TOTAL_END},$((TOTAL_END-TOTAL_START))" >> "${TIMING}"
F_DPS=$(awk -F= '$1=="f_dps"{print $2}' "${TEMPLATE_DIR}/summary.txt")
{
  echo "scope=full_source_event_bootstrap_primary_2d_fdps"
  echo "replica_id=${REPLICA_ID}"
  echo "base_seed=${BASE_SEED}"
  echo "hostname=$(hostname)"
  echo "root_version=$(root-config --version)"
  echo "single_j_mode=${SINGLE_J_MODE}"
  if [ "${SINGLE_J_MODE}" = "global_calibrated" ]; then
    echo "category_variable=none"
    echo "category_split=none"
  else
    echo "category_variable=abs_y"
    echo "category_split=1.2"
  fi
  echo "mixing=all_cross_event_pairs"
  echo "binning_mode=${BINNING_MODE}"
  echo "pp_refit=${PP_REFIT_LABEL}"
  echo "pp_fits_expected=${EXPECTED_FITS}"
  echo "control_region=abs_delta_y_ge_1p8_and_abs_delta_phi_le_pi_over_2"
  echo "primary_f_dps=${F_DPS}"
  echo "component_projection_fits=false"
  echo "cr_systematic_variations=false"
  echo "wall_seconds=$((TOTAL_END-TOTAL_START))"
  echo "status=complete"
} > "${OUTPUT_DIR}/metadata.txt"
echo "FULL_BOOTSTRAP_REPLICA_COMPLETE replica=${REPLICA_ID} f_dps=${F_DPS} wall_seconds=$((TOTAL_END-TOTAL_START))"
