#!/bin/bash
set -euo pipefail

OUTPUT_DIR="${1:-Data_driven/results/fdps_systematics_final_v1/dps_mc_full_workflow_closure_v1}"
SINGLE_J_MODE="${2:-absy1p2}"
REPO="/eos/home-l/leyao/26JJ/JPsiJPsi"
source "${REPO}/Data_driven/require_root640.sh"
MODEL="${REPO}/Data_driven/inputs/Model_4D_tot.root"
REFERENCE="${4:-${REPO}/Data_driven/results/dps_reference/WeightDPS_current.root}"
NOMINAL_FDPS="${3:-0.164184277213}"
case "${SINGLE_J_MODE}" in
  absy1p2|global_calibrated) ;;
  *) echo "Unsupported single-J mode: ${SINGLE_J_MODE}" >&2; exit 64 ;;
esac

mkdir -p "${OUTPUT_DIR}/splot"
cd "${REPO}"

root -l -b -q \
  "Data_driven/build_jpsi12_candidates.cpp+(\"${OUTPUT_DIR}/jpsi12_candidates.root\",\"dps_mc\")" \
  > "${OUTPUT_DIR}/build_candidates.log" 2>&1
for slot in 1 2; do
  if [ "${SINGLE_J_MODE}" = "global_calibrated" ]; then
    root -l -b -q \
      "Data_driven/fit_single_jpsi_global_calibrated_splot.cpp+(${slot},\"${OUTPUT_DIR}/jpsi12_candidates.root\",\"${MODEL}\",\"${OUTPUT_DIR}/splot/jpsi${slot}_sweights.root\")" \
      > "${OUTPUT_DIR}/fit_jpsi${slot}.log" 2>&1
    grep -q '^fit_mode=single_global_calibrated_two_stage$' \
      "${OUTPUT_DIR}/splot/global_calibrated_splot_summary_jpsi${slot}.txt"
  else
    root -l -b -q \
      "Data_driven/fit_single_jpsi_category_splot.cpp+(${slot},\"absy\",1.2,\"${OUTPUT_DIR}/jpsi12_candidates.root\",\"${MODEL}\",\"${OUTPUT_DIR}/splot/jpsi${slot}_sweights.root\")" \
      > "${OUTPUT_DIR}/fit_jpsi${slot}.log" 2>&1
    grep -q '^calibration_mode=independent_factorized_categories$' \
      "${OUTPUT_DIR}/splot/category_splot_summary_jpsi${slot}.txt"
  fi
done
root -l -b -q \
  "Data_driven/build_crossslot_splot_mixed.cpp+(0,0,\"${OUTPUT_DIR}/splot/jpsi1_sweights.root\",\"${OUTPUT_DIR}/splot/jpsi2_sweights.root\",\"${OUTPUT_DIR}/mixed_dps_allpairs.root\")" \
  > "${OUTPUT_DIR}/build_mixing.log" 2>&1
root -l -b -q \
  "Data_driven/fit_pp_2d_adaptive.cpp+(\"${OUTPUT_DIR}/pp_data_2d\",\"${REFERENCE}\",true,false,true,1.8,1.5707963267948966,\"${MODEL}\")" \
  > "${OUTPUT_DIR}/fit_pp.log" 2>&1
grep -q '^accepted_fits=12$' "${OUTPUT_DIR}/pp_data_2d/summary.txt"
root -l -b -q \
  "Data_driven/build_2d_templates_adaptive.cpp+(\"${OUTPUT_DIR}/templates_2d\",\"${OUTPUT_DIR}/pp_data_2d/fit_results.csv\",\"${OUTPUT_DIR}/mixed_dps_allpairs.root\",1.8,1.5707963267948966,0.0,true)" \
  > "${OUTPUT_DIR}/build_templates.log" 2>&1

F_DPS=$(awk -F= '$1=="f_dps"{print $2}' "${OUTPUT_DIR}/templates_2d/summary.txt")
ABS_SHIFT=$(awk -v f="${F_DPS}" 'BEGIN{d=f-1; if(d<0)d=-d; printf "%.12g",d}')
ABS_FDPS_SYST=$(awk -v f0="${NOMINAL_FDPS}" -v d="${ABS_SHIFT}" 'BEGIN{printf "%.12g",f0*d}')
{
  echo "artifact=pure_dps_mc_full_nominal_workflow_closure"
  echo "status=complete"
  echo "pseudo_data=pure_dps_mc_same_event_pairs"
  echo "single_jpsi_extraction=${SINGLE_J_MODE}_two_stage_splot"
  echo "mixing=nominal_deterministic_all_cross_event_Jpsi1_x_Jpsi2"
  echo "pp_extraction=nominal_12_cell_4D_fit"
  echo "control_region=abs_delta_y_ge_1p8_and_abs_delta_phi_le_pi_over_2"
  echo "expected_f_dps=1"
  echo "observed_f_dps=${F_DPS}"
  echo "absolute_closure_shift=${ABS_SHIFT}"
  echo "nominal_data_f_dps=${NOMINAL_FDPS}"
  echo "absolute_f_dps_systematic=${ABS_FDPS_SYST}"
  echo "root_version=$(root-config --version)"
  echo "git_commit=$(git rev-parse HEAD)"
} > "${OUTPUT_DIR}/metadata.txt"
touch "${OUTPUT_DIR}/complete.marker"
echo "DPS_MC_FULL_WORKFLOW_CLOSURE_COMPLETE f_DPS=${F_DPS} abs_shift=${ABS_SHIFT}"
