#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <confirmed-central-tag> <dps-reference-tag> <unique-systematics-tag>" >&2
  exit 64
fi
CENTRAL_TAG="$1"
DPS_REFERENCE_TAG="$2"
SYSTEMATICS_TAG="$3"
for tag in "${CENTRAL_TAG}" "${DPS_REFERENCE_TAG}" "${SYSTEMATICS_TAG}"; do
  [[ "${tag}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || {
    echo "Unsafe tag: ${tag}" >&2
    exit 64
  }
done

REPO="/eos/home-l/leyao/26JJ/JPsiJPsi"
source "${REPO}/Data_driven/require_root640.sh"
CENTRAL="${REPO}/Data_driven/results/${CENTRAL_TAG}"
RESULT="${REPO}/Data_driven/results/${SYSTEMATICS_TAG}"
MODEL="${REPO}/Data_driven/inputs/Model_4D_tot.root"
DATA="${CENTRAL}/input/WeightData_mJJ7p5.root"
CANDIDATES="${CENTRAL}/input/jpsi12_candidates.root"
FIT_NOMINAL="${CENTRAL}/pp_data_2d/fit_results.csv"
GLOBAL="${CENTRAL}/mixed_dps_allpairs.root"
DPS_REFERENCE="${REPO}/Data_driven/results/${DPS_REFERENCE_TAG}/WeightDPS_current.root"
SPS="${REPO}/Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/WeightSPSstar.root"
PHI70="1.2217304763960306"
PHI90="1.5707963267948966"
PHI110="1.9198621771937625"

for required in "${MODEL}" "${DATA}" "${CANDIDATES}" "${FIT_NOMINAL}" \
  "${GLOBAL}" "${DPS_REFERENCE}" "${SPS}"; do
  test -s "${required}" || { echo "Missing required input: ${required}" >&2; exit 66; }
done
test -f "${CENTRAL}/complete.marker" || {
  echo "Central tag has no completion marker: ${CENTRAL}" >&2
  exit 66
}
if [ -e "${RESULT}" ]; then
  echo "Refusing to overwrite existing tag: ${RESULT}" >&2
  exit 68
fi
mkdir -p "${RESULT}/logs" "${RESULT}/baseline" "${RESULT}/cr" \
  "${RESULT}/category/pt15/splot" "${RESULT}/category/absy1p2/splot"
cd "${REPO}"

NOMINAL=$(awk -F= '$1=="f_dps"{print $2}' "${CENTRAL}/templates_2d/summary.txt")
test -n "${NOMINAL}"
root -l -b -q \
  "Data_driven/build_2d_templates_adaptive.cpp+(\"${RESULT}/baseline\",\"${FIT_NOMINAL}\",\"${GLOBAL}\",1.8,${PHI90},0.0,true)" \
  > "${RESULT}/logs/baseline.log" 2>&1
BASELINE=$(awk -F= '$1=="f_dps"{print $2}' "${RESULT}/baseline/summary.txt")
awk -v a="${BASELINE}" -v b="${NOMINAL}" 'BEGIN{exit !(a-b<5e-10 && b-a<5e-10)}'

run_category() {
  local label="$1" variable="$2" split="$3"
  local output="${RESULT}/category/${label}"
  for slot in 1 2; do
    root -l -b -q \
      "Data_driven/fit_single_jpsi_category_splot.cpp+(${slot},\"${variable}\",${split},\"${CANDIDATES}\",\"${MODEL}\",\"${output}/splot/jpsi${slot}_sweights.root\")" \
      > "${RESULT}/logs/category_${label}_jpsi${slot}.log" 2>&1
    grep -q '^calibration_mode=independent_factorized_categories$' \
      "${output}/splot/category_splot_summary_jpsi${slot}.txt"
  done
  root -l -b -q \
    "Data_driven/build_crossslot_splot_mixed.cpp+(0,0,\"${output}/splot/jpsi1_sweights.root\",\"${output}/splot/jpsi2_sweights.root\",\"${output}/mixed_dps_allpairs.root\")" \
    > "${RESULT}/logs/category_${label}_mixing.log" 2>&1
  root -l -b -q \
    "Data_driven/build_2d_templates_adaptive.cpp+(\"${output}\",\"${FIT_NOMINAL}\",\"${output}/mixed_dps_allpairs.root\",1.8,${PHI90},0.0,true)" \
    > "${RESULT}/logs/category_${label}_templates.log" 2>&1
  grep -q '^f_dps=' "${output}/summary.txt"
}
run_category pt15 pt 15.0
run_category absy1p2 absy 1.2

root -l -b -q \
  "Data_driven/measure_sps_control_fraction.cpp+(\"${RESULT}/cr/nominal_leakage.txt\",1.8,${PHI90},\"${SPS}\")" \
  > "${RESULT}/logs/nominal_leakage.log" 2>&1
NOMINAL_LEAKAGE=$(awk -F= '$1=="sps_control_fraction"{print $2}' "${RESULT}/cr/nominal_leakage.txt")
root -l -b -q \
  "Data_driven/build_2d_templates_adaptive.cpp+(\"${RESULT}/cr/nominal_purity\",\"${FIT_NOMINAL}\",\"${GLOBAL}\",1.8,${PHI90},${NOMINAL_LEAKAGE},true)" \
  > "${RESULT}/logs/nominal_purity.log" 2>&1

run_cr() {
  local label="$1" dy="$2" phi="$3" pp="${RESULT}/cr/${1}_pp"
  root -l -b -q \
    "Data_driven/fit_pp_2d_adaptive.cpp+(\"${pp}\",\"${DATA}\",true,false,true,${dy},${phi},\"${MODEL}\")" \
    > "${RESULT}/logs/${label}_pp.log" 2>&1
  local cells leakage
  cells=$(awk -F= '$1=="cells"{print $2}' "${pp}/summary.txt")
  grep -q "^accepted_fits=${cells}$" "${pp}/summary.txt"
  root -l -b -q \
    "Data_driven/measure_sps_control_fraction.cpp+(\"${RESULT}/cr/${label}_leakage.txt\",${dy},${phi},\"${SPS}\")" \
    > "${RESULT}/logs/${label}_leakage.log" 2>&1
  leakage=$(awk -F= '$1=="sps_control_fraction"{print $2}' "${RESULT}/cr/${label}_leakage.txt")
  root -l -b -q \
    "Data_driven/build_2d_templates_adaptive.cpp+(\"${RESULT}/cr/${label}_zero\",\"${pp}/fit_results.csv\",\"${GLOBAL}\",${dy},${phi},0.0,true)" \
    > "${RESULT}/logs/${label}_zero.log" 2>&1
  root -l -b -q \
    "Data_driven/build_2d_templates_adaptive.cpp+(\"${RESULT}/cr/${label}_purity\",\"${pp}/fit_results.csv\",\"${GLOBAL}\",${dy},${phi},${leakage},true)" \
    > "${RESULT}/logs/${label}_purity.log" 2>&1
}
run_cr cr_phi70 1.8 "${PHI70}"
run_cr cr_phi110 1.8 "${PHI110}"
run_cr cr_dy2p4 2.4 "${PHI90}"

"${REPO}/Data_driven/run_dps_mc_full_workflow_closure.sh" \
  "${RESULT}/dps_mc_full_workflow_closure" global_calibrated "${NOMINAL}" "${DPS_REFERENCE}" \
  > "${RESULT}/logs/dps_mc_full_workflow_closure.log" 2>&1
grep -q '^status=complete$' "${RESULT}/dps_mc_full_workflow_closure/metadata.txt"

{
  echo "status=variations_complete_pending_bootstrap_projection_and_user_confirmation"
  echo "central_tag=${CENTRAL_TAG}"
  echo "dps_reference_tag=${DPS_REFERENCE_TAG}"
  echo "systematics_tag=${SYSTEMATICS_TAG}"
  echo "nominal_f_dps=${NOMINAL}"
  echo "root_version=$(root-config --version)"
  echo "git_head=$(git rev-parse HEAD)"
  sha256sum "${MODEL}" "${DATA}" "${CANDIDATES}" "${FIT_NOMINAL}" \
    "${GLOBAL}" "${DPS_REFERENCE}" "${SPS}" \
    Data_driven/run_fdps_global_nominal_systematics.sh
} > "${RESULT}/run.metadata.txt"
printf 'variations_complete_pending_bootstrap_projection_and_user_confirmation\n' \
  > "${RESULT}/complete.marker"
echo "DATA_DRIVEN_SYSTEMATIC_VARIATIONS_COMPLETE tag=${SYSTEMATICS_TAG} nominal_f_DPS=${NOMINAL}"
