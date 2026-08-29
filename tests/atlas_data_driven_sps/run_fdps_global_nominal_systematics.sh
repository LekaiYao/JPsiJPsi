#!/bin/bash
set -euo pipefail

if [ "${GLOBAL_SYSTEMATICS_INSIDE_EL7:-0}" != "1" ]; then
  exec /cvmfs/cms.cern.ch/common/cmssw-el7 --command-to-run \
    "export GLOBAL_SYSTEMATICS_INSIDE_EL7=1; exec $0 ${1:-tests/atlas_data_driven_sps/results/fdps_systematics_global_nominal_v1}"
fi

REPO="/eos/home-l/leyao/26JJ/JPsiJPsi"
RESULT="${1:-${REPO}/tests/atlas_data_driven_sps/results/fdps_systematics_global_nominal_v1}"
DATA="${REPO}/tests/atlas_data_driven_sps/results/unified_input/WeightData_mJJ7p5.root"
FIT_NOMINAL="${REPO}/tests/atlas_data_driven_sps/results/pp_data_2d_crsingle12_mergeddy/fit_results.csv"
GLOBAL="${REPO}/tests/atlas_data_driven_sps/results/fdps_systematics_final_v1/category/global_calibrated"
PT15="${REPO}/tests/atlas_data_driven_sps/results/fdps_systematics_final_v1/category/pt15"
ABSY="${REPO}/tests/atlas_data_driven_sps/results/fdps_systematics_final_v1/category/absy1p2_rerun"
SPS="${REPO}/Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/WeightSPSstar.root"
PHI70="1.2217304763960306"
PHI90="1.5707963267948966"
PHI110="1.9198621771937625"

mkdir -p "${RESULT}/logs" "${RESULT}/baseline" "${RESULT}/cr" \
  "${RESULT}/category/pt15" "${RESULT}/category/absy1p2"
cd "${REPO}/Data/ULntuple16/CMSSW_10_6_20/src"
eval "$(scramv1 runtime -sh)"
cd "${REPO}"
test "$(root-config --version)" = "6.14/09"

root -l -b -q \
  "tests/atlas_data_driven_sps/build_2d_templates_adaptive.cpp+(\"${RESULT}/baseline\",\"${FIT_NOMINAL}\",\"${GLOBAL}/mixed_dps_allpairs.root\",1.8,${PHI90},0.0,true)" \
  > "${RESULT}/logs/baseline.log" 2>&1
grep -q '^f_dps=0.170423232645$' "${RESULT}/baseline/summary.txt"

root -l -b -q \
  "tests/atlas_data_driven_sps/build_2d_templates_adaptive.cpp+(\"${RESULT}/category/pt15\",\"${FIT_NOMINAL}\",\"${PT15}/mixed_dps_allpairs.root\",1.8,${PHI90},0.0,true)" \
  > "${RESULT}/logs/category_pt15.log" 2>&1
root -l -b -q \
  "tests/atlas_data_driven_sps/build_2d_templates_adaptive.cpp+(\"${RESULT}/category/absy1p2\",\"${FIT_NOMINAL}\",\"${ABSY}/mixed_dps_allpairs.root\",1.8,${PHI90},0.0,true)" \
  > "${RESULT}/logs/category_absy1p2.log" 2>&1

root -l -b -q \
  "tests/atlas_data_driven_sps/measure_sps_control_fraction.cpp+(\"${RESULT}/cr/nominal_leakage.txt\",1.8,${PHI90},\"${SPS}\")" \
  > "${RESULT}/logs/nominal_leakage.log" 2>&1
NOMINAL_LEAKAGE=$(awk -F= '$1=="sps_control_fraction"{print $2}' "${RESULT}/cr/nominal_leakage.txt")
root -l -b -q \
  "tests/atlas_data_driven_sps/build_2d_templates_adaptive.cpp+(\"${RESULT}/cr/nominal_purity\",\"${FIT_NOMINAL}\",\"${GLOBAL}/mixed_dps_allpairs.root\",1.8,${PHI90},${NOMINAL_LEAKAGE},true)" \
  > "${RESULT}/logs/nominal_purity.log" 2>&1

run_cr() {
  local label="$1" dy="$2" phi="$3" pp="${RESULT}/cr/${1}_pp"
  root -l -b -q \
    "tests/atlas_data_driven_sps/fit_pp_2d_adaptive.cpp+(\"${pp}\",\"${DATA}\",true,false,true,${dy},${phi})" \
    > "${RESULT}/logs/${label}_pp.log" 2>&1
  local cells leakage
  cells=$(awk -F= '$1=="cells"{print $2}' "${pp}/summary.txt")
  grep -q "^accepted_fits=${cells}$" "${pp}/summary.txt"
  root -l -b -q \
    "tests/atlas_data_driven_sps/measure_sps_control_fraction.cpp+(\"${RESULT}/cr/${label}_leakage.txt\",${dy},${phi},\"${SPS}\")" \
    > "${RESULT}/logs/${label}_leakage.log" 2>&1
  leakage=$(awk -F= '$1=="sps_control_fraction"{print $2}' "${RESULT}/cr/${label}_leakage.txt")
  root -l -b -q \
    "tests/atlas_data_driven_sps/build_2d_templates_adaptive.cpp+(\"${RESULT}/cr/${label}_zero\",\"${pp}/fit_results.csv\",\"${GLOBAL}/mixed_dps_allpairs.root\",${dy},${phi},0.0,true)" \
    > "${RESULT}/logs/${label}_zero.log" 2>&1
  root -l -b -q \
    "tests/atlas_data_driven_sps/build_2d_templates_adaptive.cpp+(\"${RESULT}/cr/${label}_purity\",\"${pp}/fit_results.csv\",\"${GLOBAL}/mixed_dps_allpairs.root\",${dy},${phi},${leakage},true)" \
    > "${RESULT}/logs/${label}_purity.log" 2>&1
}
run_cr cr_phi70 1.8 "${PHI70}"
run_cr cr_phi110 1.8 "${PHI110}"
run_cr cr_dy2p4 2.4 "${PHI90}"

"${REPO}/tests/atlas_data_driven_sps/run_dps_mc_full_workflow_closure.sh" \
  "${RESULT}/dps_mc_full_workflow_closure" global_calibrated 0.170423232645 \
  > "${RESULT}/logs/dps_mc_full_workflow_closure.log" 2>&1
grep -q '^status=complete$' "${RESULT}/dps_mc_full_workflow_closure/metadata.txt"
python3 tests/atlas_data_driven_sps/summarize_fdps_global_nominal_systematics.py "${RESULT}" \
  > "${RESULT}/logs/summarize.log" 2>&1
printf 'complete\n' > "${RESULT}/complete.marker"
echo "GLOBAL_NOMINAL_SYSTEMATICS_COMPLETE result=${RESULT}"
