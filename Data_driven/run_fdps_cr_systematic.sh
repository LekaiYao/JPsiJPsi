#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <confirmed-central-tag> <unique-cr-systematic-tag>" >&2
  exit 64
fi
CENTRAL_TAG="$1"
CR_TAG="$2"
for tag in "${CENTRAL_TAG}" "${CR_TAG}"; do
  [[ "${tag}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || {
    echo "Unsafe tag: ${tag}" >&2
    exit 64
  }
done

REPO="/eos/home-l/leyao/26JJ/JPsiJPsi"
source "${REPO}/Data_driven/require_root640.sh"
CENTRAL="${REPO}/Data_driven/results/${CENTRAL_TAG}"
RESULT="${REPO}/Data_driven/results/${CR_TAG}"
MODEL="${REPO}/Data_driven/inputs/Model_4D_tot.root"
DATA="${CENTRAL}/input/WeightData_mJJ7p5.root"
FIT_NOMINAL="${CENTRAL}/pp_data_2d/fit_results.csv"
GLOBAL="${CENTRAL}/mixed_dps_allpairs.root"
SPS="${REPO}/Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/WeightSPSstar.root"
PHI70="1.2217304763960306"
PHI90="1.5707963267948966"
PHI110="1.9198621771937625"

for required in "${MODEL}" "${DATA}" "${FIT_NOMINAL}" "${GLOBAL}" "${SPS}"   "${CENTRAL}/templates_2d/summary.txt" "${CENTRAL}/complete.marker"; do
  test -s "${required}" || {
    echo "Missing required input: ${required}" >&2
    exit 66
  }
done
if [ -e "${RESULT}" ]; then
  echo "Refusing to overwrite existing tag: ${RESULT}" >&2
  exit 68
fi
mkdir -p "${RESULT}/logs" "${RESULT}/baseline" "${RESULT}/cr"
cd "${REPO}"

NOMINAL=$(awk -F= '$1=="f_dps"{print $2}' "${CENTRAL}/templates_2d/summary.txt")
test -n "${NOMINAL}"
root -l -b -q   "Data_driven/build_2d_templates_adaptive.cpp+(\"${RESULT}/baseline\",\"${FIT_NOMINAL}\",\"${GLOBAL}\",1.8,${PHI90},0.0,true)"   > "${RESULT}/logs/baseline.log" 2>&1
BASELINE=$(awk -F= '$1=="f_dps"{print $2}' "${RESULT}/baseline/summary.txt")
awk -v a="${BASELINE}" -v b="${NOMINAL}"   'BEGIN{exit !(a-b<5e-10 && b-a<5e-10)}'

root -l -b -q   "Data_driven/measure_sps_control_fraction.cpp+(\"${RESULT}/cr/nominal_leakage.txt\",1.8,${PHI90},\"${SPS}\")"   > "${RESULT}/logs/nominal_leakage.log" 2>&1
NOMINAL_LEAKAGE=$(awk -F= '$1=="sps_control_fraction"{print $2}'   "${RESULT}/cr/nominal_leakage.txt")
root -l -b -q   "Data_driven/build_2d_templates_adaptive.cpp+(\"${RESULT}/cr/nominal_purity\",\"${FIT_NOMINAL}\",\"${GLOBAL}\",1.8,${PHI90},${NOMINAL_LEAKAGE},true)"   > "${RESULT}/logs/nominal_purity.log" 2>&1

run_cr() {
  local label="$1" dy="$2" phi="$3"
  local pp="${RESULT}/cr/${label}_pp"
  local cells leakage
  root -l -b -q     "Data_driven/fit_pp_2d_adaptive.cpp+(\"${pp}\",\"${DATA}\",true,false,true,${dy},${phi},\"${MODEL}\")"     > "${RESULT}/logs/${label}_pp.log" 2>&1
  cells=$(awk -F= '$1=="cells"{print $2}' "${pp}/summary.txt")
  grep -q "^accepted_fits=${cells}$" "${pp}/summary.txt"
  root -l -b -q     "Data_driven/measure_sps_control_fraction.cpp+(\"${RESULT}/cr/${label}_leakage.txt\",${dy},${phi},\"${SPS}\")"     > "${RESULT}/logs/${label}_leakage.log" 2>&1
  leakage=$(awk -F= '$1=="sps_control_fraction"{print $2}'     "${RESULT}/cr/${label}_leakage.txt")
  root -l -b -q     "Data_driven/build_2d_templates_adaptive.cpp+(\"${RESULT}/cr/${label}_zero\",\"${pp}/fit_results.csv\",\"${GLOBAL}\",${dy},${phi},0.0,true)"     > "${RESULT}/logs/${label}_zero.log" 2>&1
  root -l -b -q     "Data_driven/build_2d_templates_adaptive.cpp+(\"${RESULT}/cr/${label}_purity\",\"${pp}/fit_results.csv\",\"${GLOBAL}\",${dy},${phi},${leakage},true)"     > "${RESULT}/logs/${label}_purity.log" 2>&1
}

run_cr cr_phi70 1.8 "${PHI70}"
run_cr cr_phi110 1.8 "${PHI110}"
run_cr cr_dy2p4 2.4 "${PHI90}"

VARIATIONS="${RESULT}/variations.csv"
printf 'variation,note,f_dps,signed_shift,in_envelope\n' > "${VARIATIONS}"
write_row() {
  local label="$1" note="$2" path="$3" include="$4"
  local value shift
  value=$(awk -F= '$1=="f_dps"{print $2}' "${path}")
  test -n "${value}"
  shift=$(awk -v a="${value}" -v b="${NOMINAL}"     'BEGIN{printf "%.12g",a-b}')
  printf '%s,%s,%.12g,%.12g,%s\n'     "${label}" "${note}" "${value}" "${shift}" "${include}"     >> "${VARIATIONS}"
}
write_row nominal baseline "${RESULT}/baseline/summary.txt" 0
write_row nominal_purity mc_sps_leakage   "${RESULT}/cr/nominal_purity/summary.txt" 1
write_row cr_phi70_zero zero_sps_leakage   "${RESULT}/cr/cr_phi70_zero/summary.txt" 1
write_row cr_phi70_purity mc_sps_leakage   "${RESULT}/cr/cr_phi70_purity/summary.txt" 1
write_row cr_phi110_zero zero_sps_leakage   "${RESULT}/cr/cr_phi110_zero/summary.txt" 1
write_row cr_phi110_purity mc_sps_leakage   "${RESULT}/cr/cr_phi110_purity/summary.txt" 1
write_row cr_dy2p4_zero zero_sps_leakage   "${RESULT}/cr/cr_dy2p4_zero/summary.txt" 1
write_row cr_dy2p4_purity mc_sps_leakage   "${RESULT}/cr/cr_dy2p4_purity/summary.txt" 1

read -r CR_UP CR_DOWN < <(
  awk -F, 'NR>1 && $5==1 {
    shift=$4+0
    if (shift>up) up=shift
    if (-shift>down) down=-shift
  } END{printf "%.12g %.12g\n",up+0,down+0}' "${VARIATIONS}"
)
{
  echo "status=cr_systematic_complete_pending_user_confirmation"
  echo "central_tag=${CENTRAL_TAG}"
  echo "cr_systematic_tag=${CR_TAG}"
  echo "definition=envelope_across_nominal_purity_phi70_phi110_dy2p4_zero_and_mc_sps_leakage"
  echo "cr_group_up=${CR_UP}"
  echo "cr_group_down=${CR_DOWN}"
  awk -F, 'NR>1 {
    printf "%s_f_dps=%s\n",$1,$3
    printf "%s_signed_shift=%s\n",$1,$4
  }' "${VARIATIONS}"
  echo "root_version=$(root-config --version)"
  echo "git_head=$(git rev-parse HEAD)"
} > "${RESULT}/summary.txt"

{
  echo "artifact=fdps_cr_definition_and_purity_systematic"
  echo "status=complete_pending_user_confirmation"
  echo "central_tag=${CENTRAL_TAG}"
  echo "nominal_f_dps=${NOMINAL}"
  echo "cr_variations=nominal_purity,phi70_zero,phi70_purity,phi110_zero,phi110_purity,dy2p4_zero,dy2p4_purity"
  echo "pp_comb_comb_boundary_strategy=fix_zero_and_refit_on_root_range_exception"
  echo "runtime=ROOT_$(root-config --version)"
  echo "git_head=$(git rev-parse HEAD)"
  sha256sum "${CENTRAL}/run.metadata.txt" "${MODEL}" "${DATA}"     "${FIT_NOMINAL}" "${GLOBAL}" "${SPS}"     Data_driven/fit_pp_2d_adaptive.cpp     Data_driven/measure_sps_control_fraction.cpp     Data_driven/build_2d_templates_adaptive.cpp     Data_driven/run_fdps_cr_systematic.sh
} > "${RESULT}/run.metadata.txt"
printf 'cr_systematic_complete_pending_user_confirmation\n'   > "${RESULT}/complete.marker"
echo "DATA_DRIVEN_CR_SYSTEMATIC_COMPLETE tag=${CR_TAG} nominal_f_DPS=${NOMINAL} up=${CR_UP} down=${CR_DOWN}"
