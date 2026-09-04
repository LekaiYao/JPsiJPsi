#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <confirmed-central-tag> <systematics-tag>" >&2
  exit 64
fi
CENTRAL_TAG="$1"
SYSTEMATICS_TAG="$2"
for tag in "${CENTRAL_TAG}" "${SYSTEMATICS_TAG}"; do
  [[ "${tag}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || {
    echo "Unsafe tag: ${tag}" >&2
    exit 64
  }
done

REPO="/eos/home-l/leyao/26JJ/JPsiJPsi"
source "${REPO}/Data_driven/require_root640.sh"
CENTRAL="${REPO}/Data_driven/results/${CENTRAL_TAG}"
SYSTEMATICS="${REPO}/Data_driven/results/${SYSTEMATICS_TAG}"
OUTPUT="${SYSTEMATICS}/projections_1d"
TEMPLATES="${SYSTEMATICS}/baseline/templates_2d.root"
DATA="${CENTRAL}/input/WeightData_mJJ7p5.root"
MIXED="${CENTRAL}/mixed_dps_allpairs.root"
TREES="${OUTPUT}/component_trees"
FITS="${OUTPUT}/component_fits"
TEMPLATE_FITS="${OUTPUT}/template_fraction_fits"

for required in "${TEMPLATES}" "${DATA}" "${MIXED}" \
  "${CENTRAL}/templates_2d/summary.txt"; do
  test -s "${required}" || { echo "Missing required input: ${required}" >&2; exit 66; }
done
if [ -e "${OUTPUT}" ]; then
  echo "Refusing to overwrite existing projection output: ${OUTPUT}" >&2
  exit 68
fi
NOMINAL=$(awk -F= '$1=="f_dps"{print $2}' "${CENTRAL}/templates_2d/summary.txt")
test -n "${NOMINAL}"

mkdir -p "${OUTPUT}/logs" "${TREES}" "${FITS}/jobs" "${TEMPLATE_FITS}"
cd "${REPO}"

root -l -b -q \
  "Data_driven/build_component_trees.cpp+(\"${TREES}\",\"${TEMPLATES}\",\"${DATA}\",\"${MIXED}\")" \
  > "${OUTPUT}/logs/build_component_trees.log" 2>&1
grep -q '^data_unassigned_entries=0$' "${TREES}/summary.txt"
grep -q '^mix_unassigned_entries=0$' "${TREES}/summary.txt"

TASKS="${FITS}/tasks.txt"
: > "${TASKS}"
for spec in delta_y:6 delta_phi:8 evt_mass:7 evt_y:5; do
  IFS=: read -r variable bins <<< "${spec}"
  for ((bin=0; bin<bins; ++bin)); do
    for component in total sps dps; do
      printf '%s %s %s %s %s\n' "${variable}" "${bin}" "${component}" \
        "${FITS}" "${TREES}/component_data.root" >> "${TASKS}"
    done
  done
done
xargs -a "${TASKS}" -n 5 -P 6 Data_driven/run_component_fit_worker.sh
DONE=$(grep -l 'single_fit_result_written=' "${FITS}"/jobs/*/run.log | wc -l)
test "${DONE}" -eq 78

root -l -b -q \
  "Data_driven/summarize_component_fits.cpp+(\"${FITS}\")" \
  > "${OUTPUT}/logs/summarize_component_fits.log" 2>&1
grep -q '^fits_accepted=78$' "${FITS}/summary.txt"

root -l -b -q \
  "Data_driven/fit_1d_component_templates.cpp+(\"${TEMPLATE_FITS}\",\"${FITS}/fit_results.csv\")" \
  > "${OUTPUT}/logs/fit_1d_templates.log" 2>&1
ROWS=$(awk -F, 'NR>1 && $14==1 {n++} END{print n+0}' "${TEMPLATE_FITS}/template_fit_results.csv")
test "${ROWS}" -eq 4

{
  echo "status=complete_pending_user_confirmation"
  echo "central_tag=${CENTRAL_TAG}"
  echo "systematics_tag=${SYSTEMATICS_TAG}"
  echo "nominal_f_dps=${NOMINAL}"
  awk -F, -v nominal="${NOMINAL}" 'NR>1 {
    printf "%s_f_dps=%.12g\n",$1,$10;
    printf "%s_fit_error=%.12g\n",$1,$9;
    printf "%s_signed_shift=%.12g\n",$1,$10-nominal;
    printf "%s_chi2_ndf=%.12g\n",$1,$13;
  }' "${TEMPLATE_FITS}/template_fit_results.csv"
  echo "variables=delta_y,delta_phi,evt_mass,evt_y"
  echo "component_fits=78"
  echo "component_fits_accepted=78"
  echo "projection_role=formal_crosscheck"
  echo "included_in_systematic=false"
  echo "root_version=$(root-config --version)"
  echo "git_head=$(git rev-parse HEAD)"
} > "${OUTPUT}/summary.txt"

{
  echo "artifact=global_nominal_1d_projection_crosscheck_candidate"
  echo "status=complete_pending_user_confirmation"
  echo "runtime=ROOT_$(root-config --version)"
  echo "two_dimensional_nominal=${NOMINAL}"
  echo "projection_binning=delta_y6_delta_phi8_evt_mass7_evt_y5"
  echo "signal_extraction=per_bin_4D_PP_fit_for_total_SPS_weighted_DPS_weighted_data"
  echo "template_fit=two_fixed_unit_normalized_SPS_DPS_shapes_to_total_PP_yields"
  echo "component_worker_execution=parallel_interpreted_macro_no_aclic"
  echo "parallel_workers=6"
  echo "systematic_included=false"
  sha256sum Data_driven/run_global_nominal_projection_crosscheck.sh \
    Data_driven/run_component_fit_worker.sh \
    Data_driven/build_component_trees.cpp \
    Data_driven/fit_component_variables.cpp \
    Data_driven/summarize_component_fits.cpp \
    Data_driven/fit_1d_component_templates.cpp
} > "${OUTPUT}/metadata.txt"
printf 'complete_pending_user_confirmation\n' > "${OUTPUT}/complete.marker"
echo "DATA_DRIVEN_PROJECTION_CROSSCHECK_COMPLETE output=${OUTPUT}"
