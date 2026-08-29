#!/bin/bash
set -euo pipefail

if [ "${GLOBAL_PROJECTION_INSIDE_EL7:-0}" != "1" ]; then
  exec /cvmfs/cms.cern.ch/common/cmssw-el7 --command-to-run \
    "export GLOBAL_PROJECTION_INSIDE_EL7=1; exec $0 ${1:-tests/atlas_data_driven_sps/results/fdps_systematics_global_nominal_v1/projections_1d}"
fi

REPO="/eos/home-l/leyao/26JJ/JPsiJPsi"
OUTPUT="${1:-${REPO}/tests/atlas_data_driven_sps/results/fdps_systematics_global_nominal_v1/projections_1d}"
TEMPLATES="${REPO}/tests/atlas_data_driven_sps/results/fdps_systematics_global_nominal_v1/baseline/templates_2d.root"
DATA="${REPO}/tests/atlas_data_driven_sps/results/unified_input/WeightData_mJJ7p5.root"
MIXED="${REPO}/tests/atlas_data_driven_sps/results/fdps_systematics_final_v1/category/global_calibrated/mixed_dps_allpairs.root"
TREES="${OUTPUT}/component_trees"
FITS="${OUTPUT}/component_fits"
TEMPLATE_FITS="${OUTPUT}/template_fraction_fits"
NOMINAL="0.170423232645"

mkdir -p "${OUTPUT}/logs" "${TREES}" "${FITS}/jobs" "${TEMPLATE_FITS}"
cd "${REPO}/Data/ULntuple16/CMSSW_10_6_20/src"
eval "$(scramv1 runtime -sh)"
cd "${REPO}"
test "$(root-config --version)" = "6.14/09"

root -l -b -q \
  "tests/atlas_data_driven_sps/build_component_trees.cpp+(\"${TREES}\",\"${TEMPLATES}\",\"${DATA}\",\"${MIXED}\")" \
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
xargs -a "${TASKS}" -n 5 -P 6 tests/atlas_data_driven_sps/run_component_fit_worker.sh
DONE=$(grep -l 'single_fit_result_written=' "${FITS}"/jobs/*/run.log | wc -l)
test "${DONE}" -eq 78

root -l -b -q \
  "tests/atlas_data_driven_sps/summarize_component_fits.cpp+(\"${FITS}\")" \
  > "${OUTPUT}/logs/summarize_component_fits.log" 2>&1
grep -q '^fits_accepted=78$' "${FITS}/summary.txt"

root -l -b -q \
  "tests/atlas_data_driven_sps/fit_1d_component_templates.cpp+(\"${TEMPLATE_FITS}\",\"${FITS}/fit_results.csv\")" \
  > "${OUTPUT}/logs/fit_1d_templates.log" 2>&1
ROWS=$(awk -F, 'NR>1 && $14==1 {n++} END{print n+0}' "${TEMPLATE_FITS}/template_fit_results.csv")
test "${ROWS}" -eq 4

{
  echo "status=complete"
  echo "scope=accepted_global_nominal_four_1d_projection_template_fit_crosscheck"
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
  echo "reason=ATLAS_uses_delta_y_for_reported_fDPS_and_other_distributions_as_crosschecks"
  echo "root_version=$(root-config --version)"
  echo "git_commit=$(git rev-parse HEAD)"
} > "${OUTPUT}/summary.txt"

{
  echo "artifact=accepted_global_nominal_1d_projection_crosscheck"
  echo "status=complete"
  echo "runtime=CMSSW_10_6_20_ROOT_$(root-config --version)_cmssw_el7"
  echo "two_dimensional_nominal=${NOMINAL}"
  echo "projection_binning=delta_y6_delta_phi8_evt_mass7_evt_y5"
  echo "signal_extraction=per_bin_4D_PP_fit_for_total_SPS_weighted_DPS_weighted_data"
  echo "template_fit=two_fixed_unit_normalized_SPS_DPS_shapes_to_total_PP_yields"
  echo "systematic_included=false"
  sha256sum \
    tests/atlas_data_driven_sps/build_component_trees.cpp \
    tests/atlas_data_driven_sps/fit_component_variables.cpp \
    tests/atlas_data_driven_sps/summarize_component_fits.cpp \
    tests/atlas_data_driven_sps/fit_1d_component_templates.cpp
} > "${OUTPUT}/metadata.txt"
printf 'complete\n' > "${OUTPUT}/complete.marker"
echo "GLOBAL_PROJECTION_CROSSCHECK_COMPLETE output=${OUTPUT}"
