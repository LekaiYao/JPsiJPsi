#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 5 ]; then
  echo "Usage: $0 <central-pp-tag> <full-dpsmc-tag> <dps-reference-tag> <cr-pp-source-tag> <unique-output-tag>" >&2
  exit 64
fi

CENTRAL_TAG="$1"
FULL_DPSMC_TAG="$2"
DPS_REFERENCE_TAG="$3"
CR_PP_SOURCE_TAG="$4"
OUTPUT_TAG="$5"
for tag in "${CENTRAL_TAG}" "${FULL_DPSMC_TAG}" "${DPS_REFERENCE_TAG}" \
  "${CR_PP_SOURCE_TAG}" "${OUTPUT_TAG}"; do
  [[ "${tag}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || {
    echo "Unsafe tag: ${tag}" >&2
    exit 64
  }
done

REPO="/eos/home-l/leyao/26JJ/JPsiJPsi"
source "${REPO}/Data_driven/require_root640.sh"

SEED=20260902
PHI70="1.2217304763960306"
PHI90="1.5707963267948966"
PHI110="1.9198621771937625"
CENTRAL="${REPO}/Data_driven/results/${CENTRAL_TAG}"
FULL_DPSMC="${REPO}/Data_driven/results/${FULL_DPSMC_TAG}"
DPS_REFERENCE="${REPO}/Data_driven/results/${DPS_REFERENCE_TAG}"
CR_PP_SOURCE="${REPO}/Data_driven/results/${CR_PP_SOURCE_TAG}"
OUTPUT="${REPO}/Data_driven/results/${OUTPUT_TAG}"
DPS_MC_INPUT="${DPS_REFERENCE}/nominal_dps_template_input.root"
NOMINAL_FIT="${CENTRAL}/pp_data_2d/fit_results.csv"
NOMINAL_PP_SUMMARY="${CENTRAL}/pp_data_2d/summary.txt"
HALF_DIR="${OUTPUT}/half_sample_seed${SEED}"
HALF_INPUT="${HALF_DIR}/dps_mc_half.root"
HALF_SELECTION_SUMMARY="${HALF_DIR}/sampling_summary.txt"
HALF_SELECTION_CSV="${HALF_DIR}/selected_entries.csv"
HALF_TEMPLATE="${HALF_DIR}/templates_2d"
CR_DIR="${OUTPUT}/cr_variations_full_mc"
INPUT_MANIFEST="${REPO}/Data_driven/inputs/input_manifest.txt"
ACCEPTANCE="${REPO}/Data_driven/inputs/acceptance_sps_full10_v1.txt"
EFFICIENCY="${REPO}/Data_driven/inputs/efficiency_sps0p8_dps0p2_dedup60_v1.txt"

for required in \
  "${CENTRAL}/complete.marker" "${CENTRAL}/run.metadata.txt" \
  "${NOMINAL_FIT}" "${NOMINAL_PP_SUMMARY}" \
  "${FULL_DPSMC}/complete.marker" "${FULL_DPSMC}/summary.txt" \
  "${FULL_DPSMC}/run.metadata.txt" "${FULL_DPSMC}/artifact_checksums.txt" \
  "${DPS_REFERENCE}/complete.marker" "${DPS_REFERENCE}/summary.txt" \
  "${DPS_REFERENCE}/run.metadata.txt" "${DPS_REFERENCE}/artifact_checksums.txt" \
  "${DPS_MC_INPUT}" \
  "${CR_PP_SOURCE}/complete.marker" "${CR_PP_SOURCE}/summary.txt" \
  "${CR_PP_SOURCE}/run.metadata.txt" \
  "${CR_PP_SOURCE}/cr/cr_phi70_pp/fit_results.csv" \
  "${CR_PP_SOURCE}/cr/cr_phi70_pp/summary.txt" \
  "${CR_PP_SOURCE}/cr/cr_phi110_pp/fit_results.csv" \
  "${CR_PP_SOURCE}/cr/cr_phi110_pp/summary.txt" \
  "${CR_PP_SOURCE}/cr/cr_dy2p4_pp/fit_results.csv" \
  "${CR_PP_SOURCE}/cr/cr_dy2p4_pp/summary.txt" \
  "${INPUT_MANIFEST}" "${ACCEPTANCE}" "${EFFICIENCY}"; do
  test -s "${required}" || {
    echo "Missing required input: ${required}" >&2
    exit 66
  }
done

if [ -e "${OUTPUT}" ]; then
  echo "Refusing to overwrite existing tag: ${OUTPUT}" >&2
  exit 68
fi

(
  cd "${FULL_DPSMC}"
  sha256sum -c artifact_checksums.txt >/dev/null
)
(
  cd "${DPS_REFERENCE}"
  sha256sum -c artifact_checksums.txt >/dev/null
)

FULL_F_DPS=$(awk -F= '$1=="f_dps"{print $2}' "${FULL_DPSMC}/summary.txt")
FULL_F_SPS=$(awk -F= '$1=="f_sps"{print $2}' "${FULL_DPSMC}/summary.txt")
test -n "${FULL_F_DPS}" && test -n "${FULL_F_SPS}"
test "$(awk -F= '$1=="dps_reference_tag"{print $2}' \
  "${FULL_DPSMC}/run.metadata.txt")" = "${DPS_REFERENCE_TAG}"
test "$(awk -F= '$1=="dps_template_source"{print $2}' \
  "${FULL_DPSMC}/templates_2d/summary.txt")" = "dps_reconstruction_mc"

ACTUAL_ACCEPTANCE_SHA=$(sha256sum "${ACCEPTANCE}" | awk '{print $1}')
ACTUAL_EFFICIENCY_SHA=$(sha256sum "${EFFICIENCY}" | awk '{print $1}')
test "${ACTUAL_ACCEPTANCE_SHA}" = "$(awk -F= '$1=="acceptance_sha256"{print $2}' \
  "${FULL_DPSMC}/run.metadata.txt")"
test "${ACTUAL_EFFICIENCY_SHA}" = "$(awk -F= '$1=="efficiency_sha256"{print $2}' \
  "${FULL_DPSMC}/run.metadata.txt")"

validate_pp() {
  local summary="$1"
  local cells accepted failed
  cells=$(awk -F= '$1=="cells"{print $2}' "${summary}")
  accepted=$(awk -F= '$1=="accepted_fits"{print $2}' "${summary}")
  failed=$(awk -F= '$1=="failed_fits"{print $2}' "${summary}")
  test -n "${cells}" && test "${accepted}" -eq "${cells}" && test "${failed}" -eq 0
}
validate_pp "${NOMINAL_PP_SUMMARY}"
validate_pp "${CR_PP_SOURCE}/cr/cr_phi70_pp/summary.txt"
validate_pp "${CR_PP_SOURCE}/cr/cr_phi110_pp/summary.txt"
validate_pp "${CR_PP_SOURCE}/cr/cr_dy2p4_pp/summary.txt"

mkdir -p "${HALF_TEMPLATE}" "${CR_DIR}" "${OUTPUT}/logs"
cd "${REPO}"

root -l -b -q \
  "Data_driven/subsample_dps_template_fixed.cpp+(\"${DPS_MC_INPUT}\",\"${HALF_INPUT}\",\"${HALF_SELECTION_SUMMARY}\",\"${HALF_SELECTION_CSV}\",${SEED})" \
  > "${OUTPUT}/logs/build_half_sample.log" 2>&1
grep -q '^DPS_FIXED_HALF_SAMPLE_COMPLETE ' "${OUTPUT}/logs/build_half_sample.log"
test "$(awk -F= '$1=="input_entries"{print $2}' "${HALF_SELECTION_SUMMARY}")" -eq 5168
test "$(awk -F= '$1=="selected_entries"{print $2}' "${HALF_SELECTION_SUMMARY}")" -eq 2584
test "$(awk -F= '$1=="seed"{print $2}' "${HALF_SELECTION_SUMMARY}")" -eq "${SEED}"
test "$(awk -F= '$1=="priority_collisions"{print $2}' \
  "${HALF_SELECTION_SUMMARY}")" -eq 0

root -l -b -q \
  "Data_driven/build_2d_templates_adaptive.cpp+(\"${HALF_TEMPLATE}\",\"${NOMINAL_FIT}\",\"${HALF_INPUT}\",1.8,${PHI90},0.0,false,\"dps_reconstruction_mc\")" \
  > "${OUTPUT}/logs/half_nominal.log" 2>&1
test "$(awk -F= '$1=="mixed_pairs"{print $2}' "${HALF_TEMPLATE}/summary.txt")" \
  -eq 2584
test "$(awk -F= '$1=="dps_template_source"{print $2}' \
  "${HALF_TEMPLATE}/summary.txt")" = "dps_reconstruction_mc"

HALF_F_DPS=$(awk -F= '$1=="f_dps"{print $2}' "${HALF_TEMPLATE}/summary.txt")
HALF_F_SPS=$(awk -F= '$1=="f_sps"{print $2}' "${HALF_TEMPLATE}/summary.txt")
HALF_SIGNED_SHIFT=$(awk -v h="${HALF_F_DPS}" -v f="${FULL_F_DPS}" \
  'BEGIN{printf "%+.12g",h-f}')
STAT_ABS=$(awk -v h="${HALF_F_DPS}" -v f="${FULL_F_DPS}" \
  'BEGIN{d=h-f; if(d<0)d=-d; printf "%.12g",d}')

VARIATIONS="${OUTPUT}/cr_variations.csv"
printf 'variation,dy_min,phi_max,cells,f_dps,signed_shift,negative_sps_bins,pp_accepted\n' \
  > "${VARIATIONS}"
printf 'nominal,1.8,%.12g,12,%.12g,0,%s,12\n' "${PHI90}" "${FULL_F_DPS}" \
  "$(awk -F= '$1=="negative_sps_bins"{print $2}' "${FULL_DPSMC}/summary.txt")" \
  >> "${VARIATIONS}"

run_cr() {
  local label="$1"
  local dy="$2"
  local phi="$3"
  local pp_dir="$4"
  local result_dir="${CR_DIR}/${label}"
  local cells accepted value shift negative
  cells=$(awk -F= '$1=="cells"{print $2}' "${pp_dir}/summary.txt")
  accepted=$(awk -F= '$1=="accepted_fits"{print $2}' "${pp_dir}/summary.txt")
  root -l -b -q \
    "Data_driven/build_2d_templates_adaptive.cpp+(\"${result_dir}\",\"${pp_dir}/fit_results.csv\",\"${DPS_MC_INPUT}\",${dy},${phi},0.0,false,\"dps_reconstruction_mc\")" \
    > "${OUTPUT}/logs/${label}.log" 2>&1
  test "$(awk -F= '$1=="cells"{print $2}' "${result_dir}/summary.txt")" \
    -eq "${cells}"
  test "$(awk -F= '$1=="dps_template_source"{print $2}' \
    "${result_dir}/summary.txt")" = "dps_reconstruction_mc"
  value=$(awk -F= '$1=="f_dps"{print $2}' "${result_dir}/summary.txt")
  shift=$(awk -v v="${value}" -v f="${FULL_F_DPS}" \
    'BEGIN{printf "%+.12g",v-f}')
  negative=$(awk -F= '$1=="negative_sps_bins"{print $2}' \
    "${result_dir}/summary.txt")
  printf '%s,%.12g,%.12g,%s,%.12g,%.12g,%s,%s\n' \
    "${label}" "${dy}" "${phi}" "${cells}" "${value}" "${shift}" \
    "${negative}" "${accepted}" >> "${VARIATIONS}"
}

run_cr phi70 1.8 "${PHI70}" "${CR_PP_SOURCE}/cr/cr_phi70_pp"
run_cr phi110 1.8 "${PHI110}" "${CR_PP_SOURCE}/cr/cr_phi110_pp"
run_cr dy2p4 2.4 "${PHI90}" "${CR_PP_SOURCE}/cr/cr_dy2p4_pp"

read -r CR_UP CR_DOWN CR_SYMMETRIC < <(
  awk -F, 'NR>2 {
    shift=$6+0
    if (shift>up) up=shift
    if (-shift>down) down=-shift
  } END {
    symmetric=up>down ? up : down
    printf "%.12g %.12g %.12g\n",up+0,down+0,symmetric+0
  }' "${VARIATIONS}"
)

HALF_SELECTED_WEIGHT_FRACTION=$(awk -F= \
  '$1=="selected_weight_fraction"{print $2}' "${HALF_SELECTION_SUMMARY}")
HALF_SELECTED_CR_FRACTION=$(awk -F= \
  '$1=="selected_nominal_cr_fraction"{print $2}' "${HALF_SELECTION_SUMMARY}")
FULL_MC_CR_FRACTION=$(awk -F= '$1=="dps_mc_control_fraction"{print $2}' \
  "${FULL_DPSMC}/summary.txt")
TOTAL_PP_FITS=0
TOTAL_FIXED_ZERO=0
MAX_EDM=0
MIN_COVQUAL=999
for summary in "${NOMINAL_PP_SUMMARY}" \
  "${CR_PP_SOURCE}/cr/cr_phi70_pp/summary.txt" \
  "${CR_PP_SOURCE}/cr/cr_phi110_pp/summary.txt" \
  "${CR_PP_SOURCE}/cr/cr_dy2p4_pp/summary.txt"; do
  cells=$(awk -F= '$1=="cells"{print $2}' "${summary}")
  TOTAL_PP_FITS=$((TOTAL_PP_FITS + cells))
  fixed=$(awk -F= '$1=="comb_comb_fixed_zero_fits"{print $2}' "${summary}")
  TOTAL_FIXED_ZERO=$((TOTAL_FIXED_ZERO + fixed))
  csv="${summary%/summary.txt}/fit_results.csv"
  this_min=$(awk -F, 'NR==2{m=$15} NR>1 && $15<m{m=$15} END{print m}' "${csv}")
  this_max=$(awk -F, 'NR==2{m=$16} NR>1 && $16>m{m=$16} END{print m}' "${csv}")
  MIN_COVQUAL=$(awk -v a="${MIN_COVQUAL}" -v b="${this_min}" \
    'BEGIN{print (a<b ? a : b)}')
  MAX_EDM=$(awk -v a="${MAX_EDM}" -v b="${this_max}" \
    'BEGIN{print (a>b ? a : b)}')
done

{
  echo "status=candidate_complete_pending_user_confirmation"
  echo "artifact=fdps_dpsmc_fixed_half_stat_and_cr_systematic"
  echo "method=dps_reconstruction_mc_template_CR_normalized"
  echo "f_dps=${FULL_F_DPS}"
  echo "f_sps=${FULL_F_SPS}"
  echo "stat_definition=absolute_shift_fixed_seed_exact_half_DPS_MC_vs_full_DPS_MC"
  echo "stat_seed=${SEED}"
  echo "stat_engine=std_mt19937_64"
  echo "stat_sampling=exact_half_without_replacement"
  echo "stat_input_entries=5168"
  echo "stat_selected_entries=2584"
  echo "stat_selected_weight_fraction=${HALF_SELECTED_WEIGHT_FRACTION}"
  echo "stat_full_mc_cr_fraction=${FULL_MC_CR_FRACTION}"
  echo "stat_half_mc_cr_fraction=${HALF_SELECTED_CR_FRACTION}"
  echo "stat_half_f_dps=${HALF_F_DPS}"
  echo "stat_half_f_sps=${HALF_F_SPS}"
  echo "stat_signed_shift=${HALF_SIGNED_SHIFT}"
  echo "stat_uncertainty=${STAT_ABS}"
  echo "stat_scope=user_defined_single_fixed_half_sample_DPS_MC_template_shift_only"
  echo "cr_systematic_definition=zero_leakage_phi70_phi110_dy2p4_envelope_using_full_DPS_MC"
  echo "cr_systematic_up=${CR_UP}"
  echo "cr_systematic_down=${CR_DOWN}"
  echo "cr_systematic_symmetric_max=${CR_SYMMETRIC}"
  awk -F, 'NR>2 {
    printf "cr_%s_f_dps=%s\n",$1,$5
    printf "cr_%s_signed_shift=%s\n",$1,$6
    printf "cr_%s_negative_sps_bins=%s\n",$1,$7
  }' "${VARIATIONS}"
  echo "sps_leakage_fraction=0_for_central_and_all_CR_variations"
  echo "pp_fits_reused=${TOTAL_PP_FITS}"
  echo "pp_fit_failed=0"
  echo "pp_fit_min_covQual=${MIN_COVQUAL}"
  echo "pp_fit_max_edm=${MAX_EDM}"
  echo "pp_fit_comb_comb_fixed_zero=${TOTAL_FIXED_ZERO}"
  echo "new_roofit_fits_run=0"
  echo "correction_scope=frozen_current_Data_driven_inputs_not_next_nominal_full_SPS_acceptance"
} > "${OUTPUT}/summary.txt"

MODEL_SHA=$(awk -F= '$1=="model_sha256"{print $2}' \
  "${CENTRAL}/run.metadata.txt")
MANIFEST_SHA=$(sha256sum "${INPUT_MANIFEST}" | awk '{print $1}')
DPS_MC_SHA=$(sha256sum "${DPS_MC_INPUT}" | awk '{print $1}')
HALF_INPUT_SHA=$(sha256sum "${HALF_INPUT}" | awk '{print $1}')
HALF_SELECTION_SHA=$(sha256sum "${HALF_SELECTION_CSV}" | awk '{print $1}')
{
  echo "status=candidate_complete_pending_user_confirmation"
  echo "output_tag=${OUTPUT_TAG}"
  echo "central_pp_tag=${CENTRAL_TAG}"
  echo "full_dpsmc_tag=${FULL_DPSMC_TAG}"
  echo "dps_reference_tag=${DPS_REFERENCE_TAG}"
  echo "cr_pp_source_tag=${CR_PP_SOURCE_TAG}"
  echo "root_version=$(root-config --version)"
  echo "git_head=$(git rev-parse HEAD)"
  echo "model_sha256=${MODEL_SHA}"
  echo "input_manifest_sha256=${MANIFEST_SHA}"
  echo "acceptance_sha256=${ACTUAL_ACCEPTANCE_SHA}"
  echo "efficiency_sha256=${ACTUAL_EFFICIENCY_SHA}"
  echo "dps_mc_input=${DPS_MC_INPUT}"
  echo "dps_mc_input_sha256=${DPS_MC_SHA}"
  echo "half_sample_input=${HALF_INPUT}"
  echo "half_sample_input_sha256=${HALF_INPUT_SHA}"
  echo "half_sample_selected_entries_sha256=${HALF_SELECTION_SHA}"
  echo "stat_seed=${SEED}"
  echo "stat_engine=std_mt19937_64"
  echo "stat_sampling=exact_half_without_replacement"
  echo "cr_variations=phi70_phi110_dy2p4_zero_SPS_leakage"
  echo "cr_template=full_DPS_reconstruction_MC"
  echo "pp_fit_policy=reuse_accepted_existing_fits_no_refit"
  echo "correction_scope=frozen_current_Data_driven_inputs_not_next_nominal_full_SPS_acceptance"
  echo "generated_at=$(date --iso-8601=seconds)"
  sha256sum "${CENTRAL}/run.metadata.txt" "${NOMINAL_FIT}" \
    "${FULL_DPSMC}/run.metadata.txt" \
    "${DPS_REFERENCE}/run.metadata.txt" \
    "${CR_PP_SOURCE}/run.metadata.txt" \
    "${CR_PP_SOURCE}/cr/cr_phi70_pp/fit_results.csv" \
    "${CR_PP_SOURCE}/cr/cr_phi110_pp/fit_results.csv" \
    "${CR_PP_SOURCE}/cr/cr_dy2p4_pp/fit_results.csv" \
    Data_driven/subsample_dps_template_fixed.cpp \
    Data_driven/build_2d_templates_adaptive.cpp \
    Data_driven/run_fdps_dpsmc_halfstat_crsystematic.sh
} > "${OUTPUT}/run.metadata.txt"

find "${OUTPUT}" -maxdepth 4 -type f ! -name artifact_checksums.txt -print0 \
  | sort -z | xargs -0 sha256sum > "${OUTPUT}/artifact_checksums.txt"
printf 'candidate_complete_pending_user_confirmation\n' > "${OUTPUT}/complete.marker"
echo "DATA_DRIVEN_DPSMC_HALFSTAT_CRSYST_COMPLETE tag=${OUTPUT_TAG} f_DPS=${FULL_F_DPS} stat=${STAT_ABS} cr_up=${CR_UP} cr_down=${CR_DOWN}"
