#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <central-pp-tag> <dps-reference-tag> <unique-output-tag>" >&2
  exit 64
fi

CENTRAL_TAG="$1"
DPS_REFERENCE_TAG="$2"
OUTPUT_TAG="$3"
for tag in "${CENTRAL_TAG}" "${DPS_REFERENCE_TAG}" "${OUTPUT_TAG}"; do
  [[ "${tag}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || {
    echo "Unsafe tag: ${tag}" >&2
    exit 64
  }
done

REPO="/eos/home-l/leyao/26JJ/JPsiJPsi"
source "${REPO}/Data_driven/require_root640.sh"

CENTRAL="${REPO}/Data_driven/results/${CENTRAL_TAG}"
DPS_REFERENCE="${REPO}/Data_driven/results/${DPS_REFERENCE_TAG}"
OUTPUT="${REPO}/Data_driven/results/${OUTPUT_TAG}"
FIT_INPUT="${CENTRAL}/pp_data_2d/fit_results.csv"
PP_SUMMARY="${CENTRAL}/pp_data_2d/summary.txt"
CENTRAL_TEMPLATE_SUMMARY="${CENTRAL}/templates_2d/summary.txt"
DPS_MC_INPUT="${DPS_REFERENCE}/nominal_dps_template_input.root"
TEMPLATE_DIR="${OUTPUT}/templates_2d"
INPUT_MANIFEST="${REPO}/Data_driven/inputs/input_manifest.txt"
ACCEPTANCE="${REPO}/Data_driven/inputs/acceptance_sps_full10_v1.txt"
EFFICIENCY="${REPO}/Data_driven/inputs/efficiency_sps0p8_dps0p2_dedup60_v1.txt"

for required in "${CENTRAL}/complete.marker" "${CENTRAL}/run.metadata.txt" \
  "${FIT_INPUT}" "${PP_SUMMARY}" "${CENTRAL_TEMPLATE_SUMMARY}" \
  "${DPS_REFERENCE}/complete.marker" "${DPS_REFERENCE}/run.metadata.txt" \
  "${DPS_REFERENCE}/artifact_checksums.txt" "${DPS_REFERENCE}/summary.txt" \
  "${DPS_MC_INPUT}" \
  "${INPUT_MANIFEST}" "${ACCEPTANCE}" "${EFFICIENCY}"; do
  test -s "${required}" || {
    echo "Missing required input: ${required}" >&2
    exit 66
  }
done

(
  cd "${DPS_REFERENCE}"
  sha256sum -c artifact_checksums.txt >/dev/null
)

if [ -e "${OUTPUT}" ]; then
  echo "Refusing to overwrite existing tag: ${OUTPUT}" >&2
  exit 68
fi

test "$(awk -F= '$1=="accepted_fits"{print $2}' "${PP_SUMMARY}")" -eq 12
test "$(awk -F= '$1=="failed_fits"{print $2}' "${PP_SUMMARY}")" -eq 0
test "$(awk -F= '$1=="cells"{print $2}' "${PP_SUMMARY}")" -eq 12

ACTUAL_ACCEPTANCE_SHA=$(sha256sum "${ACCEPTANCE}" | awk '{print $1}')
ACTUAL_EFFICIENCY_SHA=$(sha256sum "${EFFICIENCY}" | awk '{print $1}')
CENTRAL_ACCEPTANCE_SHA=$(awk -F= '$1=="acceptance_sha256"{print $2}' \
  "${CENTRAL}/run.metadata.txt")
CENTRAL_EFFICIENCY_SHA=$(awk -F= '$1=="efficiency_sha256"{print $2}' \
  "${CENTRAL}/run.metadata.txt")
test "${ACTUAL_ACCEPTANCE_SHA}" = "${CENTRAL_ACCEPTANCE_SHA}"
test "${ACTUAL_EFFICIENCY_SHA}" = "${CENTRAL_EFFICIENCY_SHA}"

mkdir -p "${TEMPLATE_DIR}" "${OUTPUT}/logs"
cd "${REPO}"

root -l -b -q \
  "Data_driven/build_2d_templates_adaptive.cpp+(\"${TEMPLATE_DIR}\",\"${FIT_INPUT}\",\"${DPS_MC_INPUT}\",1.8,1.5707963267948966,0.0,false,\"dps_reconstruction_mc\")" \
  > "${OUTPUT}/logs/build_templates.log" 2>&1

TEMPLATE_SUMMARY="${TEMPLATE_DIR}/summary.txt"
test -s "${TEMPLATE_SUMMARY}"
test "$(awk -F= '$1=="cells"{print $2}' "${TEMPLATE_SUMMARY}")" -eq 12
test "$(awk -F= '$1=="dps_template_source"{print $2}' "${TEMPLATE_SUMMARY}")" = \
  "dps_reconstruction_mc"
DPS_MC_ENTRIES=$(awk -F= '$1=="accepted_entries"{print $2}' \
  "${DPS_REFERENCE}/summary.txt")
test "$(awk -F= '$1=="mixed_pairs"{print $2}' "${TEMPLATE_SUMMARY}")" \
  -eq "${DPS_MC_ENTRIES}"
test "$(awk -F= '$1=="sps_leakage_fraction"{print $2}' "${TEMPLATE_SUMMARY}")" = "0"
test "$(awk -F= '$1=="allow_signed_weights"{print $2}' "${TEMPLATE_SUMMARY}")" = "0"

F_DPS=$(awk -F= '$1=="f_dps"{print $2}' "${TEMPLATE_SUMMARY}")
F_SPS=$(awk -F= '$1=="f_sps"{print $2}' "${TEMPLATE_SUMMARY}")
DATA_CONTROL=$(awk -F= '$1=="data_control"{print $2}' "${TEMPLATE_SUMMARY}")
DPS_MC_CONTROL=$(awk -F= '$1=="mix_control"{print $2}' "${TEMPLATE_SUMMARY}")
ALPHA=$(awk -F= '$1=="alpha"{print $2}' "${TEMPLATE_SUMMARY}")
DATA_TOTAL=$(awk -F= '$1=="data_total"{print $2}' "${TEMPLATE_SUMMARY}")
DPS_TOTAL=$(awk -F= '$1=="dps_total"{print $2}' "${TEMPLATE_SUMMARY}")
SPS_TOTAL=$(awk -F= '$1=="sps_total"{print $2}' "${TEMPLATE_SUMMARY}")
NEGATIVE_SPS_BINS=$(awk -F= '$1=="negative_sps_bins"{print $2}' \
  "${TEMPLATE_SUMMARY}")
BASELINE_F_DPS=$(awk -F= '$1=="f_dps"{print $2}' "${CENTRAL_TEMPLATE_SUMMARY}")
BASELINE_DATA_CONTROL=$(awk -F= '$1=="data_control"{print $2}' \
  "${CENTRAL_TEMPLATE_SUMMARY}")
test -n "${F_DPS}" && test -n "${F_SPS}"
awk -v f="${F_DPS}" -v s="${F_SPS}" \
  'BEGIN{d=f+s-1; if(d<0)d=-d; exit !(f>=0 && f<=1 && s>=0 && s<=1 && d<1e-9)}'
awk -v a="${DATA_CONTROL}" -v b="${BASELINE_DATA_CONTROL}" \
  'BEGIN{d=a-b; if(d<0)d=-d; exit !(d<1e-9)}'

DPS_MC_TOTAL=$(awk -v d="${DPS_TOTAL}" -v a="${ALPHA}" \
  'BEGIN{printf "%.12g",d/a}')
DPS_MC_CR_FRACTION=$(awk -v c="${DPS_MC_CONTROL}" -v t="${DPS_MC_TOTAL}" \
  'BEGIN{printf "%.12g",c/t}')
ABS_SHIFT=$(awk -v f="${F_DPS}" -v b="${BASELINE_F_DPS}" \
  'BEGIN{printf "%+.12g",f-b}')
REL_SHIFT=$(awk -v f="${F_DPS}" -v b="${BASELINE_F_DPS}" \
  'BEGIN{printf "%+.12g",f/b-1}')
MIN_COVQUAL=$(awk -F, 'NR==2{m=$15} NR>1 && $15<m{m=$15} END{print m}' \
  "${FIT_INPUT}")
MAX_EDM=$(awk -F, 'NR==2{m=$16} NR>1 && $16>m{m=$16} END{printf "%.12g",m}' \
  "${FIT_INPUT}")
COMB_COMB_FIXED=$(awk -F, 'NR>1 && $24==1{n++} END{print n+0}' "${FIT_INPUT}")

{
  echo "status=candidate_complete_pending_user_confirmation"
  echo "artifact=fdps_dps_mc_template_cr_normalized"
  echo "method=dps_reconstruction_mc_template_singleCR_12cell_mergeddy"
  echo "f_dps=${F_DPS}"
  echo "f_sps=${F_SPS}"
  echo "central_event_mixing_f_dps=${BASELINE_F_DPS}"
  echo "absolute_shift_vs_event_mixing=${ABS_SHIFT}"
  echo "relative_shift_vs_event_mixing=${REL_SHIFT}"
  echo "control_region=abs_delta_y_ge_1p8_and_abs_delta_phi_le_pi_over_2"
  echo "sps_leakage_fraction=0"
  echo "data_control=${DATA_CONTROL}"
  echo "dps_mc_control_weight=${DPS_MC_CONTROL}"
  echo "dps_mc_total_weight=${DPS_MC_TOTAL}"
  echo "dps_mc_control_fraction=${DPS_MC_CR_FRACTION}"
  echo "alpha_data_per_dps_mc_weight=${ALPHA}"
  echo "data_total=${DATA_TOTAL}"
  echo "dps_total=${DPS_TOTAL}"
  echo "sps_total=${SPS_TOTAL}"
  echo "negative_sps_bins=${NEGATIVE_SPS_BINS}"
  echo "pp_fit_source_tag=${CENTRAL_TAG}"
  echo "pp_fit_accepted=12"
  echo "pp_fit_failed=0"
  echo "pp_fit_min_covQual=${MIN_COVQUAL}"
  echo "pp_fit_max_edm=${MAX_EDM}"
  echo "pp_fit_comb_comb_fixed_zero=${COMB_COMB_FIXED}"
  echo "new_roofit_fits_run=0"
  echo "uncertainty_analysis=paused_not_evaluated"
  echo "correction_scope=frozen_current_Data_driven_inputs_not_next_nominal_full_SPS_acceptance"
} > "${OUTPUT}/summary.txt"

MODEL_SHA=$(awk -F= '$1=="model_sha256"{print $2}' "${CENTRAL}/run.metadata.txt")
MANIFEST_SHA=$(sha256sum "${INPUT_MANIFEST}" | awk '{print $1}')
FIT_INPUT_SHA=$(sha256sum "${FIT_INPUT}" | awk '{print $1}')
DPS_MC_SHA=$(sha256sum "${DPS_MC_INPUT}" | awk '{print $1}')
{
  echo "status=candidate_complete_pending_user_confirmation"
  echo "output_tag=${OUTPUT_TAG}"
  echo "central_pp_tag=${CENTRAL_TAG}"
  echo "dps_reference_tag=${DPS_REFERENCE_TAG}"
  echo "root_version=$(root-config --version)"
  echo "git_head=$(git rev-parse HEAD)"
  echo "model_sha256=${MODEL_SHA}"
  echo "input_manifest=${INPUT_MANIFEST}"
  echo "input_manifest_sha256=${MANIFEST_SHA}"
  echo "acceptance=${ACCEPTANCE}"
  echo "acceptance_sha256=${ACTUAL_ACCEPTANCE_SHA}"
  echo "efficiency=${EFFICIENCY}"
  echo "efficiency_sha256=${ACTUAL_EFFICIENCY_SHA}"
  echo "pp_fit_input=${FIT_INPUT}"
  echo "pp_fit_input_sha256=${FIT_INPUT_SHA}"
  echo "dps_mc_input=${DPS_MC_INPUT}"
  echo "dps_mc_input_sha256=${DPS_MC_SHA}"
  echo "template_source=DPS_reconstruction_MC"
  echo "normalization=zero_SPS_leakage_single_control_region"
  echo "control_region=abs_delta_y_ge_1p8_and_abs_delta_phi_le_pi_over_2"
  echo "seed_policy=deterministic_no_random_seed"
  echo "pp_fit_policy=reuse_accepted_central_12cell_fit_no_refit"
  echo "uncertainty_analysis=paused"
  echo "correction_scope=frozen_current_Data_driven_inputs_not_next_nominal_full_SPS_acceptance"
  echo "generated_at=$(date --iso-8601=seconds)"
  sha256sum "${CENTRAL}/run.metadata.txt" "${PP_SUMMARY}" \
    "${DPS_REFERENCE}/run.metadata.txt" \
    Data_driven/build_2d_templates_adaptive.cpp \
    Data_driven/run_fdps_dps_mc_template.sh
} > "${OUTPUT}/run.metadata.txt"

find "${OUTPUT}" -maxdepth 3 -type f ! -name artifact_checksums.txt -print0 \
  | sort -z | xargs -0 sha256sum > "${OUTPUT}/artifact_checksums.txt"
printf 'candidate_complete_pending_user_confirmation\n' > "${OUTPUT}/complete.marker"

echo "DATA_DRIVEN_DPS_MC_TEMPLATE_COMPLETE tag=${OUTPUT_TAG} f_DPS=${F_DPS}"
