#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <f-dps-tag> <mass-extrapolation-tag> <unique-output-tag>" >&2
  exit 64
fi
F_TAG="$1"
MASS_TAG="$2"
OUTPUT_TAG="$3"
for tag in "${F_TAG}" "${MASS_TAG}" "${OUTPUT_TAG}"; do
  [[ "${tag}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || exit 64
done

REPO="/eos/home-l/leyao/26JJ/JPsiJPsi"
source "${REPO}/Data_driven/require_root640.sh"
F_INPUT="${REPO}/Data_driven/results/${F_TAG}"
MASS="${REPO}/Data_driven/results/${MASS_TAG}"
OUTPUT="${REPO}/Data_driven/results/${OUTPUT_TAG}"
DATA_TAG="newdata167_accmix23_effmix19_sps0p85_dps0p15_lumi36p684_pairsym_cseed50_root640_strategy2_nominal_v4_20260902"
DATA_DIR="${REPO}/Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/fit_results/${DATA_TAG}"
TOTAL_SUMMARY="${DATA_DIR}/total/total_fit_summary.csv"
TOTAL_METADATA="${DATA_DIR}/total/run.metadata.txt"
SYS3_CSV="${DATA_DIR}/inputs/provenance/sys3.csv"
REFERENCE_JSON="${REPO}/Data_driven/postprocess/inputs.json"
PAIR_BUILDER="${REPO}/Data_driven/postprocess/scripts/build_total_only_pair_cross_section.py"
CALCULATOR="${REPO}/Data_driven/postprocess/scripts/calculate_current_nominal_sigma_eff.py"

for required in \
  "${F_INPUT}/summary.txt" "${F_INPUT}/run.metadata.txt" \
  "${F_INPUT}/artifact_checksums.txt" "${F_INPUT}/complete.marker" \
  "${MASS}/mass_extrapolation.txt" "${MASS}/run.metadata.txt" \
  "${MASS}/artifact_checksums.txt" "${MASS}/complete.marker" \
  "${TOTAL_SUMMARY}" "${TOTAL_METADATA}" "${SYS3_CSV}" \
  "${REFERENCE_JSON}" "${PAIR_BUILDER}" "${CALCULATOR}"; do
  test -s "${required}"
done
if [ -e "${OUTPUT}" ]; then
  echo "Refusing to overwrite existing tag: ${OUTPUT}" >&2
  exit 68
fi
F_MARKER="$(cat "${F_INPUT}/complete.marker")"
case "${F_MARKER}" in
  candidate_complete_pending_user_confirmation|nominal_complete_user_confirmed) ;;
  *) echo "Unsupported f_DPS marker: ${F_MARKER}" >&2; exit 65 ;;
esac
test "$(cat "${MASS}/complete.marker")" = \
  "candidate_complete_pending_user_confirmation"
(
  cd "${F_INPUT}"
  sha256sum -c artifact_checksums.txt >/dev/null
)
(
  cd "${MASS}"
  sha256sum -c artifact_checksums.txt >/dev/null
)

mkdir -p "${OUTPUT}/logs"
PAIR_JSON="${OUTPUT}/pair_cross_section_total_only.json"
python3 "${PAIR_BUILDER}" \
  "${TOTAL_SUMMARY}" "${TOTAL_METADATA}" "${SYS3_CSV}" "${PAIR_JSON}" \
  > "${OUTPUT}/logs/build_pair_cross_section.log" 2>&1
python3 "${CALCULATOR}" \
  "${F_INPUT}/summary.txt" "${MASS}/mass_extrapolation.txt" \
  "${PAIR_JSON}" "${REFERENCE_JSON}" "${OUTPUT}" \
  > "${OUTPUT}/logs/calculate.log" 2>&1

RESULT_STATUS="candidate_complete_included_sources_pending_pair_systematic_and_user_confirmation"
grep -q "^status=${RESULT_STATUS}$" "${OUTPUT}/summary.txt"
grep -q '^pair_total_pending_systematic_sources=fitter_stability$' \
  "${OUTPUT}/summary.txt"
grep -q '^epsilon_m_uncertainty_status=pending_not_included_not_zero$' \
  "${OUTPUT}/summary.txt"
grep -q '^cov_f_dps_pair_total_stat=0_user_assumption$' \
  "${OUTPUT}/summary.txt"
grep -q '^cov_f_dps_pair_total_syst=0_user_assumption$' \
  "${OUTPUT}/summary.txt"

{
  echo "status=${RESULT_STATUS}"
  echo "artifact=sps0p85_dps0p15_sigma_DPS_sigma_eff_total_only_pair_candidate"
  echo "output_tag=${OUTPUT_TAG}"
  echo "f_dps_tag=${F_TAG}"
  echo "f_dps_marker=${F_MARKER}"
  echo "mass_extrapolation_tag=${MASS_TAG}"
  echo "pair_cross_section_scope=matching_generation_total_only_temporary"
  echo "pair_cross_section_status=temporary_pending_Data_cross_section_and_user_confirmation"
  echo "pair_fitter_stability=pending_not_included_not_zero"
  echo "cov_f_dps_pair_total_stat=0_user_assumption"
  echo "cov_f_dps_pair_total_syst=0_user_assumption"
  echo "epsilon_m_uncertainty_status=pending_not_included_not_zero"
  echo "result_status=candidate_pending_user_confirmation_and_Data_cross_section"
  echo "root_environment_version=$(root-config --version)"
  echo "python_version=$(python3 --version | awk '{print $2}')"
  echo "git_head=$(git -C "${REPO}" rev-parse HEAD)"
  echo "model_sha256=$(awk -F= '$1=="model_sha256"{print $2}' "${F_INPUT}/run.metadata.txt")"
  echo "input_manifest_sha256=$(sha256sum "${REPO}/Data_driven/inputs/input_manifest.txt" | awk '{print $1}')"
  echo "new_roofit_fits_run=0"
  sha256sum "${F_INPUT}/summary.txt" "${F_INPUT}/run.metadata.txt" \
    "${MASS}/mass_extrapolation.txt" "${MASS}/run.metadata.txt" \
    "${TOTAL_SUMMARY}" "${TOTAL_METADATA}" "${SYS3_CSV}" \
    "${PAIR_JSON}" "${REFERENCE_JSON}" "${PAIR_BUILDER}" "${CALCULATOR}" \
    "${REPO}/Data_driven/run_sps0p85_dps0p15_sigma_eff.sh"
} > "${OUTPUT}/run.metadata.txt"
find "${OUTPUT}" -maxdepth 3 -type f ! -name artifact_checksums.txt -print0 \
  | sort -z | xargs -0 sha256sum > "${OUTPUT}/artifact_checksums.txt"
printf '%s\n' "${RESULT_STATUS}" > "${OUTPUT}/complete.marker"
echo "SPS0P85_DPS0P15_SIGMA_EFF_COMPLETE tag=${OUTPUT_TAG}"
