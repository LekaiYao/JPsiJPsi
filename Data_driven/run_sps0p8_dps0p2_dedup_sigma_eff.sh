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
PAIR_JSON="${REPO}/Data_driven/postprocess/pair_cross_section_newdata167_diff35_v1_20260902.json"
PAIR_DIR="${REPO}/Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/fit_results/newdata167_accmix23_effmix19_lumi36p684_root640_diff35_v1_20260902"
PAIR_CSV="${PAIR_DIR}/cross_sections.csv"
PAIR_METADATA="${PAIR_DIR}/run.metadata.txt"
REFERENCE_JSON="${REPO}/Data_driven/postprocess/inputs.json"
CALCULATOR="${REPO}/Data_driven/postprocess/scripts/calculate_current_nominal_sigma_eff.py"

for required in \
  "${F_INPUT}/summary.txt" "${F_INPUT}/run.metadata.txt" \
  "${F_INPUT}/artifact_checksums.txt" "${F_INPUT}/complete.marker" \
  "${MASS}/mass_extrapolation.txt" "${MASS}/run.metadata.txt" \
  "${MASS}/artifact_checksums.txt" "${MASS}/complete.marker" \
  "${PAIR_JSON}" "${PAIR_CSV}" "${PAIR_METADATA}" \
  "${REFERENCE_JSON}" "${CALCULATOR}"; do
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
test "$(sha256sum "${PAIR_CSV}" | awk '{print $1}')" = \
  "$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["source_cross_sections_sha256"])' "${PAIR_JSON}")"
test "$(sha256sum "${PAIR_METADATA}" | awk '{print $1}')" = \
  "$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["source_metadata_sha256"])' "${PAIR_JSON}")"

mkdir -p "${OUTPUT}/logs"
python3 "${CALCULATOR}" \
  "${F_INPUT}/summary.txt" "${MASS}/mass_extrapolation.txt" \
  "${PAIR_JSON}" "${REFERENCE_JSON}" "${OUTPUT}" \
  > "${OUTPUT}/logs/calculate.log" 2>&1
RESULT_STATUS="candidate_complete_pending_user_confirmation"
grep -q "^status=${RESULT_STATUS}$" "${OUTPUT}/summary.txt"
grep -q '^pair_total_pending_systematic_sources=none$' "${OUTPUT}/summary.txt"
grep -q '^epsilon_m_uncertainty_status=pending_not_included_not_zero$' \
  "${OUTPUT}/summary.txt"
grep -q '^cov_f_dps_pair_total_stat=0_user_assumption$' \
  "${OUTPUT}/summary.txt"
grep -q '^cov_f_dps_pair_total_syst=0_user_assumption$' \
  "${OUTPUT}/summary.txt"

{
  echo "status=${RESULT_STATUS}"
  echo "artifact=sps0p8_dps0p2_dedup60_sigma_DPS_sigma_eff"
  echo "output_tag=${OUTPUT_TAG}"
  echo "f_dps_tag=${F_TAG}"
  echo "f_dps_marker=${F_MARKER}"
  echo "mass_extrapolation_tag=${MASS_TAG}"
  echo "pair_cross_section_scope=matching_sps0p8_dps0p2_generation_accepted_nominal"
  echo "pair_cross_section_status=accepted_nominal_user_confirmed"
  echo "cov_f_dps_pair_total_stat=0_user_assumption"
  echo "cov_f_dps_pair_total_syst=0_user_assumption"
  echo "epsilon_m_uncertainty_status=pending_not_included_not_zero"
  echo "result_status=candidate_pending_user_confirmation"
  echo "root_environment_version=$(root-config --version)"
  echo "python_version=$(python3 --version | awk '{print $2}')"
  echo "git_head=$(git -C "${REPO}" rev-parse HEAD)"
  echo "model_sha256=$(awk -F= '$1=="model_sha256"{print $2}' "${F_INPUT}/run.metadata.txt")"
  echo "input_manifest_sha256=$(sha256sum "${REPO}/Data_driven/inputs/input_manifest.txt" | awk '{print $1}')"
  echo "new_roofit_fits_run=0"
  sha256sum "${F_INPUT}/summary.txt" "${F_INPUT}/run.metadata.txt" \
    "${MASS}/mass_extrapolation.txt" "${MASS}/run.metadata.txt" \
    "${PAIR_JSON}" "${PAIR_CSV}" "${PAIR_METADATA}" "${REFERENCE_JSON}" \
    "${CALCULATOR}" "${REPO}/Data_driven/run_sps0p8_dps0p2_dedup_sigma_eff.sh"
} > "${OUTPUT}/run.metadata.txt"
find "${OUTPUT}" -maxdepth 3 -type f ! -name artifact_checksums.txt -print0 \
  | sort -z | xargs -0 sha256sum > "${OUTPUT}/artifact_checksums.txt"
printf '%s\n' "${RESULT_STATUS}" > "${OUTPUT}/complete.marker"
echo "SPS0P8_DPS0P2_DEDUP_SIGMA_EFF_COMPLETE tag=${OUTPUT_TAG}"
