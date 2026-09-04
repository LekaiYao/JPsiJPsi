#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <central-tag> <dps-reference-tag> <unique-comparison-tag>" >&2
  exit 64
fi
CENTRAL_TAG="$1"
REFERENCE_TAG="$2"
COMPARISON_TAG="$3"
for tag in "${CENTRAL_TAG}" "${REFERENCE_TAG}" "${COMPARISON_TAG}"; do
  [[ "${tag}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || {
    echo "Unsafe tag: ${tag}" >&2
    exit 64
  }
done

REPO="/eos/home-l/leyao/26JJ/JPsiJPsi"
source "${REPO}/Data_driven/require_root640.sh"
CENTRAL="${REPO}/Data_driven/results/${CENTRAL_TAG}"
REFERENCE="${REPO}/Data_driven/results/${REFERENCE_TAG}"
OUTPUT="${REPO}/Data_driven/results/${COMPARISON_TAG}"
MIXING="${CENTRAL}/mixed_dps_allpairs.root"
DPS_MC="${REFERENCE}/nominal_dps_template_input.root"

for required in "${MIXING}" "${DPS_MC}" "${CENTRAL}/run.metadata.txt" \
  "${REFERENCE}/run.metadata.txt" "${REFERENCE}/complete.marker"; do
  test -s "${required}" || { echo "Missing required input: ${required}" >&2; exit 66; }
done
if [ -e "${OUTPUT}" ]; then
  echo "Refusing to overwrite existing comparison tag: ${OUTPUT}" >&2
  exit 68
fi
mkdir -p "${OUTPUT}/logs"
cd "${REPO}"

root -l -b -q \
  "Data_driven/compare_mixed_dps_to_mc_1d.cpp+(\"${OUTPUT}\",\"${MIXING}\",\"${DPS_MC}\")" \
  > "${OUTPUT}/logs/compare.log" 2>&1
while IFS= read -r -d '' pdf; do
  stem=${pdf%.pdf}
  pdftoppm -png -singlefile -r 120 "${pdf}" "${stem}"
done < <(find "${OUTPUT}/plots" -maxdepth 1 -type f -name '*.pdf' -print0 | sort -z)

grep -q '^DPS_SHAPE_COMPARISON_COMPLETE ' "${OUTPUT}/logs/compare.log"
grep -q '^status=complete_pending_user_confirmation$' "${OUTPUT}/summary.txt"
grep -q '^mix_invalid_weights=0$' "${OUTPUT}/summary.txt"
grep -q '^mc_invalid_weights=0$' "${OUTPUT}/summary.txt"
test "$(awk -F, 'NR>1{n++} END{print n+0}' "${OUTPUT}/shape_metrics.csv")" -eq 9
test "$(find "${OUTPUT}/plots" -maxdepth 1 -type f -name '*.pdf' | wc -l)" -eq 9
test "$(find "${OUTPUT}/plots" -maxdepth 1 -type f -name '*.png' | wc -l)" -eq 9
awk -F, 'NR>1 && ($5!=0 || $12!=0){bad++} END{exit (bad+0)!=0}' \
  "${OUTPUT}/shape_metrics.csv"

{
  echo "artifact=normalized_1d_event_mixing_vs_dps_mc_shape_validation"
  echo "status=complete_pending_user_confirmation"
  echo "comparison_tag=${COMPARISON_TAG}"
  echo "central_tag=${CENTRAL_TAG}"
  echo "dps_reference_tag=${REFERENCE_TAG}"
  echo "root_version=$(root-config --version)"
  echo "git_head=$(git rev-parse HEAD)"
  echo "png_rendering=pdftoppm_from_pdf_at_120_dpi"
  echo "selection=mJJ_ge_7p5_pt_each_10_to_40_trigger_match_samePV"
  echo "corrections=same_frozen_acceptance_and_efficiency_tables_for_mixing_and_dps_mc"
  echo "mixing=deterministic_all_cross_event_Jpsi1_x_Jpsi2_saved_prompt_sweights"
  echo "mixing_uncertainty=delete_one_original_data_source_event_cluster_jackknife"
  echo "dps_mc_uncertainty=independent_weighted_event_normalized_delta_covariance"
  echo "normalization=unit_area_per_variable_full_listed_range"
  echo "systematic_interpretation=false"
  sha256sum "${CENTRAL}/run.metadata.txt" "${MIXING}" \
    "${REFERENCE}/run.metadata.txt" "${DPS_MC}" \
    Data_driven/inputs/input_manifest.txt \
    Data_driven/inputs/acceptance_sps_full10_v1.txt \
    Data_driven/inputs/efficiency_sps0p8_dps0p2_dedup60_v1.txt \
    Data_driven/compare_mixed_dps_to_mc_1d.cpp \
    Data_driven/run_compare_mixed_dps_to_mc_1d.sh
} > "${OUTPUT}/run.metadata.txt"

find "${OUTPUT}" -maxdepth 2 -type f ! -name artifact_checksums.txt -print0 \
  | sort -z | xargs -0 sha256sum > "${OUTPUT}/artifact_checksums.txt"
printf 'complete_pending_user_confirmation\n' > "${OUTPUT}/complete.marker"
echo "DPS_SHAPE_COMPARISON_TAG_COMPLETE tag=${COMPARISON_TAG}"
