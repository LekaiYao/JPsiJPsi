#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <completed-systematics-tag> <unique-postprocess-tag>" >&2
  exit 64
fi
SYSTEMATICS_TAG="$1"
POSTPROCESS_TAG="$2"
for tag in "${SYSTEMATICS_TAG}" "${POSTPROCESS_TAG}"; do
  [[ "${tag}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || exit 64
done

REPO="/eos/home-l/leyao/26JJ/JPsiJPsi"
OUTPUT="${REPO}/Data_driven/results/${POSTPROCESS_TAG}"
cd "${REPO}"
python3 Data_driven/postprocess/scripts/prepare_inputs.py \
  "${SYSTEMATICS_TAG}" "${POSTPROCESS_TAG}"
python3 Data_driven/postprocess/scripts/calculate_results.py \
  "${OUTPUT}/inputs.json" "${OUTPUT}"
python3 Data_driven/postprocess/scripts/plot_figures.py \
  "${OUTPUT}/inputs.json" "${OUTPUT}"
{
  echo "status=candidate_complete_pending_user_confirmation"
  echo "systematics_tag=${SYSTEMATICS_TAG}"
  echo "postprocess_tag=${POSTPROCESS_TAG}"
  echo "git_head=$(git rev-parse HEAD)"
  sha256sum "${OUTPUT}/inputs.json" "${OUTPUT}/results.json" \
    Data_driven/postprocess/scripts/prepare_inputs.py \
    Data_driven/postprocess/current_pair_cross_section.json \
    Data_driven/postprocess/scripts/calculate_results.py \
    Data_driven/postprocess/scripts/plot_figures.py
} > "${OUTPUT}/run.metadata.txt"
printf 'candidate_complete_pending_user_confirmation\n' > "${OUTPUT}/complete.marker"
echo "DATA_DRIVEN_POSTPROCESS_COMPLETE tag=${POSTPROCESS_TAG}"
