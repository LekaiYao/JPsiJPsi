#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <unique-reference-tag>" >&2
  exit 64
fi
TAG="$1"
[[ "${TAG}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || exit 64

REPO="/eos/home-l/leyao/26JJ/JPsiJPsi"
source "${REPO}/Data_driven/require_root640.sh"
OUTPUT="${REPO}/Data_driven/results/${TAG}"
WEIGHT="${OUTPUT}/WeightDPS_current.root"
TEMPLATE="${OUTPUT}/nominal_dps_template_input.root"
if [ -e "${OUTPUT}" ]; then
  echo "Refusing to overwrite existing tag: ${OUTPUT}" >&2
  exit 68
fi
mkdir -p "${OUTPUT}/logs"
cd "${REPO}/Data_driven/inputs"
root -l -b -q "${REPO}/Data_driven/build_dps_reference.cpp+(\"${WEIGHT}\")" \
  > "${OUTPUT}/logs/build_weight.log" 2>&1
test -s "${WEIGHT}"
cd "${REPO}"
root -l -b -q \
  "Data_driven/prepare_dps_template_input.cpp+(\"${WEIGHT}\",\"${TEMPLATE}\")" \
  > "${OUTPUT}/logs/prepare_template.log" 2>&1
test -s "${TEMPLATE}"
{
  echo "status=complete"
  echo "tag=${TAG}"
  echo "root_version=$(root-config --version)"
  echo "git_head=$(git rev-parse HEAD)"
  sha256sum "${WEIGHT}" "${TEMPLATE}" \
    Data_driven/build_dps_reference.cpp \
    Data_driven/prepare_dps_template_input.cpp
} > "${OUTPUT}/run.metadata.txt"
echo "DPS_REFERENCE_COMPLETE tag=${TAG}"
