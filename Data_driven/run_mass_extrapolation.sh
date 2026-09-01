#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <unique-mass-extrapolation-tag>" >&2
  exit 64
fi
TAG="$1"
[[ "${TAG}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || exit 64
REPO="/eos/home-l/leyao/26JJ/JPsiJPsi"
source "${REPO}/Data_driven/require_root640.sh"
OUTPUT="${REPO}/Data_driven/results/${TAG}"
SUMMARY="${OUTPUT}/mass_extrapolation.txt"
if [ -e "${OUTPUT}" ]; then
  echo "Refusing to overwrite existing tag: ${OUTPUT}" >&2
  exit 68
fi
mkdir -p "${OUTPUT}"
cd "${REPO}"
root -l -b -q \
  "Data_driven/calculate_mass_extrapolation.cpp+(\"${SUMMARY}\")" \
  > "${OUTPUT}/run.log" 2>&1
grep -q '^status=complete$' "${SUMMARY}"
{
  echo "status=candidate_complete_pending_user_confirmation"
  echo "run_tag=${TAG}"
  echo "root_version=$(root-config --version)"
  echo "git_head=$(git rev-parse HEAD)"
  sha256sum Data_driven/calculate_mass_extrapolation.cpp \
    Data_driven/inputs/acceptance_sps_full10_v1.txt
} > "${OUTPUT}/run.metadata.txt"
printf 'candidate_complete_pending_user_confirmation\n' > "${OUTPUT}/complete.marker"
echo "DATA_DRIVEN_MASS_EXTRAPOLATION_COMPLETE tag=${TAG}"
