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
DPS_BASE="${REPO}/Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/DPS_ntuple"
WEIGHT="${OUTPUT}/WeightDPS_current.root"
TEMPLATE="${OUTPUT}/nominal_dps_template_input.root"
INVENTORY="${OUTPUT}/input/dps_ntuple_inventory.tsv"
if [ -e "${OUTPUT}" ]; then
  echo "Refusing to overwrite existing tag: ${OUTPUT}" >&2
  exit 68
fi
mkdir -p "${OUTPUT}/logs" "${OUTPUT}/input"
{
  printf 'path\tsize_bytes\tmtime_epoch\tsha256\n'
  while IFS= read -r -d '' file; do
    hash=$(sha256sum "${file}" | awk '{print $1}')
    printf '%s\t%s\t%s\t%s\n' "${file}" "$(stat -c %s "${file}")" \
      "$(stat -c %Y "${file}")" "${hash}"
  done < <(find "${DPS_BASE}" -maxdepth 1 -type f \
    -name 'Ntuple_2016_DPS_*.root' -print0 | sort -z -V)
} > "${INVENTORY}"
INPUT_FILES=$(awk 'NR>1{n++} END{print n+0}' "${INVENTORY}")
test "${INPUT_FILES}" -eq 65
INPUT_BYTES=$(awk -F '\t' 'NR>1{s+=$2} END{printf "%.0f",s}' "${INVENTORY}")
INVENTORY_SHA=$(sha256sum "${INVENTORY}" | awk '{print $1}')
cd "${REPO}/Data_driven/inputs"
root -l -b -q "${REPO}/Data_driven/build_dps_reference.cpp+(\"${WEIGHT}\")" \
  > "${OUTPUT}/logs/build_weight.log" 2>&1
test -s "${WEIGHT}"
cd "${REPO}"
root -l -b -q \
  "Data_driven/prepare_dps_template_input.cpp+(\"${WEIGHT}\",\"${TEMPLATE}\")" \
  > "${OUTPUT}/logs/prepare_template.log" 2>&1
test -s "${TEMPLATE}"
SOURCE_ENTRIES=$(awk -F= '$1=="source_entries"{print $2}' "${OUTPUT}/logs/prepare_template.log")
ACCEPTED_ENTRIES=$(awk -F= '$1=="accepted_entries"{print $2}' "${OUTPUT}/logs/prepare_template.log")
REJECTED_ENTRIES=$(awk -F= '$1=="rejected_entries"{print $2}' "${OUTPUT}/logs/prepare_template.log")
test -n "${SOURCE_ENTRIES}" && test -n "${ACCEPTED_ENTRIES}" && test -n "${REJECTED_ENTRIES}"
{
  echo "status=complete"
  echo "tag=${TAG}"
  echo "dps_ntuple_files=${INPUT_FILES}"
  echo "source_entries=${SOURCE_ENTRIES}"
  echo "accepted_entries=${ACCEPTED_ENTRIES}"
  echo "rejected_entries=${REJECTED_ENTRIES}"
  echo "selection=mJJ_ge_7p5_pt_each_10_to_40_trigger_match_samePV"
  echo "root_version=$(root-config --version)"
} > "${OUTPUT}/summary.txt"

{
  echo "status=complete"
  echo "tag=${TAG}"
  echo "root_version=$(root-config --version)"
  echo "git_head=$(git rev-parse HEAD)"
  echo "dps_ntuple_base=${DPS_BASE}"
  echo "dps_ntuple_files=${INPUT_FILES}"
  echo "dps_ntuple_total_size_bytes=${INPUT_BYTES}"
  echo "dps_ntuple_inventory=${INVENTORY}"
  echo "dps_ntuple_inventory_sha256=${INVENTORY_SHA}"
  echo "selection=mJJ_ge_7p5_pt_each_10_to_40_trigger_match_samePV"
  echo "seed_policy=no_random_seed_deterministic_rephrase"
  echo "corrections=frozen_nominal_acceptance_and_efficiency"
  sha256sum "${WEIGHT}" "${TEMPLATE}" \
    "${REPO}/Data_driven/inputs/input_manifest.txt" \
    "${REPO}/Data_driven/inputs/acceptance_sps_full10_v1.txt" \
    "${REPO}/Data_driven/inputs/efficiency_sps0p8_dps0p2_dedup60_v1.txt" \
    "${REPO}/Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/rephrase.cpp" \
    Data_driven/build_dps_reference.cpp \
    Data_driven/prepare_dps_template_input.cpp \
    Data_driven/run_build_dps_reference.sh
} > "${OUTPUT}/run.metadata.txt"
sha256sum "${OUTPUT}/summary.txt" "${OUTPUT}/run.metadata.txt" \
  "${WEIGHT}" "${TEMPLATE}" > "${OUTPUT}/artifact_checksums.txt"
printf 'complete\n' > "${OUTPUT}/complete.marker"
echo "DPS_REFERENCE_COMPLETE tag=${TAG}"
