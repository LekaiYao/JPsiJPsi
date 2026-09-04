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
GEN_BASE="/eos/user/c/chensh/JPsiJPsi/GEN_nofilter/DPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc/direct"
INVENTORY="${OUTPUT}/input/gen_dps_inventory.tsv"
ACCEPTANCE="${REPO}/Data_driven/inputs/acceptance_sps_full10_v1.txt"
if [ -e "${OUTPUT}" ]; then
  echo "Refusing to overwrite existing tag: ${OUTPUT}" >&2
  exit 68
fi
test -s "${ACCEPTANCE}"
mkdir -p "${OUTPUT}/input"
{
  printf 'path\tsize_bytes\tmtime_epoch\tsha256\n'
  for index in $(seq 1 10); do
    file="${GEN_BASE}/DPS_2016_JJ_${index}.root"
    test -s "${file}"
    printf '%s\t%s\t%s\t%s\n' "${file}" "$(stat -c %s "${file}")" \
      "$(stat -c %Y "${file}")" "$(sha256sum "${file}" | awk '{print $1}')"
  done
} > "${INVENTORY}"
cd "${REPO}"
root -l -b -q \
  "Data_driven/calculate_mass_extrapolation.cpp+(\"${SUMMARY}\")" \
  > "${OUTPUT}/run.log" 2>&1
grep -q '^status=complete$' "${SUMMARY}"
test "$(awk -F= '$1=="all_cross_event_pairs"{print $2}' "${SUMMARY}")" -gt 0
test "$(awk -F= '$1=="pass_mJJ_gt_7p5"{print $2}' "${SUMMARY}")" -gt 0
{
  echo "status=candidate_complete_pending_user_confirmation"
  echo "artifact=gen_only_dps_mass_extrapolation"
  echo "run_tag=${TAG}"
  echo "root_version=$(root-config --version)"
  echo "git_head=$(git rev-parse HEAD)"
  echo "model=none_GEN_only_counting"
  echo "selection_per_jpsi=10_le_pt_le_40_abs_y_le_2"
  echo "mixing=GEN_pair_id_slot1_x_slot2_different_source_events"
  echo "epsilon_definition=unweighted_pass_mJJ_ge_7p5_over_all_cross_event_pairs"
  echo "acceptance_role=closure_diagnostic_only_not_epsilon_numerator_or_denominator"
  echo "epsilon_m_uncertainty_status=pending_not_included_not_zero"
  echo "gen_inventory_sha256=$(sha256sum "${INVENTORY}" | awk '{print $1}')"
  echo "acceptance_sha256=$(sha256sum "${ACCEPTANCE}" | awk '{print $1}')"
  sha256sum Data_driven/calculate_mass_extrapolation.cpp \
    Data_driven/run_mass_extrapolation.sh "${INVENTORY}" "${ACCEPTANCE}"
} > "${OUTPUT}/run.metadata.txt"
find "${OUTPUT}" -maxdepth 3 -type f ! -name artifact_checksums.txt -print0 \
  | sort -z | xargs -0 sha256sum > "${OUTPUT}/artifact_checksums.txt"
printf 'candidate_complete_pending_user_confirmation\n' > "${OUTPUT}/complete.marker"
echo "DATA_DRIVEN_MASS_EXTRAPOLATION_COMPLETE tag=${TAG}"
