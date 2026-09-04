#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <nominal-half-dpsmc-tag> <dps-reference-tag> <unique-output-tag>" >&2
  exit 64
fi

NOMINAL_TAG="$1"
REFERENCE_TAG="$2"
OUTPUT_TAG="$3"
for tag in "${NOMINAL_TAG}" "${REFERENCE_TAG}" "${OUTPUT_TAG}"; do
  [[ "${tag}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || {
    echo "Unsafe tag: ${tag}" >&2
    exit 64
  }
done

REPO="/eos/home-l/leyao/26JJ/JPsiJPsi"
source "${REPO}/Data_driven/require_root640.sh"
NOMINAL="${REPO}/Data_driven/results/${NOMINAL_TAG}"
REFERENCE="${REPO}/Data_driven/results/${REFERENCE_TAG}"
OUTPUT="${REPO}/Data_driven/results/${OUTPUT_TAG}"
SEED=$(awk -F= '$1=="stat_seed"{print $2}' "${NOMINAL}/summary.txt")
CENTRAL_TAG=$(awk -F= '$1=="central_pp_tag"{print $2}' \
  "${NOMINAL}/run.metadata.txt")
CENTRAL="${REPO}/Data_driven/results/${CENTRAL_TAG}"
HALF_DIR="${NOMINAL}/half_sample_seed${SEED}"
DIRECT_HALF="${HALF_DIR}/dps_mc_half.root"
SELECTED_ENTRIES="${HALF_DIR}/selected_entries.csv"
FIT_INPUT="${CENTRAL}/pp_data_2d/fit_results.csv"
FIT_SUMMARY="${CENTRAL}/pp_data_2d/summary.txt"
DPS_WEIGHT="${REFERENCE}/WeightDPS_current.root"
REFERENCE_INVENTORY="${REFERENCE}/input/dps_ntuple_inventory.tsv"
DPS_BASE="${REPO}/Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/DPS_ntuple"
CANDIDATES="${OUTPUT}/dps_mc_candidates.root"
MIXED="${OUTPUT}/dps_mc_half_crossslot_mixed.root"
MIX_SUMMARY="${OUTPUT}/mixing_summary.txt"
TEMPLATES="${OUTPUT}/templates_2d"
COMPARISON="${OUTPUT}/shape_comparison_1d"
CURRENT_INVENTORY="${OUTPUT}/input/dps_ntuple_inventory.tsv"
INPUT_MANIFEST="${REPO}/Data_driven/inputs/input_manifest.txt"
ACCEPTANCE="${REPO}/Data_driven/inputs/acceptance_sps_full10_v1.txt"
EFFICIENCY="${REPO}/Data_driven/inputs/efficiency_sps0p8_dps0p2_dedup60_v1.txt"

for required in \
  "${NOMINAL}/complete.marker" "${NOMINAL}/summary.txt" \
  "${NOMINAL}/run.metadata.txt" "${NOMINAL}/artifact_checksums.txt" \
  "${DIRECT_HALF}" "${SELECTED_ENTRIES}" \
  "${REFERENCE}/complete.marker" "${REFERENCE}/run.metadata.txt" \
  "${REFERENCE}/artifact_checksums.txt" "${DPS_WEIGHT}" \
  "${REFERENCE_INVENTORY}" "${FIT_INPUT}" "${FIT_SUMMARY}" \
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

test "$(cat "${NOMINAL}/complete.marker")" = \
  "nominal_complete_user_confirmed"
test "${SEED}" = "20260902"
test "$(awk -F= '$1=="stat_selected_entries"{print $2}' \
  "${NOMINAL}/summary.txt")" -eq 2584
test "$(awk -F= '$1=="stat_input_entries"{print $2}' \
  "${NOMINAL}/summary.txt")" -eq 5168
(
  cd "${NOMINAL}"
  sha256sum -c artifact_checksums.txt >/dev/null
)
(
  cd "${REFERENCE}"
  sha256sum -c artifact_checksums.txt >/dev/null
)

mkdir -p "${OUTPUT}/logs" "${OUTPUT}/input" "${TEMPLATES}"
{
  printf 'path\tsize_bytes\tmtime_epoch\tsha256\n'
  while IFS= read -r -d '' file; do
    hash=$(sha256sum "${file}" | awk '{print $1}')
    printf '%s\t%s\t%s\t%s\n' "${file}" "$(stat -c %s "${file}")" \
      "$(stat -c %Y "${file}")" "${hash}"
  done < <(find "${DPS_BASE}" -maxdepth 1 -type f \
    -name 'Ntuple_2016_DPS_*.root' -print0 | sort -z -V)
} > "${CURRENT_INVENTORY}"
cmp -s "${CURRENT_INVENTORY}" "${REFERENCE_INVENTORY}"

cd "${REPO}"
root -l -b -q \
  "Data_driven/build_jpsi12_candidates.cpp+(\"${CANDIDATES}\",\"dps_mc\",\"\")" \
  > "${OUTPUT}/logs/build_candidates.log" 2>&1
FULL_CANDIDATE_EVENTS=$(awk -F= '$1=="selected_events"{print $2}' \
  "${OUTPUT}/logs/build_candidates.log")
test "${FULL_CANDIDATE_EVENTS}" -eq 5168

root -l -b -q \
  "Data_driven/mix_dps_mc_half_crossslot.cpp+(\"${CANDIDATES}\",\"${SELECTED_ENTRIES}\",\"${DIRECT_HALF}\",\"${MIXED}\",\"${MIX_SUMMARY}\")" \
  > "${OUTPUT}/logs/build_mixing.log" 2>&1
grep -q '^DPS_MC_HALF_CROSSSLOT_MIXING_COMPLETE$' \
  "${OUTPUT}/logs/build_mixing.log"
grep -q '^status=complete$' "${MIX_SUMMARY}"
test "$(awk -F= '$1=="full_candidate_events"{print $2}' \
  "${MIX_SUMMARY}")" -eq 5168
test "$(awk -F= '$1=="selected_pair_events"{print $2}' \
  "${MIX_SUMMARY}")" -eq 2584
test "$(awk -F= '$1=="direct_half_validation_mismatches"{print $2}' \
  "${MIX_SUMMARY}")" -eq 0

root -l -b -q \
  "Data_driven/build_2d_templates_adaptive.cpp+(\"${TEMPLATES}\",\"${FIT_INPUT}\",\"${MIXED}\",1.8,1.5707963267948966,0.0,false,\"dps_mc_fixed_half_crossslot_event_mixing\")" \
  > "${OUTPUT}/logs/build_templates.log" 2>&1
test "$(awk -F= '$1=="dps_template_source"{print $2}' \
  "${TEMPLATES}/summary.txt")" = \
  "dps_mc_fixed_half_crossslot_event_mixing"

root -l -b -q \
  "Data_driven/compare_mixed_dps_to_mc_1d.cpp+(\"${COMPARISON}\",\"${MIXED}\",\"${DIRECT_HALF}\",\"dps_mc\",\"DPS MC event mixing\",\"DPS MC direct\",\"DPS MC mixing closure\",\"shared_dps_mc_source_events_cross_covariance_not_propagated\")" \
  > "${OUTPUT}/logs/compare_1d.log" 2>&1
grep -q '^DPS_SHAPE_COMPARISON_COMPLETE ' "${OUTPUT}/logs/compare_1d.log"
grep -q '^mix_source_events=2584$' "${COMPARISON}/summary.txt"
grep -q '^mix_invalid_weights=0$' "${COMPARISON}/summary.txt"
grep -q '^mc_invalid_weights=0$' "${COMPARISON}/summary.txt"
while IFS= read -r -d '' pdf; do
  stem=${pdf%.pdf}
  pdftoppm -png -singlefile -r 120 "${pdf}" "${stem}"
done < <(find "${COMPARISON}/plots" -maxdepth 1 -type f \
  -name '*.pdf' -print0 | sort -z)
test "$(find "${COMPARISON}/plots" -maxdepth 1 -type f \
  -name '*.pdf' | wc -l)" -eq 9
test "$(find "${COMPARISON}/plots" -maxdepth 1 -type f \
  -name '*.png' | wc -l)" -eq 9
test "$(awk -F, 'NR>1{n++} END{print n+0}' \
  "${COMPARISON}/shape_metrics.csv")" -eq 9

NOMINAL_F_DPS=$(awk -F= '$1=="f_dps"{print $2}' "${NOMINAL}/summary.txt")
MIXED_F_DPS=$(awk -F= '$1=="f_dps"{print $2}' "${TEMPLATES}/summary.txt")
SIGNED_SHIFT=$(awk -v mixed="${MIXED_F_DPS}" -v nominal="${NOMINAL_F_DPS}" \
  'BEGIN{printf "%+.12g",mixed-nominal}')
ABS_SHIFT=$(awk -v mixed="${MIXED_F_DPS}" -v nominal="${NOMINAL_F_DPS}" \
  'BEGIN{d=mixed-nominal; if(d<0)d=-d; printf "%.12g",d}')
RELATIVE_SHIFT_PERCENT=$(awk -v shift="${SIGNED_SHIFT}" \
  -v nominal="${NOMINAL_F_DPS}" \
  'BEGIN{printf "%+.12g",100*shift/nominal}')
ACCEPTED_PAIRS=$(awk -F= '$1=="accepted_pairs"{print $2}' "${MIX_SUMMARY}")
MASS_REJECTED=$(awk -F= '$1=="mass_rejected"{print $2}' "${MIX_SUMMARY}")
CORRECTION_REJECTED=$(awk -F= '$1=="correction_rejected"{print $2}' \
  "${MIX_SUMMARY}")
NEGATIVE_SPS_BINS=$(awk -F= '$1=="negative_sps_bins"{print $2}' \
  "${TEMPLATES}/summary.txt")
MIN_COVQUAL=$(awk -F, 'NR==2{m=$15} NR>1&&$15<m{m=$15} END{print m}' \
  "${FIT_INPUT}")
MAX_EDM=$(awk -F, 'NR==2{m=$16} NR>1&&$16>m{m=$16} END{print m}' \
  "${FIT_INPUT}")
FIXED_ZERO=$(awk -F= '$1=="comb_comb_fixed_zero_fits"{print $2}' \
  "${FIT_SUMMARY}")

{
  echo "status=complete_pending_user_confirmation"
  echo "artifact=dps_mc_half_event_mixing_systematic_and_1d_shape_check"
  echo "nominal_tag=${NOMINAL_TAG}"
  echo "nominal_pair_events=2584"
  echo "nominal_single_jpsi_candidates=5168"
  echo "full_dps_mc_pair_events=5168"
  echo "mixing=deterministic_all_cross_Jpsi1_A_x_Jpsi2_B_A_not_equal_B"
  echo "mixing_seed=none"
  echo "mixed_accepted_pairs=${ACCEPTED_PAIRS}"
  echo "mixed_mass_rejected=${MASS_REJECTED}"
  echo "mixed_correction_rejected=${CORRECTION_REJECTED}"
  echo "nominal_f_dps=${NOMINAL_F_DPS}"
  echo "mixed_dps_mc_f_dps=${MIXED_F_DPS}"
  echo "signed_shift_mixed_minus_nominal=${SIGNED_SHIFT}"
  echo "relative_shift_percent=${RELATIVE_SHIFT_PERCENT}"
  echo "dps_mc_usage_systematic=${ABS_SHIFT}"
  echo "systematic_definition=absolute_f_dps_shift_DPS_MC_crossslot_event_mixing_vs_direct_nominal"
  echo "negative_sps_bins=${NEGATIVE_SPS_BINS}"
  echo "shape_comparison_plots=9"
  echo "shape_metrics=max_CDF_Jensen_Shannon_and_diagnostic_covariance_chi2"
  echo "shape_shared_source_cross_covariance=not_propagated"
  echo "pp_fits_reused=12"
  echo "pp_fit_failed=0"
  echo "pp_fit_min_covQual=${MIN_COVQUAL}"
  echo "pp_fit_max_edm=${MAX_EDM}"
  echo "pp_fit_comb_comb_fixed_zero=${FIXED_ZERO}"
  echo "new_roofit_fits_run=0"
  echo "correction_scope=frozen_current_Data_driven_inputs_not_next_pairconditioned_acceptance"
} > "${OUTPUT}/summary.txt"

MODEL_SHA=$(awk -F= '$1=="model_sha256"{print $2}' \
  "${NOMINAL}/run.metadata.txt")
MANIFEST_SHA=$(sha256sum "${INPUT_MANIFEST}" | awk '{print $1}')
ACCEPTANCE_SHA=$(sha256sum "${ACCEPTANCE}" | awk '{print $1}')
EFFICIENCY_SHA=$(sha256sum "${EFFICIENCY}" | awk '{print $1}')
INVENTORY_SHA=$(sha256sum "${CURRENT_INVENTORY}" | awk '{print $1}')
{
  echo "status=complete_pending_user_confirmation"
  echo "output_tag=${OUTPUT_TAG}"
  echo "nominal_tag=${NOMINAL_TAG}"
  echo "dps_reference_tag=${REFERENCE_TAG}"
  echo "central_pp_tag=${CENTRAL_TAG}"
  echo "root_version=$(root-config --version)"
  echo "git_head=$(git rev-parse HEAD)"
  echo "model_sha256=${MODEL_SHA}"
  echo "input_manifest_sha256=${MANIFEST_SHA}"
  echo "acceptance_sha256=${ACCEPTANCE_SHA}"
  echo "efficiency_sha256=${EFFICIENCY_SHA}"
  echo "dps_ntuple_inventory_sha256=${INVENTORY_SHA}"
  echo "nominal_seed=${SEED}"
  echo "mixing_seed_policy=deterministic_all_cross_no_seed"
  echo "selection=mJJ_ge_7p5_pt_each_10_to_40_trigger_match_samePV"
  echo "pair_weight=recomputed_frozen_acceptance_efficiency_correction"
  echo "shape_normalization=unit_area_per_variable_full_listed_range"
  echo "shape_mix_uncertainty=delete_one_DPS_MC_source_event_cluster_jackknife"
  echo "shape_reference_uncertainty=independent_weighted_event_normalized_delta_covariance"
  echo "shape_shared_source_cross_covariance=not_propagated_metrics_diagnostic_only"
  echo "generated_at=$(date --iso-8601=seconds)"
  sha256sum "${NOMINAL}/run.metadata.txt" "${NOMINAL}/summary.txt" \
    "${DIRECT_HALF}" "${SELECTED_ENTRIES}" \
    "${REFERENCE}/run.metadata.txt" "${DPS_WEIGHT}" \
    "${REFERENCE_INVENTORY}" "${FIT_INPUT}" \
    "${INPUT_MANIFEST}" "${ACCEPTANCE}" "${EFFICIENCY}" \
    Data_driven/build_jpsi12_candidates.cpp \
    Data_driven/mix_dps_mc_half_crossslot.cpp \
    Data_driven/build_2d_templates_adaptive.cpp \
    Data_driven/compare_mixed_dps_to_mc_1d.cpp \
    Data_driven/run_dpsmc_half_eventmix_systematic.sh
} > "${OUTPUT}/run.metadata.txt"

find "${OUTPUT}" -maxdepth 4 -type f ! -name artifact_checksums.txt -print0 \
  | sort -z | xargs -0 sha256sum > "${OUTPUT}/artifact_checksums.txt"
printf 'complete_pending_user_confirmation\n' > "${OUTPUT}/complete.marker"
echo "DPS_MC_HALF_EVENTMIX_SYSTEMATIC_COMPLETE tag=${OUTPUT_TAG} nominal=${NOMINAL_F_DPS} mixed=${MIXED_F_DPS} systematic=${ABS_SHIFT}"
