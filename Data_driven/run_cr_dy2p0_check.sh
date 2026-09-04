#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <unique-output-tag>" >&2
  exit 64
fi
OUTPUT_TAG="$1"
[[ "${OUTPUT_TAG}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || exit 64

REPO="/eos/home-l/leyao/26JJ/JPsiJPsi"
source "${REPO}/Data_driven/require_root640.sh"
NOMINAL_TAG="newdata167_accmix23_effmix19_fixedhalf_seed20260902_combcomb_sigcombfix0_root640_v2_20260902"
NOMINAL="${REPO}/Data_driven/results/${NOMINAL_TAG}"
OUTPUT="${REPO}/Data_driven/results/${OUTPUT_TAG}"
MODEL="${NOMINAL}/input/Model_4D_tot.root"
DATA="${NOMINAL}/input/WeightData_newdata167_accmix23_effmix19.root"
HALF="${NOMINAL}/nominal_fixed_half_seed20260902/dps_mc_half.root"
NOMINAL_SUMMARY="${NOMINAL}/summary.txt"
PP="${OUTPUT}/pp_data_dy2p0"
TEMPLATE="${OUTPUT}/template_fixed_half"
PHI90="1.5707963267948966"

for required in "${NOMINAL}/complete.marker" "${NOMINAL}/artifact_checksums.txt" \
  "${NOMINAL_SUMMARY}" "${NOMINAL}/run.metadata.txt" "${MODEL}" "${DATA}" \
  "${HALF}" "${REPO}/Data_driven/fit_pp_2d_adaptive.cpp" \
  "${REPO}/Data_driven/build_2d_templates_adaptive.cpp"; do
  test -s "${required}"
done
test "$(cat "${NOMINAL}/complete.marker")" = "nominal_complete_user_confirmed"
test "$(awk -F= '$1=="status"{print $2}' "${NOMINAL_SUMMARY}")" = \
  "nominal_complete_user_confirmed"
test "$(sha256sum "${MODEL}" | awk '{print $1}')" = \
  "9987bf734f9308f4177b290828cff82c7cb470d035cca8bbdba337534af61238"
test "$(sha256sum "${DATA}" | awk '{print $1}')" = \
  "8cd1abae2a896f5ebf82a805e62009e96a7f7e4102350c1f6e18c2568c1e8957"
(
  cd "${NOMINAL}"
  sha256sum -c artifact_checksums.txt >/dev/null
)
if [ -e "${OUTPUT}" ]; then
  echo "Refusing to overwrite existing tag: ${OUTPUT}" >&2
  exit 68
fi

mkdir -p "${OUTPUT}/logs"
cd "${REPO}"
root -l -b -q \
  "Data_driven/fit_pp_2d_adaptive.cpp+(\"${PP}\",\"${DATA}\",true,false,true,2.0,${PHI90},\"${MODEL}\")" \
  > "${OUTPUT}/logs/fit_pp_dy2p0.log" 2>&1
root -l -b -q \
  "Data_driven/build_2d_templates_adaptive.cpp+(\"${TEMPLATE}\",\"${PP}/fit_results.csv\",\"${HALF}\",2.0,${PHI90},0.0,false,\"dps_reconstruction_mc_fixed_seed_half\")" \
  > "${OUTPUT}/logs/build_template_dy2p0.log" 2>&1

python3 - "${NOMINAL_SUMMARY}" "${PP}/summary.txt" \
  "${PP}/fit_results.csv" "${TEMPLATE}/summary.txt" \
  "${TEMPLATE}/cells.csv" "${OUTPUT}/summary.txt" <<'PY'
import csv
import math
import pathlib
import sys


def read_kv(path):
    values = {}
    for line in pathlib.Path(path).read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key] = value
    return values


nominal_path, pp_path, fit_path, template_path, cells_path, output_path = sys.argv[1:]
nominal = read_kv(nominal_path)
pp = read_kv(pp_path)
template = read_kv(template_path)
fits = list(csv.DictReader(open(fit_path, encoding="utf-8")))
cells = list(csv.DictReader(open(cells_path, encoding="utf-8")))
if len(fits) != 15 or int(pp["cells"]) != 15:
    raise RuntimeError("dy2p0 geometry is not 15 cells")
if any(int(row["accepted"]) != 1 or int(row["status"]) != 0
       or int(row["covQual"]) < 2 or float(row["edm"]) >= 0.01
       for row in fits):
    raise RuntimeError("one or more dy2p0 PP fits failed the QA gate")
control_fits = [
    row for row in fits
    if math.isclose(float(row["dy_low"]), 2.0)
    and math.isclose(float(row["dy_high"]), 4.0)
    and math.isclose(float(row["dphi_low"]), 0.0)
    and math.isclose(float(row["dphi_high"]), math.pi / 2.0)
]
control_cells = [row for row in cells if int(row["in_control"]) == 1]
if len(control_fits) != 1 or len(control_cells) != 1:
    raise RuntimeError("dy2p0 control cell is not unique")
nominal_f = float(nominal["f_dps"])
variation_f = float(template["f_dps"])
shift = variation_f - nominal_f
relative = 100.0 * shift / nominal_f
control = control_fits[0]
lines = [
    "status=candidate_complete_pending_user_confirmation",
    "artifact=newdata_fixed_half_dpsmc_dy2p0_cr_check",
    "definition=abs_delta_y_ge_2p0_and_abs_delta_phi_lt_pi_over_2",
    "boundary_note=continuous_observables_existing_fit_convention_uses_delta_y_ge_lower_edge",
    "sps_leakage_fraction=0",
    f"nominal_f_dps={nominal_f:.12g}",
    f"dy2p0_f_dps={variation_f:.12g}",
    f"signed_shift={shift:.12g}",
    f"absolute_shift={abs(shift):.12g}",
    f"relative_shift_percent={relative:.12g}",
    f"absolute_relative_shift_percent={abs(relative):.12g}",
    f"data_control={float(template['data_control']):.12g}",
    f"data_control_error={float(template['data_control_error']):.12g}",
    f"dps_mc_control={float(template['mix_control']):.12g}",
    f"alpha={float(template['alpha']):.12g}",
    f"alpha_pair_level_error={float(template['alpha_pair_level_error']):.12g}",
    f"data_total={float(template['data_total']):.12g}",
    f"dps_total={float(template['dps_total']):.12g}",
    f"negative_sps_bins={int(template['negative_sps_bins'])}",
    f"control_tree_entries={int(control['tree_entries'])}",
    f"control_fit_entries={int(control['fit_entries'])}",
    f"control_effective_entries={float(control['effective_entries']):.12g}",
    f"control_pp_yield={float(control['pp_yield']):.12g}",
    f"control_pp_error={float(control['pp_error']):.12g}",
    f"pp_fits={len(fits)}",
    f"pp_accepted_fits={sum(int(row['accepted']) for row in fits)}",
    f"pp_min_covQual={min(int(row['covQual']) for row in fits)}",
    f"pp_max_edm={max(float(row['edm']) for row in fits):.12g}",
    f"pp_comb_comb_fixed_zero={sum(int(row['comb_comb_fixed_zero']) for row in fits)}",
    f"pp_sig_comb_fixed_zero={sum(int(row['sig_comb_fixed_zero']) for row in fits)}",
    "pp_fit_error_convention=ROOT_native_AsymptoticError_true",
    "pp_comb_comb_boundary_strategy=fix_zero_and_refit_on_root_range_exception",
    "pp_sig_comb_boundary_strategy=fix_shared_zero_for_Sig_Comb_and_Comb_Sig_and_refit_on_root_range_exception",
    "nominal_replacement_status=not_replaced_pending_user_confirmation",
    "new_roofit_fits_run=15",
]
pathlib.Path(output_path).write_text("\n".join(lines) + "\n", encoding="utf-8")
PY

grep -q '^status=candidate_complete_pending_user_confirmation$' \
  "${OUTPUT}/summary.txt"
grep -q '^pp_accepted_fits=15$' "${OUTPUT}/summary.txt"

{
  echo "status=candidate_complete_pending_user_confirmation"
  echo "artifact=newdata_fixed_half_dpsmc_dy2p0_cr_check"
  echo "output_tag=${OUTPUT_TAG}"
  echo "nominal_tag=${NOMINAL_TAG}"
  echo "root_version=$(root-config --version)"
  echo "git_head=$(git -C "${REPO}" rev-parse HEAD)"
  echo "selection=abs_delta_y_ge_2p0_abs_delta_phi_lt_pi_over_2"
  echo "sps_leakage_fraction=0"
  echo "dps_template=fixed_seed_exact_half_DPS_reconstruction_MC"
  echo "seed=20260902"
  echo "input_entries=2584"
  echo "model_sha256=$(sha256sum "${MODEL}" | awk '{print $1}')"
  echo "weightdata_sha256=$(sha256sum "${DATA}" | awk '{print $1}')"
  echo "half_dps_mc_sha256=$(sha256sum "${HALF}" | awk '{print $1}')"
  echo "nominal_summary_sha256=$(sha256sum "${NOMINAL_SUMMARY}" | awk '{print $1}')"
  echo "fit_pp_source_sha256=$(sha256sum "${REPO}/Data_driven/fit_pp_2d_adaptive.cpp" | awk '{print $1}')"
  echo "template_source_sha256=$(sha256sum "${REPO}/Data_driven/build_2d_templates_adaptive.cpp" | awk '{print $1}')"
  echo "runner_sha256=$(sha256sum "${REPO}/Data_driven/run_cr_dy2p0_check.sh" | awk '{print $1}')"
  echo "nominal_replacement_status=not_replaced_pending_user_confirmation"
} > "${OUTPUT}/run.metadata.txt"
printf 'candidate_complete_pending_user_confirmation\n' > "${OUTPUT}/complete.marker"
find "${OUTPUT}" -maxdepth 4 -type f ! -name artifact_checksums.txt -print0 \
  | sort -z | xargs -0 sha256sum > "${OUTPUT}/artifact_checksums.txt"
echo "CR_DY2P0_CHECK_COMPLETE tag=${OUTPUT_TAG}"
