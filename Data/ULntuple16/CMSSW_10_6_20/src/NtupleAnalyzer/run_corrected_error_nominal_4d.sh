#!/usr/bin/env bash
set -euo pipefail

run_tag="${1:-corrected_error_nominal_root640_unseeded_20260901}"
result_dir="fit_results/${run_tag}"
result_csv="${result_dir}/cross_sections.csv"
metadata_file="${result_dir}/run.metadata.txt"

if [[ -e "${result_csv}" || -e "${metadata_file}" ]]; then
    echo "Refusing to overwrite an existing run: ${result_dir}" >&2
    exit 2
fi

mkdir -p "${result_dir}/logs"

bins=(
    "delta_y 0 0.5"
    "delta_y 0.5 1"
    "delta_y 1 1.5"
    "delta_y 1.5 2"
    "delta_y 2 2.5"
    "delta_y 2.5 4"
    "delta_phi 0 0.3927"
    "delta_phi 0.3927 0.7854"
    "delta_phi 0.7854 1.1781"
    "delta_phi 1.1781 1.5708"
    "delta_phi 1.5708 1.9635"
    "delta_phi 1.9635 2.3562"
    "delta_phi 2.3562 2.7489"
    "delta_phi 2.7489 3.1416"
    "evt_mass 7.5 17.5"
    "evt_mass 17.5 27.5"
    "evt_mass 27.5 37.5"
    "evt_mass 37.5 47.5"
    "evt_mass 47.5 57.5"
    "evt_mass 57.5 67.5"
    "evt_mass 67.5 107.5"
    "evt_y 0 0.4"
    "evt_y 0.4 0.8"
    "evt_y 0.8 1.2"
    "evt_y 1.2 1.6"
    "evt_y 1.6 2.0"
    "evt_pt 0 5"
    "evt_pt 5 10"
    "evt_pt 10 15"
    "evt_pt 15 20"
    "evt_pt 20 25"
    "evt_pt 25 30"
    "evt_pt 30 35"
    "evt_pt 35 40"
    "evt_pt 40 80"
)

printf '%s\n' 'scope,variable,bin_min,bin_max,yield,stat_yield,reference_yield,fitter_relative,sigma_pb_per_unit,stat_pb_per_unit,fitter_pb_per_unit,br_relative,luminosity_relative,correction_relative,lifetime_relative,total_systematic_relative,total_systematic_pb_per_unit,nominal_n_Sig_Comb_n_Comb_Sig_shared_fallback,nominal_n_Comb_Comb_fallback,reference_n_Sig_Comb_n_Comb_Sig_shared_fallback,reference_n_Comb_Comb_fallback,status,covQual,edm,reference_status,reference_covQual,reference_edm' > "${result_csv}"
root -l -b -q 'Fit_valid.cpp("total")' | tail -n 1 >> "${result_csv}"

completed=0
for bin_spec in "${bins[@]}"; do
    read -r variable minimum maximum <<< "${bin_spec}"
    log_stem="${variable}_${minimum//./p}_${maximum//./p}"
    echo "[$((completed + 1))/35] ${variable}: ${minimum} to ${maximum}"
    root -l -b -q "Fit_4D_diff.cpp(\"${variable}\",${minimum},${maximum},false,true,\"${run_tag}\")" > "${result_dir}/logs/${log_stem}_nominal.log" 2>&1
    root -l -b -q "Fit_4D_diff.cpp(\"${variable}\",${minimum},${maximum},true,true,\"${run_tag}\")" > "${result_dir}/logs/${log_stem}_reference.log" 2>&1
    root -l -b -q "Fit_valid.cpp(\"${variable}\",${minimum},${maximum},\"${run_tag}\")" | tail -n 1 >> "${result_csv}"
    completed=$((completed + 1))
done

{
    echo "status=complete"
    echo "run_tag=${run_tag}"
    echo "completed_bins=${completed}"
    echo "error_convention=ROOT native AsymptoticError(true), no external seed"
    echo "fit_gate=status=0,covQual=3,EDM<0.01; maximum 3 attempts per fit"
    echo "boundary_fallback=if the shared n_Sig_Comb/n_Comb_Sig yield or n_Comb_Comb reaches 0 and corrected covariance throws, fix the affected independent yield to 0 and refit; separate flags are stored per bin"
    echo "luminosity_fb=36.31"
    echo "branching_fraction_jpsi_to_mumu=0.05961"
    echo "branching_fraction_systematic_relative=0.011"
    echo "luminosity_systematic_relative=0.012"
    echo "lifetime_systematic_relative=0.003 (historical reused input)"
    echo "total_correction_systematic_relative=0.061 (temporary historical input)"
    echo "differential_correction_systematics=temporary historical sys3 arrays from Tpl_Fit.cpp"
    echo "fitter_systematics=absolute nominal-reference yield differences from this run"
    echo "formula_total_pb=Ncorr*1e-3/(luminosity_fb*branching_fraction^2)"
    echo "formula_differential_pb_per_unit=formula_total_pb/bin_width"
    echo "root_version=$(root-config --version)"
    echo "git_head=$(git rev-parse HEAD)"
    echo "generated_at=$(date --iso-8601=seconds)"
    sha256sum WeightData.root Model_4D_tot.root Model_4D_tot_ref.root \
        Fit_4D_tot_native_asymptotic_unseeded.root \
        Fit_4D_tot_ref_native_asymptotic_unseeded.root \
        Fit_4D_diff.cpp Fit_valid.cpp run_corrected_error_nominal_4d.sh
} > "${metadata_file}"

echo "Completed ${completed} bins: ${result_csv}"
