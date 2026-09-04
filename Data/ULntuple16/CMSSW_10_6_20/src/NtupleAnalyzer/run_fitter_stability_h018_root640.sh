#!/usr/bin/env bash
set -euo pipefail

run_tag="${1:-h018_fitter_stability_accmix23_effmix19_sps0p85_dps0p15_root640_v8_20260903}"
run_mode="${2:-fresh}"
run_scope="${3:-full}"
differential_reference_choice="${4:-auto}"
if [[ ! "${run_tag}" =~ ^[A-Za-z0-9._-]+$ ]]; then
    echo "Invalid run tag: ${run_tag}" >&2
    exit 2
fi
if [[ "${run_mode}" != "fresh" && "${run_mode}" != "resume" ]]; then
    echo "Run mode must be fresh or resume: ${run_mode}" >&2
    exit 2
fi
if [[ "${run_scope}" != "full" && "${run_scope}" != "total-only" && \
      "${run_scope}" != "collect-only" ]]; then
    echo "Run scope must be full, total-only, or collect-only: ${run_scope}" >&2
    exit 2
fi
if [[ ! "${differential_reference_choice}" =~ ^[A-Za-z0-9._-]+$ ]]; then
    echo "Invalid differential reference choice: ${differential_reference_choice}" >&2
    exit 2
fi

repository="/eos/home-l/leyao/26JJ/JPsiJPsi"
analysis="${repository}/Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer"
nominal_total_tag="newdata167_accmix23_effmix19_sps0p85_dps0p15_lumi36p684_pairsym_cseed50_root640_strategy2_nominal_v4_20260902"
nominal_diff_tag="newdata167_accmix23_effmix19_sps0p85_dps0p15_lumi36p684_pairsym_cseed50_root640_strategy2_diff35_v1_20260902"
weight_name="WeightData_newdata167_accmix23_effmix19_sps0p85_dps0p15_pairsym_cseed50_v1_20260902.root"
nominal_total_dir="${analysis}/fit_results/${nominal_total_tag}"
nominal_diff_dir="${analysis}/fit_results/${nominal_diff_tag}"
weight_data="${nominal_total_dir}/inputs/${weight_name}"
nominal_model="${nominal_total_dir}/total/Model_4D_tot.root"
nominal_fit="${nominal_total_dir}/total/Fit_4D_tot_native_asymptotic_unseeded.root"
nominal_csv="${nominal_diff_dir}/cross_sections.csv"
result_dir="${analysis}/fit_results/${run_tag}"
answer_dir="${repository}/docs/answers/H018_recalculate_fitter_stability_with_legacy_variations"
base_total_source="Fit_4D_tot.cpp"
diff_source="Fit_4D_diff.cpp"
plot_header="Plot_4D.hpp"
generator="make_fitter_stability_h018_source.py"
collector="collect_fitter_stability_h018.py"

if [[ "$(root-config --version)" != "6.40.02" ]]; then
    echo "H018 requires ROOT 6.40.02" >&2
    exit 3
fi
for required in \
    "${weight_data}" "${nominal_model}" "${nominal_fit}" "${nominal_csv}" \
    "${analysis}/${base_total_source}" "${analysis}/${diff_source}" \
    "${analysis}/${plot_header}" "${analysis}/${generator}" \
    "${analysis}/${collector}"; do
    test -s "${required}"
done
if [[ "${run_mode}" == "fresh" ]]; then
    if [[ -e "${result_dir}" || -e "${answer_dir}" ]]; then
        echo "Refusing to reuse an H018 result or answer directory" >&2
        exit 17
    fi
else
    if [[ ! -d "${result_dir}" ]]; then
        echo "H018 resume requires an existing partial result" >&2
        exit 17
    fi
    if [[ -e "${answer_dir}" ]]; then
        if [[ "${run_scope}" != "collect-only" || \
              -n "$(find "${answer_dir}" -mindepth 1 -print -quit)" ]]; then
            echo "H018 resume found a non-empty or unexpected answer directory" >&2
            exit 17
        fi
    fi
fi

[[ "$(sha256sum "${weight_data}" | awk '{print $1}')" == \
    "bc73bd9411c4e84ce343e183a68767377280e905cb042ee6bca61580af6ad7ff" ]]
[[ "$(sha256sum "${nominal_model}" | awk '{print $1}')" == \
    "822b924c25b6f97121e6bf66adb708d41b4a2acbb08737a0718bdb03b4d00b7b" ]]
[[ "$(sha256sum "${nominal_fit}" | awk '{print $1}')" == \
    "27cd0530309dd8fb501586c9ecbf84319372ce5e3a4b47295cdd08635a54b0c6" ]]
[[ "$(sha256sum "${nominal_csv}" | awk '{print $1}')" == \
    "cd1461abb266d3d66b2cf1407b7a67b8500ca4557a2496a11540f05f6a9ce3bc" ]]

mkdir -p "${result_dir}/code" "${result_dir}/inputs" \
    "${result_dir}/total" "${result_dir}/differential" "${result_dir}/logs"
cp "${analysis}/${base_total_source}" "${analysis}/${diff_source}" \
    "${analysis}/${plot_header}" "${analysis}/${generator}" \
    "${analysis}/${collector}" "${result_dir}/code/"
if [[ "${run_mode}" == "fresh" ]]; then
    ln -s "${weight_data}" "${result_dir}/inputs/${weight_name}"
    ln -s "${nominal_model}" "${result_dir}/inputs/Model_4D_tot.root"
    ln -s "${nominal_fit}" "${result_dir}/inputs/Fit_4D_tot_native_asymptotic_unseeded.root"
    ln -s "${nominal_csv}" "${result_dir}/inputs/h017_cross_sections.csv"
else
    test -s "${result_dir}/inputs/${weight_name}"
    test -s "${result_dir}/inputs/Model_4D_tot.root"
    test -s "${result_dir}/inputs/Fit_4D_tot_native_asymptotic_unseeded.root"
    test -s "${result_dir}/inputs/h017_cross_sections.csv"
fi

start_epoch="$(date +%s)"
variations=(
    f1_double_gaussian
    f1_double_crystal_ball
    f1_triple_gaussian
    f2_linear
    f3_single_gaussian
    f3_triple_gaussian
    f4_exp_double_gaussian
    f4_double_exp_gaussian
)

run_total_attempt() {
    local variation="$1"
    local is_reference="$2"
    local offset_enabled="$3"
    local attempt_name="$4"
    local output_dir="$5"
    local fit_result_name="$6"
    local model_name="$7"
    local fix_nx1="$8"
    local relax_prefit_gate="$9"
    local allow_post_covariance_edm="${10}"
    local fix_sigma4="${11:-false}"
    local fix_sigma5="${12:-false}"
    local fix_coef2="${13:-false}"
    local attempt_dir="${output_dir}/${attempt_name}"
    local log_file="${result_dir}/logs/total_${variation}_${attempt_name}.log"
    local raw_gate_row
    local optimizer_edm
    local yield_value
    local yield_error
    local fit_status
    local covariance_quality
    local post_asymptotic_edm
    local selected_entries

    mkdir -p "${attempt_dir}"
    python3 "${analysis}/${generator}" "${variation}" \
        "${analysis}/${base_total_source}" "${attempt_dir}/Fit_4D_tot.cpp"
    cp "${analysis}/${plot_header}" "${attempt_dir}/"
    ln -s "${weight_data}" "${attempt_dir}/WeightData.root"
    if [[ "${offset_enabled}" == "true" ]]; then
        perl -0pi -e \
            's/Strategy\(minimizerStrategy\)\)\;/Strategy(minimizerStrategy), Offset(kTRUE));/g' \
            "${attempt_dir}/Fit_4D_tot.cpp"
    fi
    if [[ "${fix_nx1}" == "true" ]]; then
        perl -0pi -e \
            's@(    RooRealVar Jpsi_nx1\([^\n]+\);\n)@$1    Jpsi_nx1.setVal(50.0);\n    Jpsi_nx1.setConstant(kTRUE);\n@ or die "Jpsi_nx1 declaration not found\n"' \
            "${attempt_dir}/Fit_4D_tot.cpp"
    fi
    if [[ "${fix_sigma4}" == "true" ]]; then
        if [[ "${variation}" != "f3_triple_gaussian" ]]; then
            echo "Jpsi_sigma4 fallback is only defined for f3_triple_gaussian" >&2
            return 1
        fi
        perl -0pi -e \
            's@(    RooRealVar Jpsi_sigma4\([^\n]+\);\n)@$1    Jpsi_sigma4.setVal(0.01);\n    Jpsi_sigma4.setConstant(kTRUE);\n@ or die "Jpsi_sigma4 declaration not found\n"' \
            "${attempt_dir}/Fit_4D_tot.cpp"
    fi
    if [[ "${fix_sigma5}" == "true" ]]; then
        if [[ "${variation}" != "f4_exp_double_gaussian" ]]; then
            echo "Jpsi_sigma5 fallback is only defined for f4_exp_double_gaussian" >&2
            return 1
        fi
        perl -0pi -e \
            's@(    RooRealVar Jpsi_sigma5\([^\n]+\);\n)@$1    Jpsi_sigma5.setVal(0.02);\n    Jpsi_sigma5.setConstant(kTRUE);\n@ or die "Jpsi_sigma5 declaration not found\n"' \
            "${attempt_dir}/Fit_4D_tot.cpp"
    fi
    if [[ "${fix_coef2}" == "true" ]]; then
        if [[ "${variation}" != "f4_double_exp_gaussian" ]]; then
            echo "Jpsi_coef2 fallback is only defined for f4_double_exp_gaussian" >&2
            return 1
        fi
        perl -0pi -e \
            's@(    RooRealVar Jpsi_coef2\([^\n]+\);\n)@$1    Jpsi_coef2.setVal(0.02);\n    Jpsi_coef2.setConstant(kTRUE);\n@ or die "Jpsi_coef2 declaration not found\n"' \
            "${attempt_dir}/Fit_4D_tot.cpp"
    fi
    if [[ "${relax_prefit_gate}" == "true" ]]; then
        perl -0pi -e \
            's@if\(prefitResult && !prefitResult->status\(\) &&\n               prefitResult->edm\(\) < 0\.01\) break;@if(prefitResult && prefitResult->covQual() == 3 &&\n               std::isfinite(prefitResult->edm()) &&\n               prefitResult->edm() < 0.01) break;@ or die "prefit accept gate not found\n"; s@if\(!prefitResult \|\| prefitResult->status\(\) \|\|\n           prefitResult->edm\(\) >= 0\.01\)@if(!prefitResult || prefitResult->covQual() != 3 ||\n           !std::isfinite(prefitResult->edm()) ||\n           prefitResult->edm() >= 0.01)@ or die "prefit failure gate not found\n"' \
            "${attempt_dir}/Fit_4D_tot.cpp"
    fi
    if [[ "${allow_post_covariance_edm}" == "true" ]]; then
        perl -0pi -e \
            's@if\(res && !res->status\(\) && res->edm\(\)<0\.01\) break;@if(res && !res->status() && res->covQual() == 3) break;@ or die "final accept gate not found\n"; s@if\(!res \|\| res->status\(\) \|\| res->edm\(\)>=0\.01\)@if(!res || res->status() || res->covQual() != 3)@ or die "final failure gate not found\n"' \
            "${attempt_dir}/Fit_4D_tot.cpp"
    fi
    cd "${attempt_dir}"
    if ! root -l -b -q \
        "Fit_4D_tot.cpp(${is_reference},false,false,true,false,2)" \
        > "${log_file}" 2>&1; then
        return 1
    fi
    if [[ ! -s "${fit_result_name}" || ! -s "${model_name}" ]]; then
        return 1
    fi
    if ! raw_gate_row="$(root -l -b -q -e \
        "TFile f(\"${fit_result_name}\",\"READ\"); auto r=dynamic_cast<RooFitResult *>(f.Get(\"fit_result_native_asymptotic\")); auto n=r ? dynamic_cast<RooRealVar *>(r->floatParsFinal().find(\"n_P_P\")) : nullptr; if(f.IsZombie() || !r || !n || r->status()!=0 || r->covQual()!=3 || !std::isfinite(r->edm())) gSystem->Exit(5); printf(\"H018ROW,%.17g,%.17g,%d,%d,%.17g,2456\\n\",n->getVal(),n->getError(),r->status(),r->covQual(),r->edm()); gSystem->Exit(0);" \
        2>/dev/null | sed -n 's/^H018ROW,//p')"; then
        return 1
    fi
    if [[ -z "${raw_gate_row}" ]]; then
        return 1
    fi
    optimizer_edm="$(awk -v allow="${allow_post_covariance_edm}" '
        /^Edm   = / { last_edm = $3 }
        /Calculating covariance matrix according/ { before_covariance = last_edm }
        /Fit attempt [0-9]+:.*status=0, covQual=3/ {
            post_edm = $NF
            sub(/^EDM=/, "", post_edm)
            if (allow == "true" || post_edm + 0 < 0.01) {
                print before_covariance
                exit
            }
        }
    ' "${log_file}")"
    if [[ -z "${optimizer_edm}" ]] ||
        ! awk -v edm="${optimizer_edm}" 'BEGIN { exit !(edm + 0 >= 0 && edm + 0 < 0.01) }'; then
        return 1
    fi
    IFS=, read -r yield_value yield_error fit_status covariance_quality \
        post_asymptotic_edm selected_entries <<< "${raw_gate_row}"
    if [[ "${allow_post_covariance_edm}" == "true" ]]; then
        if ! awk -v edm="${post_asymptotic_edm}" 'BEGIN { exit !(edm + 0 >= 0.01) }'; then
            return 1
        fi
    elif ! awk -v edm="${post_asymptotic_edm}" \
        'BEGIN { exit !(edm + 0 >= 0 && edm + 0 < 0.01) }'; then
        return 1
    fi
    printf '%s,%s,%s,%s,%s,%s,%s\n' \
        "${yield_value}" "${yield_error}" "${fit_status}" \
        "${covariance_quality}" "${optimizer_edm}" \
        "${post_asymptotic_edm}" "${selected_entries}" > gate_row.csv
}

run_sigma5_boundary_sequence() {
    local variation="$1"
    local is_reference="$2"
    local output_dir="$3"
    local fit_result_name="$4"
    local model_name="$5"
    local fix_nx1="$6"
    local attempt_key="sigma5_fixed"
    local recovery_key="sigma5_fixed_0p02"
    if [[ "${fix_nx1}" == "true" ]]; then
        attempt_key="nx1_sigma5_fixed"
        recovery_key="nx1_fixed_50_sigma5_fixed_0p02"
    fi

    recovery="${recovery_key}"
    accepted_attempt="attempt_default_${attempt_key}"
    if run_total_attempt "${variation}" "${is_reference}" false \
        "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" \
        "${fix_nx1}" false false false true; then
        return 0
    fi
    recovery="offset_true_${recovery_key}"
    accepted_attempt="attempt_offset_${attempt_key}"
    if run_total_attempt "${variation}" "${is_reference}" true \
        "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" \
        "${fix_nx1}" false false false true; then
        return 0
    fi
    recovery="${recovery_key}_prefit_relaxed"
    accepted_attempt="attempt_default_${attempt_key}_prefit_relaxed"
    if run_total_attempt "${variation}" "${is_reference}" false \
        "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" \
        "${fix_nx1}" true false false true; then
        return 0
    fi
    recovery="offset_true_${recovery_key}_prefit_relaxed"
    accepted_attempt="attempt_offset_${attempt_key}_prefit_relaxed"
    if run_total_attempt "${variation}" "${is_reference}" true \
        "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" \
        "${fix_nx1}" true false false true; then
        return 0
    fi
    recovery="${recovery_key}_prefit_relaxed_post_covariance_edm"
    accepted_attempt="attempt_default_${attempt_key}_prefit_relaxed_post_covariance_edm"
    if run_total_attempt "${variation}" "${is_reference}" false \
        "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" \
        "${fix_nx1}" true true false true; then
        return 0
    fi
    recovery="offset_true_${recovery_key}_prefit_relaxed_post_covariance_edm"
    accepted_attempt="attempt_offset_${attempt_key}_prefit_relaxed_post_covariance_edm"
    run_total_attempt "${variation}" "${is_reference}" true \
        "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" \
        "${fix_nx1}" true true false true
}

recover_sigma5_boundary() {
    local variation="$1"
    local is_reference="$2"
    local output_dir="$3"
    local fit_result_name="$4"
    local model_name="$5"
    if run_sigma5_boundary_sequence "${variation}" "${is_reference}" \
        "${output_dir}" "${fit_result_name}" "${model_name}" false; then
        return 0
    fi
    if rg -q -g "total_${variation}_attempt_*sigma5_fixed*.log" \
        'Jpsi_nx1.*50|outside.*Jpsi_nx1|variable "Jpsi_nx1"' \
        "${result_dir}/logs"; then
        run_sigma5_boundary_sequence "${variation}" "${is_reference}" \
            "${output_dir}" "${fit_result_name}" "${model_name}" true
        return $?
    fi
    return 1
}

run_coef2_boundary_sequence() {
    local variation="$1"
    local is_reference="$2"
    local output_dir="$3"
    local fit_result_name="$4"
    local model_name="$5"
    local fix_nx1="$6"
    local attempt_key="coef2_fixed"
    local recovery_key="coef2_fixed_0p02"
    if [[ "${fix_nx1}" == "true" ]]; then
        attempt_key="nx1_coef2_fixed"
        recovery_key="nx1_fixed_50_coef2_fixed_0p02"
    fi

    recovery="${recovery_key}"
    accepted_attempt="attempt_default_${attempt_key}"
    if run_total_attempt "${variation}" "${is_reference}" false \
        "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" \
        "${fix_nx1}" false false false false true; then
        return 0
    fi
    recovery="offset_true_${recovery_key}"
    accepted_attempt="attempt_offset_${attempt_key}"
    if run_total_attempt "${variation}" "${is_reference}" true \
        "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" \
        "${fix_nx1}" false false false false true; then
        return 0
    fi
    recovery="${recovery_key}_prefit_relaxed"
    accepted_attempt="attempt_default_${attempt_key}_prefit_relaxed"
    if run_total_attempt "${variation}" "${is_reference}" false \
        "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" \
        "${fix_nx1}" true false false false true; then
        return 0
    fi
    recovery="offset_true_${recovery_key}_prefit_relaxed"
    accepted_attempt="attempt_offset_${attempt_key}_prefit_relaxed"
    if run_total_attempt "${variation}" "${is_reference}" true \
        "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" \
        "${fix_nx1}" true false false false true; then
        return 0
    fi
    recovery="${recovery_key}_prefit_relaxed_post_covariance_edm"
    accepted_attempt="attempt_default_${attempt_key}_prefit_relaxed_post_covariance_edm"
    if run_total_attempt "${variation}" "${is_reference}" false \
        "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" \
        "${fix_nx1}" true true false false true; then
        return 0
    fi
    recovery="offset_true_${recovery_key}_prefit_relaxed_post_covariance_edm"
    accepted_attempt="attempt_offset_${attempt_key}_prefit_relaxed_post_covariance_edm"
    run_total_attempt "${variation}" "${is_reference}" true \
        "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" \
        "${fix_nx1}" true true false false true
}

recover_coef2_boundary() {
    local variation="$1"
    local is_reference="$2"
    local output_dir="$3"
    local fit_result_name="$4"
    local model_name="$5"
    if run_coef2_boundary_sequence "${variation}" "${is_reference}" \
        "${output_dir}" "${fit_result_name}" "${model_name}" false; then
        return 0
    fi
    if rg -q -g "total_${variation}_attempt_*coef2_fixed*.log" \
        'Jpsi_nx1.*50|outside.*Jpsi_nx1|variable "Jpsi_nx1"' \
        "${result_dir}/logs"; then
        run_coef2_boundary_sequence "${variation}" "${is_reference}" \
            "${output_dir}" "${fit_result_name}" "${model_name}" true
        return $?
    fi
    return 1
}

for variation in "${variations[@]}"; do
    output_dir="${result_dir}/total/${variation}"
    mkdir -p "${output_dir}"
    if [[ "${run_mode}" == "resume" && -s "${output_dir}/summary.csv" ]]; then
        [[ "$(wc -l < "${output_dir}/summary.csv")" -eq 2 ]]
        continue
    fi
    if [[ "${variation}" == "f1_double_gaussian" ]]; then
        is_reference=true
        fit_result_name="Fit_4D_tot_ref_native_asymptotic_unseeded.root"
        model_name="Model_4D_tot_ref.root"
    else
        is_reference=false
        fit_result_name="Fit_4D_tot_native_asymptotic_unseeded.root"
        model_name="Model_4D_tot.root"
    fi
    recovery="none"
    accepted_attempt="attempt_default"
    skip_standard_attempts=false
    if [[ "${run_mode}" == "resume" && \
          "${variation}" == "f4_exp_double_gaussian" ]] && \
        rg -q -g "total_${variation}_attempt_*.log" \
            'outside.*Jpsi_sigma5|variable "Jpsi_sigma5"' \
            "${result_dir}/logs"; then
        if ! recover_sigma5_boundary "${variation}" "${is_reference}" \
            "${output_dir}" "${fit_result_name}" "${model_name}"; then
            echo "H018 total variation failed after approved sigma5 and nx1 recoveries: ${variation}" >&2
            exit 6
        fi
        skip_standard_attempts=true
    elif [[ "${run_mode}" == "resume" && \
            "${variation}" == "f4_double_exp_gaussian" ]] && \
        rg -q -g "total_${variation}_attempt_*.log" \
            'outside.*Jpsi_coef2|variable "Jpsi_coef2"' \
            "${result_dir}/logs"; then
        if ! recover_coef2_boundary "${variation}" "${is_reference}" \
            "${output_dir}" "${fit_result_name}" "${model_name}"; then
            echo "H018 total variation failed after approved coef2 and nx1 recoveries: ${variation}" >&2
            exit 6
        fi
        skip_standard_attempts=true
    fi
    if [[ "${skip_standard_attempts}" == "false" ]] && \
        ! run_total_attempt "${variation}" "${is_reference}" false \
        "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" false false false; then
        recovery="offset_true"
        accepted_attempt="attempt_offset"
        if ! run_total_attempt "${variation}" "${is_reference}" true \
            "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" false false false; then
            if rg -q 'Jpsi_nx1.*50|outside.*Jpsi_nx1|variable "Jpsi_nx1"' \
                "${result_dir}/logs/total_${variation}_attempt_offset.log"; then
                recovery="nx1_fixed_50"
                accepted_attempt="attempt_default_nx1_fixed"
                if ! run_total_attempt "${variation}" "${is_reference}" false \
                    "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" true false false; then
                    recovery="offset_true_nx1_fixed_50"
                    accepted_attempt="attempt_offset_nx1_fixed"
                    if ! run_total_attempt "${variation}" "${is_reference}" true \
                        "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" true false false; then
                        recovery="nx1_fixed_50_prefit_relaxed"
                        accepted_attempt="attempt_default_nx1_fixed_prefit_relaxed"
                        if ! run_total_attempt "${variation}" "${is_reference}" false \
                            "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" true true false; then
                            recovery="offset_true_nx1_fixed_50_prefit_relaxed"
                            accepted_attempt="attempt_offset_nx1_fixed_prefit_relaxed"
                            if ! run_total_attempt "${variation}" "${is_reference}" true \
                                "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" true true false; then
                                recovery="nx1_fixed_50_prefit_relaxed_post_covariance_edm"
                                accepted_attempt="attempt_default_nx1_fixed_prefit_relaxed_post_covariance_edm"
                                if ! run_total_attempt "${variation}" "${is_reference}" false \
                                    "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" true true true; then
                                    recovery="offset_true_nx1_fixed_50_prefit_relaxed_post_covariance_edm"
                                    accepted_attempt="attempt_offset_nx1_fixed_prefit_relaxed_post_covariance_edm"
                                    if ! run_total_attempt "${variation}" "${is_reference}" true \
                                        "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" true true true; then
                                        if rg -q 'outside.*Jpsi_sigma4|variable "Jpsi_sigma4"' \
                                            "${result_dir}/logs/total_${variation}_attempt_default_nx1_fixed_prefit_relaxed.log" \
                                            "${result_dir}/logs/total_${variation}_attempt_default_nx1_fixed_prefit_relaxed_post_covariance_edm.log"; then
                                            recovery="nx1_fixed_50_sigma4_fixed_0p01"
                                            accepted_attempt="attempt_default_nx1_sigma4_fixed"
                                            if ! run_total_attempt "${variation}" "${is_reference}" false \
                                                "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" true false false true; then
                                                recovery="offset_true_nx1_fixed_50_sigma4_fixed_0p01"
                                                accepted_attempt="attempt_offset_nx1_sigma4_fixed"
                                                if ! run_total_attempt "${variation}" "${is_reference}" true \
                                                    "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" true false false true; then
                                                    recovery="nx1_fixed_50_sigma4_fixed_0p01_prefit_relaxed"
                                                    accepted_attempt="attempt_default_nx1_sigma4_fixed_prefit_relaxed"
                                                    if ! run_total_attempt "${variation}" "${is_reference}" false \
                                                        "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" true true false true; then
                                                        recovery="offset_true_nx1_fixed_50_sigma4_fixed_0p01_prefit_relaxed"
                                                        accepted_attempt="attempt_offset_nx1_sigma4_fixed_prefit_relaxed"
                                                        if ! run_total_attempt "${variation}" "${is_reference}" true \
                                                            "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" true true false true; then
                                                            recovery="nx1_fixed_50_sigma4_fixed_0p01_prefit_relaxed_post_covariance_edm"
                                                            accepted_attempt="attempt_default_nx1_sigma4_fixed_prefit_relaxed_post_covariance_edm"
                                                            if ! run_total_attempt "${variation}" "${is_reference}" false \
                                                                "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" true true true true; then
                                                                recovery="offset_true_nx1_fixed_50_sigma4_fixed_0p01_prefit_relaxed_post_covariance_edm"
                                                                accepted_attempt="attempt_offset_nx1_sigma4_fixed_prefit_relaxed_post_covariance_edm"
                                                                if ! run_total_attempt "${variation}" "${is_reference}" true \
                                                                    "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" true true true true; then
                                                                    echo "H018 total variation failed after approved nx1, sigma4, prefit, and post-covariance EDM recoveries: ${variation}" >&2
                                                                    exit 6
                                                                fi
                                                            fi
                                                        fi
                                                    fi
                                                fi
                                            fi
                                        else
                                            echo "H018 total variation failed after approved nx1, prefit, and post-covariance EDM recoveries: ${variation}" >&2
                                            exit 6
                                        fi
                                    fi
                                fi
                            fi
                        fi
                    fi
                fi
            else
                recovery="prefit_relaxed"
                accepted_attempt="attempt_default_prefit_relaxed"
                if ! run_total_attempt "${variation}" "${is_reference}" false \
                    "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" false true false; then
                    recovery="offset_true_prefit_relaxed"
                    accepted_attempt="attempt_offset_prefit_relaxed"
                    if ! run_total_attempt "${variation}" "${is_reference}" true \
                        "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" false true false; then
                        recovery="prefit_relaxed_post_covariance_edm"
                        accepted_attempt="attempt_default_prefit_relaxed_post_covariance_edm"
                        if ! run_total_attempt "${variation}" "${is_reference}" false \
                            "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" false true true; then
                            recovery="offset_true_prefit_relaxed_post_covariance_edm"
                            accepted_attempt="attempt_offset_prefit_relaxed_post_covariance_edm"
                            if ! run_total_attempt "${variation}" "${is_reference}" true \
                                "${accepted_attempt}" "${output_dir}" "${fit_result_name}" "${model_name}" false true true; then
                                if [[ "${variation}" == "f4_exp_double_gaussian" ]] && \
                                    rg -q -g "total_${variation}_attempt_*.log" \
                                        'outside.*Jpsi_sigma5|variable "Jpsi_sigma5"' \
                                        "${result_dir}/logs"; then
                                    if ! recover_sigma5_boundary "${variation}" "${is_reference}" \
                                        "${output_dir}" "${fit_result_name}" "${model_name}"; then
                                        echo "H018 total variation failed after approved sigma5 and nx1 recoveries: ${variation}" >&2
                                        exit 6
                                    fi
                                elif [[ "${variation}" == "f4_double_exp_gaussian" ]] && \
                                    rg -q -g "total_${variation}_attempt_*.log" \
                                        'outside.*Jpsi_coef2|variable "Jpsi_coef2"' \
                                        "${result_dir}/logs"; then
                                    if ! recover_coef2_boundary "${variation}" "${is_reference}" \
                                        "${output_dir}" "${fit_result_name}" "${model_name}"; then
                                        echo "H018 total variation failed after approved coef2 and nx1 recoveries: ${variation}" >&2
                                        exit 6
                                    fi
                                else
                                    echo "H018 total variation failed after approved prefit and post-covariance EDM recoveries: ${variation}" >&2
                                    exit 6
                                fi
                            fi
                        fi
                    fi
                fi
            fi
        fi
    fi
    cd "${output_dir}"
    ln -s "${accepted_attempt}/${fit_result_name}" "${fit_result_name}"
    ln -s "${accepted_attempt}/${model_name}" "${model_name}"
    {
        echo "variation,yield,yield_error,status,covQual,optimizer_edm,post_asymptotic_edm,selected_entries,recovery"
        printf '%s,' "${variation}"
        tr -d '\n' < "${accepted_attempt}/gate_row.csv"
        printf ',%s\n' "${recovery}"
    } > summary.csv
    [[ "$(wc -l < summary.csv)" -eq 2 ]]
done

summary_files=()
for variation in "${variations[@]}"; do
    summary_files+=("${result_dir}/total/${variation}/summary.csv")
    [[ "$(wc -l < "${result_dir}/total/${variation}/summary.csv")" -eq 2 ]]
done
if [[ "${run_scope}" == "total-only" ]]; then
    end_epoch="$(date +%s)"
    {
        echo "status=total_complete_pending_user_differential_reference"
        echo "request=H018_recalculate_fitter_stability_with_legacy_variations"
        echo "run_tag=${run_tag}"
        echo "run_mode=${run_mode}"
        echo "run_scope=${run_scope}"
        echo "completed_total_variations=${#variations[@]}"
        echo "root_version=$(root-config --version)"
        echo "resume_start_epoch=${start_epoch}"
        echo "stage_end_epoch=${end_epoch}"
        echo "resume_wall_seconds=$((end_epoch - start_epoch))"
        echo "differential_status=not_started_pending_user_reference_choice"
    } > "${result_dir}/total_stage.metadata.txt"
    echo "H018_TOTAL_COMPLETE_PENDING_REFERENCE result=${result_dir}"
    exit 0
fi
envelope_variation="$(awk -F, -v nominal=5687.8051561631491 '
    FNR == 2 {
        relative_shift = ($2 - nominal) / nominal
        if (relative_shift < 0) relative_shift = -relative_shift
        if (relative_shift > maximum || selected == "") {
            maximum = relative_shift
            selected = $1
        }
    }
    END { print selected }
' "${summary_files[@]}")"
if [[ -z "${envelope_variation}" ]]; then
    echo "Cannot determine the total fitter-envelope variation" >&2
    exit 7
fi
if [[ "${differential_reference_choice}" == "auto" ]]; then
    differential_reference="${envelope_variation}"
    differential_reference_selection="automatic_total_envelope"
else
    differential_reference="${differential_reference_choice}"
    differential_reference_selection="user_confirmed"
fi
if [[ "${differential_reference}" != "${envelope_variation}" ]]; then
    echo "Chosen differential reference is not the largest total variation: " \
         "${differential_reference} != ${envelope_variation}" >&2
    exit 7
fi
if [[ "${differential_reference}" == "f1_double_gaussian" ]]; then
    reference_model="${result_dir}/total/${differential_reference}/Model_4D_tot_ref.root"
else
    reference_model="${result_dir}/total/${differential_reference}/Model_4D_tot.root"
fi
test -s "${reference_model}"
reference_model_directory="${result_dir}/differential_reference"
reference_model_alias="${reference_model_directory}/Model_4D_tot_ref.root"
mkdir -p "${reference_model_directory}"
if [[ -e "${reference_model_alias}" || -L "${reference_model_alias}" ]]; then
    [[ "$(readlink -f "${reference_model_alias}")" == \
       "$(readlink -f "${reference_model}")" ]]
else
    ln -s "${reference_model}" "${reference_model_alias}"
fi
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
if [[ "${run_scope}" == "full" ]]; then
    cd "${result_dir}/code"
    root -l -b -q -e \
        'if(!gSystem->CompileMacro("Fit_4D_diff.cpp","kO")) gSystem->Exit(1);' \
        > "${result_dir}/logs/compile_diff.log" 2>&1
    diff_library="${result_dir}/code/Fit_4D_diff_cpp.so"
    test -s "${diff_library}"
    for index in "${!bins[@]}"; do
        read -r variable minimum maximum <<< "${bins[${index}]}"
        printf -v suffix '%02d' "${index}"
        cd "${analysis}"
        root -l -b -q -e \
            "gSystem->Load(\"${diff_library}\"); Fit_4D_diff(\"${variable}\",${minimum},${maximum},true,true,\"${run_tag}/differential\",\"${weight_data}\",\"${reference_model_directory}\",false,true);" \
            > "${result_dir}/logs/diff_${suffix}_${variable}.log" 2>&1
        relative_result="$(sed -n 's/^Accepted differential result: \([^,]*\),.*/\1/p' "${result_dir}/logs/diff_${suffix}_${variable}.log" | tail -n 1)"
        test -n "${relative_result}"
        output_root="${analysis}/${relative_result}"
        test -s "${output_root}"
    done
fi

mapfile -t differential_roots < <(
    find "${result_dir}/differential" -maxdepth 1 -type f -name '*.root' | sort
)
[[ "${#differential_roots[@]}" -eq 35 ]]
for output_root in "${differential_roots[@]}"; do
    temporary_csv="${output_root}.csv.tmp"
    {
        echo "variable,bin_min,bin_max,yield,yield_error,status,covQual,edm,selected_entries,n_Sig_Comb_fallback,n_Comb_Comb_fallback,strategy,offset"
        root -l -b -q -e \
            "TFile f(\"${output_root}\",\"READ\"); auto r=dynamic_cast<RooFitResult *>(f.Get(\"fit_result_native_asymptotic\")); auto n=r ? dynamic_cast<RooRealVar *>(r->floatParsFinal().find(\"n_P_P\")) : nullptr; auto v=dynamic_cast<TObjString *>(f.Get(\"variable\")); auto lo=dynamic_cast<TParameter<double> *>(f.Get(\"bin_min\")); auto hi=dynamic_cast<TParameter<double> *>(f.Get(\"bin_max\")); auto e=dynamic_cast<TParameter<int> *>(f.Get(\"selected_entries\")); auto fs=dynamic_cast<TParameter<int> *>(f.Get(\"boundary_fallback_n_Sig_Comb_and_n_Comb_Sig_shared_fixed_zero\")); auto fc=dynamic_cast<TParameter<int> *>(f.Get(\"boundary_fallback_n_Comb_Comb_fixed_zero\")); auto st=dynamic_cast<TParameter<int> *>(f.Get(\"minimizer_strategy\")); auto of=dynamic_cast<TParameter<int> *>(f.Get(\"offset_enabled\")); if(f.IsZombie() || !r || !n || !v || !lo || !hi || !e || !fs || !fc || !st || !of || r->status()!=0 || r->covQual()!=3 || !std::isfinite(r->edm()) || r->edm()>=0.01) gSystem->Exit(5); printf(\"H018ROW,%s,%.17g,%.17g,%.17g,%.17g,%d,%d,%.17g,%d,%d,%d,%d,%d\\n\",v->GetString().Data(),lo->GetVal(),hi->GetVal(),n->getVal(),n->getError(),r->status(),r->covQual(),r->edm(),e->GetVal(),fs->GetVal(),fc->GetVal(),st->GetVal(),of->GetVal()); gSystem->Exit(0);" \
            2>/dev/null | sed -n 's/^H018ROW,//p'
    } > "${temporary_csv}"
    [[ "$(wc -l < "${temporary_csv}")" -eq 2 ]]
    mv "${temporary_csv}" "${output_root}.csv"
done

mkdir -p "${answer_dir}"
python3 "${analysis}/${collector}" \
    --run-dir "${result_dir}" \
    --nominal-csv "${nominal_csv}" \
    --answer-dir "${answer_dir}" \
    --differential-reference "${differential_reference}" \
    > "${result_dir}/logs/collect.log" 2>&1

end_epoch="$(date +%s)"
{
    echo "status=complete"
    echo "request=H018_recalculate_fitter_stability_with_legacy_variations"
    echo "run_tag=${run_tag}"
    echo "run_mode=${run_mode}"
    echo "run_scope=${run_scope}"
    echo "nominal_total_tag=${nominal_total_tag}"
    echo "nominal_diff_tag=${nominal_diff_tag}"
    echo "total_prescription=maximum absolute relative n_P_P shift among eight legacy one-at-a-time PDF variations"
    echo "differential_prescription=per-bin absolute relative n_P_P shift using the total-envelope variation as the single reference"
    echo "differential_reference=${differential_reference}"
    echo "differential_reference_selection=${differential_reference_selection}"
    echo "differential_reference_model=$(readlink -f "${reference_model}")"
    echo "differential_reference_model_sha256=$(sha256sum "${reference_model}" | awk '{print $1}')"
    echo "fit_error=ROOT native AsymptoticError(true), no external seed"
    echo "total_minimizer=Strategy(2); default Offset(false), recover with Offset(true) only after failure; across all eight models, when a present parameter causes a boundary failure, fix Jpsi_nx1=50, Jpsi_sigma4=0.01, Jpsi_sigma5=0.02, or Jpsi_coef2=0.02 at the corresponding boundary and refit; final recovery may accept a discarded central prefit with covQual=3 and optimizer EDM<0.01 regardless of prefit status; a dedicated last recovery accepts status=0,covQual=3 with optimizer EDM<0.01 when native asymptotic covariance recalculation makes the stored post-asymptotic EDM>=0.01; native corrected covariance is retained and both EDM values are recorded"
    echo "differential_recovery=Strategy(1), then Strategy(2)+Offset(true); approved zero-boundary fallback for n_Sig_Comb and n_Comb_Comb"
    echo "fit_gate=status=0,covQual=3,optimizer_EDM<0.01; normally post_asymptotic_EDM<0.01, with an explicit recorded exception only for post-covariance EDM anomaly"
    echo "root_version=$(root-config --version)"
    echo "git_head=$(git -C "${repository}" rev-parse HEAD)"
    echo "git_dirty=$(git -C "${repository}" status --porcelain | wc -l) entries"
    echo "start_epoch=${start_epoch}"
    echo "end_epoch=${end_epoch}"
    echo "wall_seconds=$((end_epoch - start_epoch))"
    echo "weightdata_sha256=$(sha256sum "${weight_data}" | awk '{print $1}')"
    echo "nominal_model_sha256=$(sha256sum "${nominal_model}" | awk '{print $1}')"
    echo "nominal_fit_sha256=$(sha256sum "${nominal_fit}" | awk '{print $1}')"
    echo "nominal_table_sha256=$(sha256sum "${nominal_csv}" | awk '{print $1}')"
    echo "legacy_variation_source=git b6a7526 nopre_fit.cpp"
    echo "legacy_variation_source_sha256=$(git -C "${repository}" show b6a7526:Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/nopre_fit.cpp | sha256sum | awk '{print $1}')"
} > "${result_dir}/run.metadata.txt"

find "${result_dir}" -type f ! -name checksums.sha256 -print0 \
    | sort -z | xargs -0 sha256sum > "${result_dir}/checksums.sha256"
find "${answer_dir}" -type f ! -name checksums.sha256 -print0 \
    | sort -z | xargs -0 sha256sum > "${answer_dir}/checksums.sha256"

echo "H018_COMPLETE result=${result_dir} answer=${answer_dir}"
