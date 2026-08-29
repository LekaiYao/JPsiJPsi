#!/usr/bin/env python3
import csv
import hashlib
import math
import pathlib
import subprocess
import sys


def read_kv(path):
    values = {}
    for line in pathlib.Path(path).read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key] = value
    return values


def sha256(path):
    digest = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def envelope(shifts):
    return max([0.0] + shifts), max([0.0] + [-value for value in shifts])


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: summarize_fdps_global_nominal_systematics.py RESULT_DIR")
    base = pathlib.Path(sys.argv[1])
    repo = pathlib.Path(__file__).resolve().parents[2]
    nominal = float(read_kv(base / "baseline" / "summary.txt")["f_dps"])
    expected = 0.170423232645
    if abs(nominal - expected) > 5e-10:
        raise RuntimeError(f"global nominal f_DPS changed: {nominal} vs {expected}")
    bootstrap_path = (repo / "tests/atlas_data_driven_sps/results/"
                      "route9p3_global_calibrated_allpairs/"
                      "bootstrap_fixed500_global_nominal_v1/summary.txt")
    bootstrap = read_kv(bootstrap_path)
    if bootstrap.get("status") != "complete" or bootstrap.get("failures") != "0":
        raise RuntimeError("global nominal fixed-500 bootstrap is not complete")
    statistical_uncertainty = float(bootstrap["statistical_uncertainty"])

    rows = []
    def add(group, label, value, note):
        rows.append({
            "group": group,
            "variation": label,
            "f_dps": f"{value:.12g}",
            "signed_shift": f"{value - nominal:.12g}",
            "note": note,
        })

    add("baseline", "nominal", nominal,
        "calibrated_global_allpairs_singleCR_12cell_mergeddy")
    category_shifts = []
    for label, note in [
        ("pt15", "two_categories_split_at_pt15"),
        ("absy1p2", "two_categories_split_at_abs_y_1p2"),
    ]:
        value = float(read_kv(base / "category" / label / "summary.txt")["f_dps"])
        add("single_j_category", label, value, note)
        category_shifts.append(value - nominal)

    cr_shifts = []
    for label, note in [
        ("nominal_purity", "purity"),
        ("cr_phi70_zero", "zero_sps_leakage"),
        ("cr_phi70_purity", "mc_sps_leakage"),
        ("cr_phi110_zero", "zero_sps_leakage"),
        ("cr_phi110_purity", "mc_sps_leakage"),
        ("cr_dy2p4_zero", "zero_sps_leakage"),
        ("cr_dy2p4_purity", "mc_sps_leakage"),
    ]:
        value = float(read_kv(base / "cr" / label / "summary.txt")["f_dps"])
        add("cr_definition_and_purity", label, value, note)
        cr_shifts.append(value - nominal)

    closure = read_kv(base / "dps_mc_full_workflow_closure" / "metadata.txt")
    closure_value = float(closure["observed_f_dps"])
    closure_relative = closure_value - 1.0
    closure_abs = float(closure["absolute_f_dps_systematic"])
    add("dps_mc_full_workflow_closure", "pure_dps_pseudodata",
        nominal + closure_abs,
        f"symmetric_abs_uncertainty;closure_f_dps={closure_value:.12g};"
        f"relative_bias={closure_relative:.12g}")

    projection_path = base / "projections_1d" / "summary.txt"
    projections = read_kv(projection_path)
    if (projections.get("status") != "complete" or
            projections.get("included_in_systematic") != "false"):
        raise RuntimeError("1D projection cross-check is incomplete or marked as systematic")
    projection_variables = ["delta_y", "delta_phi", "evt_mass", "evt_y"]
    for variable in projection_variables:
        value = float(projections[f"{variable}_f_dps"])
        add("projection_crosscheck", variable, value,
            "formal_crosscheck_not_in_systematic")

    category_up, category_down = envelope(category_shifts)
    cr_up, cr_down = envelope(cr_shifts)
    total_up = math.sqrt(category_up**2 + cr_up**2 + closure_abs**2)
    total_down = math.sqrt(category_down**2 + cr_down**2 + closure_abs**2)

    with (base / "variations.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    summary = [
        "status=global_nominal_stat_and_systematics_complete",
        f"nominal_f_dps={nominal:.12g}",
        "nominal_definition=calibrated_global_allpairs_singleCR_12cell_mergeddy",
        f"statistical_uncertainty={statistical_uncertainty:.12g}",
        "statistical_uncertainty_status=frozen_fixed500_global_bootstrap",
        f"bootstrap_mean={float(bootstrap['mean']):.12g}",
        f"bootstrap_mean_shift={float(bootstrap['mean_shift_from_nominal']):.12g}",
        f"bootstrap_sd_monte_carlo_precision_approx={float(bootstrap['sd_monte_carlo_precision_approx']):.12g}",
        f"category_group_up={category_up:.12g}",
        f"category_group_down={category_down:.12g}",
        "category_group_approved=true",
        f"cr_group_up={cr_up:.12g}",
        f"cr_group_down={cr_down:.12g}",
        "cr_group_approved=true_definition_recalculated_for_global_nominal",
        f"dps_mc_full_workflow_closure_f_dps={closure_value:.12g}",
        f"dps_mc_full_workflow_closure_relative_bias={closure_relative:.12g}",
        f"dps_mc_full_workflow_closure_symmetric={closure_abs:.12g}",
        "dps_mc_full_workflow_closure_approved_definition_recalculated_for_global_nominal",
        f"combined_method_systematic_up={total_up:.12g}",
        f"combined_method_systematic_down={total_down:.12g}",
        "combination=quadrature_across_category_CR_DPS_MC_closure",
        "binning_formal_systematic=false",
        "projection_crosscheck_status=complete_not_in_systematic",
        *[
            f"projection_{variable}_f_dps={float(projections[f'{variable}_f_dps']):.12g}"
            for variable in projection_variables
        ],
        *[
            f"projection_{variable}_signed_shift={float(projections[f'{variable}_signed_shift']):.12g}"
            for variable in projection_variables
        ],
        "scope=f_DPS_specific_method_systematics_only",
    ]
    (base / "summary.txt").write_text("\n".join(summary) + "\n")

    sources = [
        repo / "tests/atlas_data_driven_sps/fit_single_jpsi_category_splot.cpp",
        repo / "tests/atlas_data_driven_sps/fit_single_jpsi_global_calibrated_splot.cpp",
        repo / "tests/atlas_data_driven_sps/fit_pp_2d_adaptive.cpp",
        repo / "tests/atlas_data_driven_sps/build_2d_templates_adaptive.cpp",
        repo / "tests/atlas_data_driven_sps/run_dps_mc_full_workflow_closure.sh",
        repo / "tests/atlas_data_driven_sps/run_fdps_global_nominal_systematics.sh",
        repo / "tests/atlas_data_driven_sps/run_global_nominal_projection_crosscheck.sh",
        repo / "tests/atlas_data_driven_sps/build_component_trees.cpp",
        repo / "tests/atlas_data_driven_sps/fit_component_variables.cpp",
        repo / "tests/atlas_data_driven_sps/summarize_component_fits.cpp",
        repo / "tests/atlas_data_driven_sps/fit_1d_component_templates.cpp",
        repo / "tests/atlas_data_driven_sps/summarize_fdps_global_nominal_systematics.py",
    ]
    metadata = [
        "artifact=fdps_global_nominal_method_systematics",
        "runtime=CMSSW_10_6_20_ROOT_6_14_09_cmssw_el7",
        f"git_commit={subprocess.check_output(['git','rev-parse','HEAD'], cwd=str(repo), universal_newlines=True).strip()}",
        "nominal_definition=calibrated_global_allpairs_singleCR_12cell_mergeddy",
        f"nominal_f_dps={nominal:.12g}",
        "category_systematic_definition=pt15_and_absy1p2_variations_about_global",
        "category_systematic_approved=true",
        "bootstrap=fixed500_global_complete_500_of_500_no_replacement",
        f"bootstrap_statistical_uncertainty={statistical_uncertainty:.12g}",
        f"sha256_bootstrap_summary={sha256(bootstrap_path)}",
        "projection_crosscheck=complete_four_variables_not_in_systematic",
        f"sha256_projection_summary={sha256(projection_path)}",
    ]
    for source in sources:
        metadata.append(f"sha256_{source.name}={sha256(source)}")
    metadata.append("real_newlines=true")
    (base / "metadata.txt").write_text("\n".join(metadata) + "\n")
    print("\n".join(summary))


if __name__ == "__main__":
    main()
