#!/usr/bin/env python3
"""Propagate the confirmed fixed-half-DPS-MC f_DPS to sigma_DPS and sigma_eff."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from pathlib import Path


def read_kv(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        if key in values:
            raise RuntimeError(f"duplicate key {key} in {path}")
        values[key] = value
    return values


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def quadrature(*values: float) -> float:
    return math.sqrt(sum(value * value for value in values))


def percent(value: float, central: float) -> float:
    return 100.0 * value / central


def sigma_eff(
    f_dps: float,
    pair_pb: float,
    single_visible_nb: float,
    epsilon_m: float,
    branching: float,
) -> float:
    sigma_dps_no_mass_pb = f_dps * pair_pb / epsilon_m
    visible_dps_pb = branching * branching * sigma_dps_no_mass_pb
    return 0.001 * 0.5 * single_visible_nb * single_visible_nb / visible_dps_pb


def main() -> None:
    if len(sys.argv) != 6:
        raise SystemExit(
            "usage: calculate_current_nominal_sigma_eff.py "
            "F_DPS_SUMMARY MASS_SUMMARY PAIR_JSON REFERENCE_INPUTS_JSON OUTPUT_DIR"
        )
    f_path, mass_path, pair_path, reference_path, output_path = map(Path, sys.argv[1:])
    output_path.mkdir(parents=True, exist_ok=True)
    targets = [
        output_path / "inputs.json",
        output_path / "results.json",
        output_path / "results.csv",
        output_path / "uncertainty_components.csv",
        output_path / "summary.txt",
    ]
    if any(path.exists() for path in targets):
        raise RuntimeError(f"refusing to overwrite calculation files in {output_path}")

    f_values = read_kv(f_path)
    mass_values = read_kv(mass_path)
    pair = json.loads(pair_path.read_text(encoding="utf-8"))
    reference = json.loads(reference_path.read_text(encoding="utf-8"))
    f_status = f_values.get("status", "")
    if f_status not in {
        "nominal_complete_user_confirmed",
        "candidate_complete_pending_user_confirmation",
    }:
        raise RuntimeError("f_DPS input is neither a confirmed nominal nor a complete candidate")
    if mass_values.get("status") != "complete":
        raise RuntimeError("mass extrapolation input is incomplete")
    pair_status = pair.get("status", "")
    if pair_status not in {
        "accepted_nominal_user_confirmed",
        "temporary_total_only_pending_Data_cross_section_and_user_confirmation",
    }:
        raise RuntimeError("unsupported pair cross section status")
    pair_pending_sources = pair.get("pending_systematic_sources", {})
    if not isinstance(pair_pending_sources, dict):
        raise RuntimeError("pending pair systematic sources must be a mapping")
    for name, source in pair_pending_sources.items():
        if source.get("relative") is not None:
            raise RuntimeError(f"pending pair source {name} unexpectedly has a value")
        if source.get("status") != "pending_not_included_not_zero":
            raise RuntimeError(f"pending pair source {name} lacks explicit nonzero-pending status")
    result_status = (
        "candidate_complete_included_sources_pending_pair_systematic_and_user_confirmation"
        if pair_pending_sources
        else "candidate_complete_pending_user_confirmation"
    )
    result_scope = (
        "matching_generation_total_only_pair_candidate_included_systematics"
        if pair_pending_sources
        else "newdata_matching_generation_current_confirmed_method"
    )

    f_dps = float(f_values["f_dps"])
    f_stat = float(f_values["stat_uncertainty"])
    f_cr_up = float(f_values["cr_systematic_up"])
    f_cr_down = float(f_values["cr_systematic_down"])
    f_usage = float(f_values["dps_mc_usage_systematic"])
    f_syst_up = quadrature(f_cr_up, f_usage)
    f_syst_down = quadrature(f_cr_down, f_usage)

    pair_pb = float(pair["central_pb"])
    pair_stat = float(pair["stat_pb"])
    pair_sources = {
        name: float(source["relative"])
        for name, source in pair["systematic_sources_relative"].items()
    }
    pair_syst_rel = quadrature(*pair_sources.values())
    if not math.isclose(
        pair_syst_rel, float(pair["reconstructed_relative_systematic"]), rel_tol=1e-12
    ):
        raise RuntimeError("pair systematic sources do not reconstruct the recorded total")
    pair_syst_pb = pair_pb * pair_syst_rel
    pair_non_br_sources = {
        name: value for name, value in pair_sources.items() if name != "branching_fraction"
    }
    pair_non_br_syst_rel = quadrature(*pair_non_br_sources.values())

    n_all = int(mass_values["all_cross_event_pairs"])
    n_pass = int(mass_values["pass_mJJ_gt_7p5"])
    epsilon_m = n_pass / n_all
    if not 0.0 < epsilon_m <= 1.0:
        raise RuntimeError("invalid epsilon_m")

    single = reference["single_jpsi_cross_section"]
    branching_input = reference["branching_fraction"]
    single_nb = float(single["visible_central_nb"])
    single_stat = float(single["visible_stat_nb"])
    single_bin_syst = float(single["visible_bin_syst_nb"])
    single_lumi_rel = float(single["luminosity_relative"])
    branching = float(branching_input["central"])

    sigma_dps_cut = f_dps * pair_pb
    sigma_dps_stat_f = pair_pb * f_stat
    sigma_dps_stat_pair = f_dps * pair_stat
    sigma_dps_stat = quadrature(sigma_dps_stat_f, sigma_dps_stat_pair)
    sigma_dps_syst_f_up = pair_pb * f_syst_up
    sigma_dps_syst_f_down = pair_pb * f_syst_down
    sigma_dps_syst_pair = f_dps * pair_syst_pb
    sigma_dps_syst_up = quadrature(sigma_dps_syst_f_up, sigma_dps_syst_pair)
    sigma_dps_syst_down = quadrature(sigma_dps_syst_f_down, sigma_dps_syst_pair)

    inverse_epsilon = 1.0 / epsilon_m
    sigma_dps_no_mass = sigma_dps_cut * inverse_epsilon
    sigma_dps_no_mass_stat = sigma_dps_stat * inverse_epsilon
    sigma_dps_no_mass_syst_up = sigma_dps_syst_up * inverse_epsilon
    sigma_dps_no_mass_syst_down = sigma_dps_syst_down * inverse_epsilon

    sigma_eff_central = sigma_eff(f_dps, pair_pb, single_nb, epsilon_m, branching)
    sigma_eff_stat_f = sigma_eff_central * f_stat / f_dps
    sigma_eff_stat_pair = sigma_eff_central * pair_stat / pair_pb
    sigma_eff_stat_single = sigma_eff_central * 2.0 * single_stat / single_nb
    sigma_eff_stat_f_pair = quadrature(sigma_eff_stat_f, sigma_eff_stat_pair)
    sigma_eff_stat = quadrature(sigma_eff_stat_f_pair, sigma_eff_stat_single)

    sigma_eff_syst_f_up = (
        sigma_eff(f_dps - f_syst_down, pair_pb, single_nb, epsilon_m, branching)
        - sigma_eff_central
    )
    sigma_eff_syst_f_down = (
        sigma_eff_central
        - sigma_eff(f_dps + f_syst_up, pair_pb, single_nb, epsilon_m, branching)
    )
    sigma_eff_syst_pair_up = (
        sigma_eff(
            f_dps,
            pair_pb * (1.0 - pair_non_br_syst_rel),
            single_nb,
            epsilon_m,
            branching,
        )
        - sigma_eff_central
    )
    sigma_eff_syst_pair_down = (
        sigma_eff_central
        - sigma_eff(
            f_dps,
            pair_pb * (1.0 + pair_non_br_syst_rel),
            single_nb,
            epsilon_m,
            branching,
        )
    )
    sigma_eff_syst_f_pair_up = quadrature(
        sigma_eff_syst_f_up, sigma_eff_syst_pair_up
    )
    sigma_eff_syst_f_pair_down = quadrature(
        sigma_eff_syst_f_down, sigma_eff_syst_pair_down
    )
    sigma_eff_single_bin_up = (
        sigma_eff(f_dps, pair_pb, single_nb + single_bin_syst, epsilon_m, branching)
        - sigma_eff_central
    )
    sigma_eff_single_bin_down = (
        sigma_eff_central
        - sigma_eff(f_dps, pair_pb, single_nb - single_bin_syst, epsilon_m, branching)
    )
    sigma_eff_single_lumi_up = (
        sigma_eff(
            f_dps,
            pair_pb,
            single_nb * (1.0 + single_lumi_rel),
            epsilon_m,
            branching,
        )
        - sigma_eff_central
    )
    sigma_eff_single_lumi_down = (
        sigma_eff_central
        - sigma_eff(
            f_dps,
            pair_pb,
            single_nb * (1.0 - single_lumi_rel),
            epsilon_m,
            branching,
        )
    )
    sigma_eff_syst_up = quadrature(
        sigma_eff_syst_f_pair_up,
        sigma_eff_single_bin_up,
        sigma_eff_single_lumi_up,
    )
    sigma_eff_syst_down = quadrature(
        sigma_eff_syst_f_pair_down,
        sigma_eff_single_bin_down,
        sigma_eff_single_lumi_down,
    )

    relative_results = {
        "sigma_dps_mjj_ge_7p5": {
            "stat_percent": percent(sigma_dps_stat, sigma_dps_cut),
            "syst_up_percent": percent(sigma_dps_syst_up, sigma_dps_cut),
            "syst_down_percent": percent(sigma_dps_syst_down, sigma_dps_cut),
        },
        "sigma_dps_no_mjj_cut": {
            "stat_percent": percent(sigma_dps_no_mass_stat, sigma_dps_no_mass),
            "syst_up_percent": percent(sigma_dps_no_mass_syst_up, sigma_dps_no_mass),
            "syst_down_percent": percent(sigma_dps_no_mass_syst_down, sigma_dps_no_mass),
        },
        "sigma_eff": {
            "stat_percent": percent(sigma_eff_stat, sigma_eff_central),
            "syst_up_percent": percent(sigma_eff_syst_up, sigma_eff_central),
            "syst_down_percent": percent(sigma_eff_syst_down, sigma_eff_central),
        },
    }
    sigma_eff_stat_components_percent = {
        "f_dps": percent(sigma_eff_stat_f, sigma_eff_central),
        "pair_total": percent(sigma_eff_stat_pair, sigma_eff_central),
        "atlas_single_jpsi": percent(sigma_eff_stat_single, sigma_eff_central),
    }
    sigma_eff_syst_components_percent = {
        "f_dps_up": percent(sigma_eff_syst_f_up, sigma_eff_central),
        "f_dps_down": percent(sigma_eff_syst_f_down, sigma_eff_central),
        "pair_total_non_br_up": percent(sigma_eff_syst_pair_up, sigma_eff_central),
        "pair_total_non_br_down": percent(sigma_eff_syst_pair_down, sigma_eff_central),
        "atlas_single_bin_up": percent(sigma_eff_single_bin_up, sigma_eff_central),
        "atlas_single_bin_down": percent(sigma_eff_single_bin_down, sigma_eff_central),
        "atlas_luminosity_up": percent(sigma_eff_single_lumi_up, sigma_eff_central),
        "atlas_luminosity_down": percent(sigma_eff_single_lumi_down, sigma_eff_central),
    }

    assumptions = {
        "cov_f_dps_pair_total_stat": 0.0,
        "cov_f_dps_pair_total_syst": 0.0,
        "f_dps_systematic_sources_combination": "CR_and_DPS_MC_usage_in_quadrature",
        "pair_systematic_sources_combination": "recorded_relative_sources_in_quadrature",
        "single_jpsi_uncertainties": "ATLAS_stat_bin_systematic_and_luminosity_included",
        "branching_fraction": "fully_correlated_cancellation_in_sigma_eff_only",
        "epsilon_m_uncertainty": "pending_not_included_not_zero",
        "f_dps_input_status": f_status,
        "pair_cross_section_input_status": pair_status,
        "pair_pending_systematic_sources": sorted(pair_pending_sources),
    }
    inputs = {
        "status": "frozen_snapshot_for_current_candidate",
        "source_paths": {
            "f_dps_summary": str(f_path.resolve()),
            "mass_extrapolation": str(mass_path.resolve()),
            "pair_cross_section": str(pair_path.resolve()),
            "single_jpsi_and_branching_reference": str(reference_path.resolve()),
        },
        "source_sha256": {
            "f_dps_summary": sha256(f_path),
            "mass_extrapolation": sha256(mass_path),
            "pair_cross_section": sha256(pair_path),
            "single_jpsi_and_branching_reference": sha256(reference_path),
        },
        "f_dps": {
            "central": f_dps,
            "stat": f_stat,
            "cr_syst_up": f_cr_up,
            "cr_syst_down": f_cr_down,
            "dps_mc_usage_syst": f_usage,
            "combined_syst_up": f_syst_up,
            "combined_syst_down": f_syst_down,
        },
        "pair_cross_section": pair,
        "mass_extrapolation": {
            "all_cross_event_pairs": n_all,
            "pass_mJJ_ge_7p5": n_pass,
            "epsilon_m": epsilon_m,
            "uncertainty_status": "pending_not_included_not_zero",
        },
        "single_jpsi_cross_section": single,
        "branching_fraction": branching_input,
        "assumptions": assumptions,
    }
    results = {
        "status": result_status,
        "scope": result_scope,
        "assumptions": assumptions,
        "f_dps": inputs["f_dps"],
        "mass_extrapolation": {
            "epsilon_m": epsilon_m,
            "inverse_factor": inverse_epsilon,
            "uncertainty_included": False,
            "uncertainty_status": "pending_not_included_not_zero",
        },
        "sigma_dps_mjj_ge_7p5_pb": {
            "central": sigma_dps_cut,
            "stat": sigma_dps_stat,
            "syst_up": sigma_dps_syst_up,
            "syst_down": sigma_dps_syst_down,
            **relative_results["sigma_dps_mjj_ge_7p5"],
            "components": {
                "stat_from_f_dps": sigma_dps_stat_f,
                "stat_from_pair_total": sigma_dps_stat_pair,
                "syst_from_f_dps_up": sigma_dps_syst_f_up,
                "syst_from_f_dps_down": sigma_dps_syst_f_down,
                "syst_from_pair_total": sigma_dps_syst_pair,
            },
        },
        "sigma_dps_no_mjj_cut_pb": {
            "central": sigma_dps_no_mass,
            "stat": sigma_dps_no_mass_stat,
            "syst_up": sigma_dps_no_mass_syst_up,
            "syst_down": sigma_dps_no_mass_syst_down,
            **relative_results["sigma_dps_no_mjj_cut"],
        },
        "sigma_eff_mb": {
            "central": sigma_eff_central,
            "stat": sigma_eff_stat,
            "syst_up": sigma_eff_syst_up,
            "syst_down": sigma_eff_syst_down,
            **relative_results["sigma_eff"],
            "components": {
                "stat_from_f_dps": sigma_eff_stat_f,
                "stat_from_pair_total": sigma_eff_stat_pair,
                "stat_f_dps_pair_subtotal": sigma_eff_stat_f_pair,
                "stat_from_ATLAS_single": sigma_eff_stat_single,
                "syst_from_f_dps_up": sigma_eff_syst_f_up,
                "syst_from_f_dps_down": sigma_eff_syst_f_down,
                "syst_from_pair_total_non_BR_up": sigma_eff_syst_pair_up,
                "syst_from_pair_total_non_BR_down": sigma_eff_syst_pair_down,
                "syst_f_dps_pair_subtotal_up": sigma_eff_syst_f_pair_up,
                "syst_f_dps_pair_subtotal_down": sigma_eff_syst_f_pair_down,
                "syst_from_ATLAS_single_bin_up": sigma_eff_single_bin_up,
                "syst_from_ATLAS_single_bin_down": sigma_eff_single_bin_down,
                "syst_from_ATLAS_luminosity_up": sigma_eff_single_lumi_up,
                "syst_from_ATLAS_luminosity_down": sigma_eff_single_lumi_down,
                "branching_fraction_effect": 0.0,
                "stat_percent": sigma_eff_stat_components_percent,
                "syst_percent": sigma_eff_syst_components_percent,
            },
        },
    }

    (output_path / "inputs.json").write_text(
        json.dumps(inputs, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output_path / "results.json").write_text(
        json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    with (output_path / "results.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["quantity", "phase_space", "central", "stat", "stat_percent", "syst_up", "syst_up_percent", "syst_down", "syst_down_percent", "unit"])
        writer.writerow(["sigma_DPS", "mJJ>=7.5", sigma_dps_cut, sigma_dps_stat, relative_results["sigma_dps_mjj_ge_7p5"]["stat_percent"], sigma_dps_syst_up, relative_results["sigma_dps_mjj_ge_7p5"]["syst_up_percent"], sigma_dps_syst_down, relative_results["sigma_dps_mjj_ge_7p5"]["syst_down_percent"], "pb"])
        writer.writerow(["sigma_DPS", "no_mJJ_cut", sigma_dps_no_mass, sigma_dps_no_mass_stat, relative_results["sigma_dps_no_mjj_cut"]["stat_percent"], sigma_dps_no_mass_syst_up, relative_results["sigma_dps_no_mjj_cut"]["syst_up_percent"], sigma_dps_no_mass_syst_down, relative_results["sigma_dps_no_mjj_cut"]["syst_down_percent"], "pb"])
        writer.writerow(["sigma_eff", "single_Jpsi_matched_no_pair_mass_cut", sigma_eff_central, sigma_eff_stat, relative_results["sigma_eff"]["stat_percent"], sigma_eff_syst_up, relative_results["sigma_eff"]["syst_up_percent"], sigma_eff_syst_down, relative_results["sigma_eff"]["syst_down_percent"], "mb"])

    component_rows = [
        ["sigma_DPS_mJJ_ge_7p5", "f_DPS", "stat", sigma_dps_stat_f, sigma_dps_stat_f, "pb", "included_independent_of_pair_total"],
        ["sigma_DPS_mJJ_ge_7p5", "pair_total", "stat", sigma_dps_stat_pair, sigma_dps_stat_pair, "pb", "included_independent_of_f_DPS"],
        ["sigma_DPS_mJJ_ge_7p5", "f_DPS", "syst", sigma_dps_syst_f_up, sigma_dps_syst_f_down, "pb", "CR_and_DPS_MC_usage_combined_in_quadrature"],
        ["sigma_DPS_mJJ_ge_7p5", "pair_total", "syst", sigma_dps_syst_pair, sigma_dps_syst_pair, "pb", "included_independent_of_f_DPS_accepted_newdata_sources"],
        ["sigma_eff", "f_DPS", "stat", sigma_eff_stat_f, sigma_eff_stat_f, "mb", "included_independent_of_pair_total"],
        ["sigma_eff", "pair_total", "stat", sigma_eff_stat_pair, sigma_eff_stat_pair, "mb", "included_independent_of_f_DPS"],
        ["sigma_eff", "ATLAS_single_Jpsi", "stat", sigma_eff_stat_single, sigma_eff_stat_single, "mb", "included_independent_CMS_inputs"],
        ["sigma_eff", "f_DPS", "syst_endpoint", sigma_eff_syst_f_up, sigma_eff_syst_f_down, "mb", "included_independent_of_pair_total"],
        ["sigma_eff", "pair_total_non_BR", "syst_endpoint", sigma_eff_syst_pair_up, sigma_eff_syst_pair_down, "mb", "included_independent_of_f_DPS"],
        ["sigma_eff", "pair_branching_fraction", "syst", 0.0, 0.0, "mb", "cancelled_exactly"],
        ["sigma_eff", "ATLAS_single_bin", "syst_endpoint", sigma_eff_single_bin_up, sigma_eff_single_bin_down, "mb", "included"],
        ["sigma_eff", "ATLAS_luminosity", "syst_endpoint", sigma_eff_single_lumi_up, sigma_eff_single_lumi_down, "mb", "included_independent_CMS_luminosity"],
        ["sigma_eff", "epsilon_m", "stat_or_model", "", "", "mb", "pending_not_included_not_zero"],
    ]
    for name in sorted(pair_pending_sources):
        component_rows.extend([
            ["sigma_DPS_mJJ_ge_7p5", f"pair_total_{name}", "syst",
             "", "", "pb", "pending_not_included_not_zero"],
            ["sigma_eff", f"pair_total_{name}", "syst",
             "", "", "mb", "pending_not_included_not_zero"],
        ])
    with (output_path / "uncertainty_components.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["quantity", "source", "type", "effect_up", "effect_down", "unit", "treatment", "effect_up_percent", "effect_down_percent"])
        for row in component_rows:
            central = sigma_eff_central if row[0] == "sigma_eff" else sigma_dps_cut
            up = percent(float(row[3]), central) if row[3] != "" else ""
            down = percent(float(row[4]), central) if row[4] != "" else ""
            writer.writerow(row + [up, down])

    summary_lines = [
        f"status={result_status}",
        f"scope={result_scope}",
        f"f_dps={f_dps:.12g}",
        f"f_dps_stat={f_stat:.12g}",
        f"f_dps_combined_syst_up={f_syst_up:.12g}",
        f"f_dps_combined_syst_down={f_syst_down:.12g}",
        f"pair_total_pb={pair_pb:.12g}",
        f"pair_total_stat_pb={pair_stat:.12g}",
        f"pair_total_syst_pb={pair_syst_pb:.12g}",
        f"epsilon_m={epsilon_m:.15g}",
        "epsilon_m_uncertainty_status=pending_not_included_not_zero",
        f"sigma_dps_mjj_ge_7p5_pb={sigma_dps_cut:.12g}",
        f"sigma_dps_mjj_ge_7p5_stat_pb={sigma_dps_stat:.12g}",
        f"sigma_dps_mjj_ge_7p5_stat_percent={relative_results['sigma_dps_mjj_ge_7p5']['stat_percent']:.12g}",
        f"sigma_dps_mjj_ge_7p5_syst_up_pb={sigma_dps_syst_up:.12g}",
        f"sigma_dps_mjj_ge_7p5_syst_up_percent={relative_results['sigma_dps_mjj_ge_7p5']['syst_up_percent']:.12g}",
        f"sigma_dps_mjj_ge_7p5_syst_down_pb={sigma_dps_syst_down:.12g}",
        f"sigma_dps_mjj_ge_7p5_syst_down_percent={relative_results['sigma_dps_mjj_ge_7p5']['syst_down_percent']:.12g}",
        f"sigma_dps_no_mjj_cut_pb={sigma_dps_no_mass:.12g}",
        f"sigma_dps_no_mjj_cut_stat_pb={sigma_dps_no_mass_stat:.12g}",
        f"sigma_dps_no_mjj_cut_stat_percent={relative_results['sigma_dps_no_mjj_cut']['stat_percent']:.12g}",
        f"sigma_dps_no_mjj_cut_syst_up_pb={sigma_dps_no_mass_syst_up:.12g}",
        f"sigma_dps_no_mjj_cut_syst_up_percent={relative_results['sigma_dps_no_mjj_cut']['syst_up_percent']:.12g}",
        f"sigma_dps_no_mjj_cut_syst_down_pb={sigma_dps_no_mass_syst_down:.12g}",
        f"sigma_dps_no_mjj_cut_syst_down_percent={relative_results['sigma_dps_no_mjj_cut']['syst_down_percent']:.12g}",
        f"sigma_eff_mb={sigma_eff_central:.12g}",
        f"sigma_eff_stat_mb={sigma_eff_stat:.12g}",
        f"sigma_eff_stat_percent={relative_results['sigma_eff']['stat_percent']:.12g}",
        f"sigma_eff_syst_up_mb={sigma_eff_syst_up:.12g}",
        f"sigma_eff_syst_up_percent={relative_results['sigma_eff']['syst_up_percent']:.12g}",
        f"sigma_eff_syst_down_mb={sigma_eff_syst_down:.12g}",
        f"sigma_eff_syst_down_percent={relative_results['sigma_eff']['syst_down_percent']:.12g}",
        f"sigma_eff_stat_from_f_dps_mb={sigma_eff_stat_f:.12g}",
        f"sigma_eff_stat_from_f_dps_percent={sigma_eff_stat_components_percent['f_dps']:.12g}",
        f"sigma_eff_stat_from_pair_total_mb={sigma_eff_stat_pair:.12g}",
        f"sigma_eff_stat_from_pair_total_percent={sigma_eff_stat_components_percent['pair_total']:.12g}",
        f"sigma_eff_stat_from_atlas_single_jpsi_mb={sigma_eff_stat_single:.12g}",
        f"sigma_eff_stat_from_atlas_single_jpsi_percent={sigma_eff_stat_components_percent['atlas_single_jpsi']:.12g}",
        f"sigma_eff_syst_from_f_dps_up_percent={sigma_eff_syst_components_percent['f_dps_up']:.12g}",
        f"sigma_eff_syst_from_f_dps_down_percent={sigma_eff_syst_components_percent['f_dps_down']:.12g}",
        f"sigma_eff_syst_from_pair_total_non_br_up_percent={sigma_eff_syst_components_percent['pair_total_non_br_up']:.12g}",
        f"sigma_eff_syst_from_pair_total_non_br_down_percent={sigma_eff_syst_components_percent['pair_total_non_br_down']:.12g}",
        f"sigma_eff_syst_from_atlas_single_bin_up_percent={sigma_eff_syst_components_percent['atlas_single_bin_up']:.12g}",
        f"sigma_eff_syst_from_atlas_single_bin_down_percent={sigma_eff_syst_components_percent['atlas_single_bin_down']:.12g}",
        f"sigma_eff_syst_from_atlas_luminosity_up_percent={sigma_eff_syst_components_percent['atlas_luminosity_up']:.12g}",
        f"sigma_eff_syst_from_atlas_luminosity_down_percent={sigma_eff_syst_components_percent['atlas_luminosity_down']:.12g}",
        f"sigma_eff_stat_f_pair_subtotal_mb={sigma_eff_stat_f_pair:.12g}",
        f"sigma_eff_syst_f_pair_subtotal_up_mb={sigma_eff_syst_f_pair_up:.12g}",
        f"sigma_eff_syst_f_pair_subtotal_down_mb={sigma_eff_syst_f_pair_down:.12g}",
        "cov_f_dps_pair_total_stat=0_user_assumption",
        "cov_f_dps_pair_total_syst=0_user_assumption",
        "branching_fraction_sigma_eff=cancelled_exactly",
        f"f_dps_input_status={f_status}",
        f"pair_total_systematic_status={pair_status}",
        "pair_total_systematic_scope=" +
        ("included_numeric_sources_only_pending_not_included_not_zero"
         if pair_pending_sources else "all_accepted_sources"),
        "pair_total_pending_systematic_sources=" +
        (",".join(sorted(pair_pending_sources)) if pair_pending_sources else "none"),
        "new_roofit_fits_run=0",
    ]
    (output_path / "summary.txt").write_text(
        "\n".join(summary_lines) + "\n", encoding="utf-8"
    )
    print(
        f"sigma_DPS(mJJ>=7.5) = {sigma_dps_cut:.6f} +/- {sigma_dps_stat:.6f} "
        f"(stat.) +{sigma_dps_syst_up:.6f}/-{sigma_dps_syst_down:.6f} (syst.) pb"
    )
    print(
        f"sigma_eff = {sigma_eff_central:.6f} +/- {sigma_eff_stat:.6f} "
        f"(stat.) +{sigma_eff_syst_up:.6f}/-{sigma_eff_syst_down:.6f} (syst.) mb"
    )
    if pair_pending_sources:
        print(
            "PENDING_NOT_INCLUDED_NOT_ZERO pair systematics: "
            + ",".join(sorted(pair_pending_sources))
        )


if __name__ == "__main__":
    main()
