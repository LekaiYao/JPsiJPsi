#!/usr/bin/env python3
"""Deterministic H016 propagation from frozen JSON inputs."""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path


PACKAGE = Path(__file__).resolve().parents[1]
if len(sys.argv) > 3:
    raise SystemExit("usage: calculate_results.py [INPUT_JSON] [OUTPUT_DIR]")
INPUT_PATH = Path(sys.argv[1]).resolve() if len(sys.argv) >= 2 else PACKAGE / "inputs.json"
OUTPUT = Path(sys.argv[2]).resolve() if len(sys.argv) >= 3 else PACKAGE / "output"


def rho_error(a: float, b: float, extra: float = 0.0, rho: float = 0.0) -> float:
    variance = a * a + b * b + extra * extra + 2.0 * rho * a * b
    return math.sqrt(max(0.0, variance))


def quadrature(values) -> float:
    return math.sqrt(sum(value * value for value in values))


def sigma_eff_value(f: float, pair_pb: float, x_nb: float, epsilon_m: float, branching: float) -> float:
    sigma_dps_no_mass = f * pair_pb / epsilon_m
    visible_pair_pb = branching * branching * sigma_dps_no_mass
    return 0.001 * 0.5 * x_nb * x_nb / visible_pair_pb


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    inputs = json.loads(INPUT_PATH.read_text(encoding="utf-8"))
    frac = inputs["data_driven_fraction"]
    pair = inputs["pair_cross_section"]
    single = inputs["single_jpsi_cross_section"]
    br = inputs["branching_fraction"]
    mass = inputs["mass_extrapolation"]

    f = float(frac["central"])
    f_stat = float(frac["stat"])
    f_up = float(frac["method_syst_up"])
    f_down = float(frac["method_syst_down"])
    pair_pb = float(pair["central_pb"])
    pair_stat = float(pair["stat_pb"])
    pair_sources = pair["systematic_sources_relative"]
    pair_source_rel = {name: float(source["relative"]) for name, source in pair_sources.items()}
    pair_reconstructed_rel = quadrature(pair_source_rel.values())
    pair_reconstructed_syst_pb = pair_pb * pair_reconstructed_rel

    epsilon_m = float(mass["pass_mJJ_gt_7p5"]) / float(mass["all_cross_event_pairs"])
    extrapolation = 1.0 / epsilon_m

    sigma_dps_fid = f * pair_pb
    stat_from_f = pair_pb * f_stat
    stat_from_pair = f * pair_stat
    sigma_dps_stat_scan = {
        str(rho): rho_error(stat_from_f, stat_from_pair, rho=rho)
        for rho in (-1.0, 0.0, 1.0)
    }
    f_groups = frac["method_groups"]
    sigma_dps_f_group_effects = {
        name: {"up": pair_pb * float(group["up"]), "down": pair_pb * float(group["down"])}
        for name, group in f_groups.items()
    }
    sigma_dps_pair_source_effects = {
        name: sigma_dps_fid * relative for name, relative in pair_source_rel.items()
    }
    sigma_dps_syst_up = quadrature(
        [effect["up"] for effect in sigma_dps_f_group_effects.values()]
        + list(sigma_dps_pair_source_effects.values())
    )
    sigma_dps_syst_down = quadrature(
        [effect["down"] for effect in sigma_dps_f_group_effects.values()]
        + list(sigma_dps_pair_source_effects.values())
    )

    sigma_dps_no_mass = sigma_dps_fid * extrapolation
    sigma_dps_no_mass_stat_scan = {
        rho: value * extrapolation for rho, value in sigma_dps_stat_scan.items()
    }
    sigma_dps_no_mass_syst_up = sigma_dps_syst_up * extrapolation
    sigma_dps_no_mass_syst_down = sigma_dps_syst_down * extrapolation

    x_nb = float(single["visible_central_nb"])
    x_stat = float(single["visible_stat_nb"])
    x_bin_syst = float(single["visible_bin_syst_nb"])
    atlas_lumi_rel = float(single["luminosity_relative"])
    branching = float(br["central"])
    visible_pair_pb = branching * branching * sigma_dps_no_mass
    sigma_eff_mb = sigma_eff_value(f, pair_pb, x_nb, epsilon_m, branching)

    sigma_eff_stat_from_f = sigma_eff_mb * f_stat / f
    sigma_eff_stat_from_pair = sigma_eff_mb * pair_stat / pair_pb
    sigma_eff_stat_from_single = sigma_eff_mb * 2.0 * x_stat / x_nb
    sigma_eff_stat_scan = {
        str(rho): rho_error(
            sigma_eff_stat_from_f,
            sigma_eff_stat_from_pair,
            extra=sigma_eff_stat_from_single,
            rho=rho,
        )
        for rho in (-1.0, 0.0, 1.0)
    }

    sigma_eff_f_group_effects = {}
    for name, group in f_groups.items():
        group_up = float(group["up"])
        group_down = float(group["down"])
        effect_up = (
            sigma_eff_value(f - group_down, pair_pb, x_nb, epsilon_m, branching) - sigma_eff_mb
            if group_down else 0.0
        )
        effect_down = (
            sigma_eff_mb - sigma_eff_value(f + group_up, pair_pb, x_nb, epsilon_m, branching)
            if group_up else 0.0
        )
        sigma_eff_f_group_effects[name] = {"up": effect_up, "down": effect_down}

    sigma_eff_pair_source_effects = {}
    for name, relative in pair_source_rel.items():
        if name == "branching_fraction":
            sigma_eff_pair_source_effects[name] = {"up": 0.0, "down": 0.0}
            continue
        effect_up = sigma_eff_value(
            f, pair_pb * (1.0 - relative), x_nb, epsilon_m, branching
        ) - sigma_eff_mb
        effect_down = sigma_eff_mb - sigma_eff_value(
            f, pair_pb * (1.0 + relative), x_nb, epsilon_m, branching
        )
        sigma_eff_pair_source_effects[name] = {"up": effect_up, "down": effect_down}

    sigma_eff_single_bin_effect = {
        "up": sigma_eff_value(f, pair_pb, x_nb + x_bin_syst, epsilon_m, branching) - sigma_eff_mb,
        "down": sigma_eff_mb - sigma_eff_value(f, pair_pb, x_nb - x_bin_syst, epsilon_m, branching),
    }
    sigma_eff_single_lumi_effect = {
        "up": sigma_eff_value(f, pair_pb, x_nb * (1.0 + atlas_lumi_rel), epsilon_m, branching) - sigma_eff_mb,
        "down": sigma_eff_mb - sigma_eff_value(f, pair_pb, x_nb * (1.0 - atlas_lumi_rel), epsilon_m, branching),
    }
    sigma_eff_syst_up = quadrature(
        [effect["up"] for effect in sigma_eff_f_group_effects.values()]
        + [effect["up"] for effect in sigma_eff_pair_source_effects.values()]
        + [sigma_eff_single_bin_effect["up"], sigma_eff_single_lumi_effect["up"]]
    )
    sigma_eff_syst_down = quadrature(
        [effect["down"] for effect in sigma_eff_f_group_effects.values()]
        + [effect["down"] for effect in sigma_eff_pair_source_effects.values()]
        + [sigma_eff_single_bin_effect["down"], sigma_eff_single_lumi_effect["down"]]
    )

    frozen_h016_reference = (
        frac.get("status") == "user_frozen_and_H014_consumer_accepted"
    )
    result_status = (
        "producer_ready_pending_consumer_review"
        if frozen_h016_reference
        else "candidate_complete_pending_user_confirmation"
    )
    fraction_status = (
        "accepted" if frozen_h016_reference else "candidate_pending_user_confirmation"
    )
    first_difference = (
        "No disagreement for the accepted f_DPS or the preliminary sigma_DPS(mJJ>7.5) central value."
        if frozen_h016_reference
        else "All values are recalculated from the selected inputs.json; no historical display number is hard-coded."
    )

    results = {
        "status": result_status,
        "differences_from_h016_request_baseline": [
            first_difference,
            "The current AN and earlier H013-style sigma_eff calculation omitted the pair-specific mJJ extrapolation.",
            "H016 applies epsilon_m from factorized DPS GEN-only event mixing before calculating sigma_eff.",
            "No replica-matched Cov(f_DPS,sigma_pair_total) exists; rho=0 remains a preliminary reference, not a final statistical conclusion."
        ],
        "mass_extrapolation": {
            "epsilon_m": epsilon_m,
            "inverse_factor": extrapolation,
            "uncertainty_included": False,
        },
        "f_DPS": {
            "central": f,
            "stat": f_stat,
            "method_syst_up": f_up,
            "method_syst_down": f_down,
            "status": fraction_status,
            "display": (
                f"f_DPS = {f:.3f} +/- {f_stat:.3f} (stat.) "
                f"+{f_up:.3f}/-{f_down:.3f} (method syst.)"
            ),
        },
        "sigma_DPS_mJJ_gt_7p5_pb": {
            "central": sigma_dps_fid,
            "stat_rho_scan": sigma_dps_stat_scan,
            "syst_up": sigma_dps_syst_up,
            "syst_down": sigma_dps_syst_down,
            "components": {
                "stat_from_f_DPS": stat_from_f,
                "stat_from_pair_total": stat_from_pair,
                "f_DPS_method_groups": sigma_dps_f_group_effects,
                "AN_pair_sources": sigma_dps_pair_source_effects,
                "AN_pair_reconstructed_relative": pair_reconstructed_rel,
                "AN_pair_reconstructed_syst_pb": pair_reconstructed_syst_pb,
                "AN_pair_reported_rounded_syst_pb": float(pair["reported_syst_pb_rounded"]),
            },
            "status": "preliminary; rho=0 statistical reference; AN percentages used directly; unknown shared-systematic correlations use conditional independence",
            "display": (
                f"sigma_DPS(mJJ>7.5) = {sigma_dps_fid:.2f} +/- "
                f"{sigma_dps_stat_scan['0.0']:.2f} (stat., rho=0) "
                f"+{sigma_dps_syst_up:.2f}/-{sigma_dps_syst_down:.2f} (syst.) pb"
            ),
        },
        "sigma_DPS_no_mJJ_cut_pb": {
            "central": sigma_dps_no_mass,
            "stat_rho_scan": sigma_dps_no_mass_stat_scan,
            "syst_up": sigma_dps_no_mass_syst_up,
            "syst_down": sigma_dps_no_mass_syst_down,
            "status": "derived preliminary; epsilon_m fixed without uncertainty",
            "display": (
                f"sigma_DPS(no mJJ cut) = {sigma_dps_no_mass:.2f} +/- "
                f"{sigma_dps_no_mass_stat_scan['0.0']:.2f} (stat., rho=0) "
                f"+{sigma_dps_no_mass_syst_up:.2f}/-{sigma_dps_no_mass_syst_down:.2f} (syst.) pb"
            ),
        },
        "sigma_eff_mb": {
            "central": sigma_eff_mb,
            "stat_rho_scan": sigma_eff_stat_scan,
            "syst_up": sigma_eff_syst_up,
            "syst_down": sigma_eff_syst_down,
            "components": {
                "stat_from_f_DPS": sigma_eff_stat_from_f,
                "stat_from_pair_total": sigma_eff_stat_from_pair,
                "stat_from_ATLAS_single": sigma_eff_stat_from_single,
                "f_DPS_method_endpoint_sources": sigma_eff_f_group_effects,
                "CMS_pair_endpoint_sources": sigma_eff_pair_source_effects,
                "ATLAS_single_bin_endpoint": sigma_eff_single_bin_effect,
                "ATLAS_luminosity_endpoint": sigma_eff_single_lumi_effect,
                "AN_pair_reconstructed_relative": pair_reconstructed_rel,
                "AN_pair_reconstructed_syst_pb": pair_reconstructed_syst_pb,
                "branching_fraction_effect": 0.0,
                "mass_extrapolation_uncertainty_effect": 0.0,
            },
            "visible_pair_cross_section_pb": visible_pair_pb,
            "status": "preliminary; rho=0 statistical reference; exact source-by-source systematic endpoints; epsilon_m fixed by user gate; BR cancelled",
            "display": (
                f"sigma_eff = {sigma_eff_mb:.2f} +/- "
                f"{sigma_eff_stat_scan['0.0']:.2f} (stat., rho=0) "
                f"+{sigma_eff_syst_up:.2f}/-{sigma_eff_syst_down:.2f} (syst.) mb"
            ),
        },
    }

    (OUTPUT / "results.json").write_text(
        json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    rows = [
        ["f_DPS", "mJJ>7.5 fiducial pair fraction", f, f_stat, f_up, f_down, "", fraction_status],
        ["sigma_DPS", "mJJ>7.5 GeV", sigma_dps_fid, sigma_dps_stat_scan["0.0"], sigma_dps_syst_up, sigma_dps_syst_down, "pb", "preliminary_rho0"],
        ["sigma_DPS", "no mJJ cut", sigma_dps_no_mass, sigma_dps_no_mass_stat_scan["0.0"], sigma_dps_no_mass_syst_up, sigma_dps_no_mass_syst_down, "pb", "derived_preliminary_epsilon_fixed"],
        ["sigma_eff", "single-Jpsi-matched; no pair-mass cut", sigma_eff_mb, sigma_eff_stat_scan["0.0"], sigma_eff_syst_up, sigma_eff_syst_down, "mb", "preliminary_rho0_epsilon_fixed"],
    ]
    with (OUTPUT / "results.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["quantity", "phase_space", "central", "stat_rho0", "syst_up", "syst_down", "unit", "status"])
        writer.writerows(rows)

    source_rows = [
        ["f_DPS", "fixed-500 full source-event bootstrap", "stat", f_stat, f_stat, "absolute fraction", "sample SD; covers the full f_DPS workflow", "yes", "accepted" if frozen_h016_reference else "candidate pending user confirmation"],
        ["sigma_DPS", "f_DPS bootstrap", "stat", stat_from_f, stat_from_f, "pb", "Cov with pair total unknown; rho scan supplied", "yes", "preliminary propagation"],
        ["sigma_DPS", "pair total statistical", "stat", stat_from_pair, stat_from_pair, "pb", "Cov with f_DPS unknown; rho scan supplied", "yes", "preliminary propagation"],
        ["sigma_eff", "f_DPS bootstrap", "stat", sigma_eff_stat_from_f, sigma_eff_stat_from_f, "mb", "Cov with pair total unknown; rho scan supplied", "yes", "preliminary propagation"],
        ["sigma_eff", "pair total statistical", "stat", sigma_eff_stat_from_pair, sigma_eff_stat_from_pair, "mb", "Cov with f_DPS unknown; rho scan supplied", "yes", "preliminary propagation"],
        ["sigma_eff", "ATLAS single-Jpsi statistical", "stat", sigma_eff_stat_from_single, sigma_eff_stat_from_single, "mb", "independent of CMS inputs", "yes", "preliminary propagation"],
    ]

    f_labels = {
        "single_j_category": "single-J category envelope",
        "cr_definition_and_purity": "CR definition/purity envelope",
        "pure_dps_full_workflow_closure": "pure-DPS full-workflow closure",
    }
    for name, label in f_labels.items():
        group = f_groups[name]
        source_rows.append(["f_DPS", label, "method syst", group["up"], group["down"], "absolute fraction", "included group; separate up/down quadrature", "yes", "accepted conditional method source" if frozen_h016_reference else "candidate pending user confirmation"])
        dps_effect = sigma_dps_f_group_effects[name]
        source_rows.append(["sigma_DPS", label, "method syst", dps_effect["up"], dps_effect["down"], "pb", "conditional independence from AN pair sources; matched nuisance variations unavailable", "yes", "preliminary propagation"])
        eff_effect = sigma_eff_f_group_effects[name]
        source_rows.append(["sigma_eff", label, "method syst endpoint", eff_effect["up"], eff_effect["down"], "mb", "exact inverse-denominator endpoints; conditional independence from AN pair sources", "yes", "preliminary propagation"])

    pair_correlations = {
        "branching_fraction": "does not enter f_DPS; retained for sigma_DPS; fully correlated cancellation in sigma_eff",
        "cms_luminosity": "does not enter f_DPS; independent of ATLAS luminosity",
        "acceptance_efficiency": "can affect both pair total and corrected f_DPS workflow; correlation unknown; conditional independence only",
        "lifetime_variable": "can overlap the PP fit chain used by pair total and f_DPS cells; correlation unknown; conditional independence only",
        "fitter_stability": "can overlap the PP fit chain used by pair total and f_DPS cells; correlation unknown; conditional independence only",
    }
    for name, source in pair_sources.items():
        label = source["label"]
        dps_effect = sigma_dps_pair_source_effects[name]
        source_rows.append(["sigma_DPS", label, "AN pair syst", dps_effect, dps_effect, "pb", pair_correlations[name] + "; AN percentage used directly without rescaling", "yes", "preliminary propagation"])
        eff_effect = sigma_eff_pair_source_effects[name]
        included = "cancelled" if name == "branching_fraction" else "yes"
        status = "cancelled" if name == "branching_fraction" else "preliminary propagation"
        source_rows.append(["sigma_eff", label, "AN pair syst endpoint", eff_effect["up"], eff_effect["down"], "mb", pair_correlations[name] + "; exact endpoint", included, status])

    source_rows.extend([
        ["sigma_DPS", "AN reported aggregate 2.5 pb", "rounded display", "", "", "pb", f"five percentages reconstruct {pair_reconstructed_syst_pb:.12g} pb; 2.5 pb is not used as an exact propagation input", "no separate term", "display-only comparison"],
        ["sigma_eff", "ATLAS integrated-bin systematic", "syst endpoint", sigma_eff_single_bin_effect["up"], sigma_eff_single_bin_effect["down"], "mb", "HEPData integrated uncertainty; exact numerator-squared endpoint", "yes", "preliminary propagation"],
        ["sigma_eff", "ATLAS luminosity", "syst endpoint", sigma_eff_single_lumi_effect["up"], sigma_eff_single_lumi_effect["down"], "mb", "1.13%; independent of CMS luminosity; exact numerator-squared endpoint", "yes", "preliminary propagation"],
        ["sigma_eff", "mJJ extrapolation efficiency", "model/stat", 0.0, 0.0, "mb", "fixed central by explicit 2026-08-29 user gate; future source-event jackknife and model-correlation proposal recorded", "not in current round", "preliminary gate; not zero uncertainty"],
        ["sigma_DPS and sigma_eff", "Cov(f_DPS,pair total)", "stat covariance", "", "", "", "rho=-1,0,+1 sensitivity; no replica-matched estimator; rho=0 preliminary use approved by user", "not resolved", "future covariance improvement"],
    ])
    with (OUTPUT / "uncertainty_source_map.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["quantity", "source", "type", "effect_up", "effect_down", "unit", "correlation_treatment", "included", "status"])
        writer.writerows(source_rows)

    print(results["f_DPS"]["display"])
    print(results["sigma_DPS_mJJ_gt_7p5_pb"]["display"])
    print(results["sigma_DPS_no_mJJ_cut_pb"]["display"])
    print(results["sigma_eff_mb"]["display"])


if __name__ == "__main__":
    main()
