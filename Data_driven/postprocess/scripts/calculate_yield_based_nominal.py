#!/usr/bin/env python3
"""Calculate the yield-based f_DPS, sigma_DPS, and sigma_eff candidate.

The statistical definitions are intentionally those requested on 2026-09-03:

* f_DPS: CR DPS-yield statistical error and total-cross-section statistical
  error, treated as independent;
* sigma_DPS: CR DPS-yield statistical error only;
* half-vs-full self-mixing is retained as a cross-check, not as a formal
  statistical source.

The Data fitter-stability systematic is deliberately left pending.  Any
fitter value already present in the transferred cross-section CSV is recorded
for provenance but is not propagated.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import subprocess
import sys
from datetime import datetime
from pathlib import Path


EXPECTED_ROOT_VERSION = "6.40.02"
LUMINOSITY_FB = 36.684


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


def inverse_up(central: float, relative_down: float) -> float:
    """Positive effect on an inverse quantity from a downward denominator."""
    return central / (1.0 - relative_down) - central


def inverse_down(central: float, relative_up: float) -> float:
    """Negative-effect magnitude on an inverse quantity from an upward denominator."""
    return central - central / (1.0 + relative_up)


def require_close(name: str, actual: float, expected: float) -> None:
    if not math.isclose(actual, expected, rel_tol=5e-11, abs_tol=1e-12):
        raise RuntimeError(f"{name} mismatch: {actual} vs {expected}")


def main() -> None:
    if len(sys.argv) != 12:
        raise SystemExit(
            "usage: calculate_yield_based_nominal.py PROMOTION_SUMMARY "
            "HALF_NOMINAL HALF_PHI70 HALF_PHI110 HALF_DY2P0 FULL_DIRECT "
            "DATA_CROSS_SECTIONS DATA_METADATA MASS_SUMMARY REFERENCE_JSON "
            "OUTPUT_DIR"
        )

    (
        promotion_path,
        nominal_path,
        phi70_path,
        phi110_path,
        dy2p0_path,
        full_direct_path,
        cross_sections_path,
        data_metadata_path,
        mass_path,
        reference_path,
        output_path,
    ) = map(Path, sys.argv[1:])

    input_paths = {
        "confirmed_method_promotion": promotion_path,
        "half_selfmix_nominal": nominal_path,
        "half_selfmix_phi70": phi70_path,
        "half_selfmix_phi110": phi110_path,
        "half_selfmix_dy2p0": dy2p0_path,
        "full_direct_dps_reference": full_direct_path,
        "data_cross_sections": cross_sections_path,
        "data_cross_section_metadata": data_metadata_path,
        "mass_extrapolation": mass_path,
        "single_jpsi_and_branching_reference": reference_path,
        "postprocess_script": Path(__file__).resolve(),
    }
    for label, path in input_paths.items():
        if not path.is_file():
            raise RuntimeError(f"missing {label}: {path}")
    if output_path.exists():
        raise RuntimeError(f"refusing to overwrite existing output: {output_path}")

    root_version = subprocess.check_output(
        ["root-config", "--version"], text=True
    ).strip()
    if root_version != EXPECTED_ROOT_VERSION:
        raise RuntimeError(
            f"ROOT version mismatch: {root_version} != {EXPECTED_ROOT_VERSION}"
        )
    repo_root = Path(__file__).resolve().parents[3]
    git_head = subprocess.check_output(
        ["git", "-C", str(repo_root), "rev-parse", "HEAD"], text=True
    ).strip()

    promotion = read_kv(promotion_path)
    if promotion.get("status") != "nominal_complete_user_confirmed":
        raise RuntimeError("method promotion is not confirmed")
    nominal = read_kv(nominal_path)
    variations = {
        "phi70": read_kv(phi70_path),
        "phi110": read_kv(phi110_path),
        "dy2p0": read_kv(dy2p0_path),
    }
    full_direct = read_kv(full_direct_path)
    mass = read_kv(mass_path)
    if mass.get("status") != "complete":
        raise RuntimeError("mass extrapolation input is incomplete")
    reference = json.loads(reference_path.read_text(encoding="utf-8"))
    data_metadata = read_kv(data_metadata_path)
    if data_metadata.get("status") != "complete":
        raise RuntimeError("Data cross-section campaign is incomplete")
    if data_metadata.get("promotion_status") != "accepted_by_user":
        raise RuntimeError("Data cross-section campaign is not user accepted")

    with cross_sections_path.open(newline="", encoding="utf-8") as handle:
        total_rows = [
            row for row in csv.DictReader(handle) if row["scope"] == "total"
        ]
    if len(total_rows) != 1:
        raise RuntimeError("cross_sections.csv does not contain one total row")
    total = total_rows[0]
    if int(total["status"]) != 0 or int(total["covQual"]) != 3:
        raise RuntimeError("accepted total fit fails status/covQual gate")
    if float(total["edm"]) >= 0.01:
        raise RuntimeError("accepted total fit fails EDM gate")

    data_control = float(nominal["data_control"])
    data_control_error = float(nominal["data_control_error"])
    dps_yield = float(nominal["dps_total"])
    if min(data_control, data_control_error, dps_yield) <= 0.0:
        raise RuntimeError("invalid nominal control/full DPS yield")
    variation_dps_yields = {
        name: float(values["dps_total"]) for name, values in variations.items()
    }
    full_direct_dps_yield = float(full_direct["dps_total"])

    total_yield = float(total["yield"])
    total_yield_error = float(total["stat_yield"])
    pair_pb = float(total["sigma_pb_per_unit"])
    pair_stat_pb = float(total["stat_pb_per_unit"])
    conversion_pb_per_yield = pair_pb / total_yield

    single = reference["single_jpsi_cross_section"]
    branching_input = reference["branching_fraction"]
    branching = float(branching_input["central"])
    direct_conversion = 1.0 / (
        LUMINOSITY_FB * 1000.0 * branching * branching
    )
    require_close(
        "cross-section conversion",
        conversion_pb_per_yield,
        direct_conversion,
    )
    require_close(
        "total statistical relative error",
        pair_stat_pb / pair_pb,
        total_yield_error / total_yield,
    )

    pair_sources = {
        "acceptance_efficiency": float(total["correction_relative"]),
        "branching_fraction": float(total["br_relative"]),
        "cms_luminosity": float(total["luminosity_relative"]),
        "lifetime_variable": float(total["lifetime_relative"]),
    }
    transferred_fitter_relative = float(total["fitter_relative"])
    reconstructed_transferred_total = quadrature(
        *pair_sources.values(), transferred_fitter_relative
    )
    require_close(
        "transferred total systematic",
        reconstructed_transferred_total,
        float(total["total_systematic_relative"]),
    )
    pair_known_relative = quadrature(*pair_sources.values())
    pair_non_br_relative = quadrature(
        *(
            value
            for name, value in pair_sources.items()
            if name != "branching_fraction"
        )
    )

    # Central values: use the same accepted total-fit conversion for the full
    # fiducial DPS yield and for the total pair yield.
    sigma_dps_cut = dps_yield * conversion_pb_per_yield
    f_dps = sigma_dps_cut / pair_pb
    require_close("f_DPS central", f_dps, dps_yield / total_yield)

    # Statistical uncertainties requested by the user.
    control_stat_relative = data_control_error / data_control
    total_stat_relative = pair_stat_pb / pair_pb
    dps_yield_stat = dps_yield * control_stat_relative
    f_stat_control = f_dps * control_stat_relative
    f_stat_total = f_dps * total_stat_relative
    f_stat = quadrature(f_stat_control, f_stat_total)
    sigma_dps_cut_stat = sigma_dps_cut * control_stat_relative

    # Method systematics are recomputed at fixed total cross section from the
    # full-fiducial DPS yields.  No total-cross-section systematic enters f_DPS.
    cr_signed_yield_shifts = {
        name: value - dps_yield for name, value in variation_dps_yields.items()
    }
    cr_up_yield = max(0.0, *cr_signed_yield_shifts.values())
    cr_down_yield = max(
        0.0, *(-value for value in cr_signed_yield_shifts.values())
    )
    dps_method_yield = abs(full_direct_dps_yield - dps_yield)
    cr_up_relative = cr_up_yield / dps_yield
    cr_down_relative = cr_down_yield / dps_yield
    dps_method_relative = dps_method_yield / dps_yield
    method_up_relative = quadrature(cr_up_relative, dps_method_relative)
    method_down_relative = quadrature(cr_down_relative, dps_method_relative)
    f_cr_up = f_dps * cr_up_relative
    f_cr_down = f_dps * cr_down_relative
    f_dps_method = f_dps * dps_method_relative
    f_syst_up = quadrature(f_cr_up, f_dps_method)
    f_syst_down = quadrature(f_cr_down, f_dps_method)

    # sigma_DPS gets the two method sources and all currently known pair-total
    # sources.  The new fitter-stability value remains explicitly pending.
    sigma_dps_syst_up_relative = quadrature(
        method_up_relative, pair_known_relative
    )
    sigma_dps_syst_down_relative = quadrature(
        method_down_relative, pair_known_relative
    )
    sigma_dps_cut_syst_up = sigma_dps_cut * sigma_dps_syst_up_relative
    sigma_dps_cut_syst_down = sigma_dps_cut * sigma_dps_syst_down_relative

    n_all = int(mass["all_cross_event_pairs"])
    n_pass = int(mass["pass_mJJ_gt_7p5"])
    epsilon_m = n_pass / n_all
    if not 0.0 < epsilon_m <= 1.0:
        raise RuntimeError("invalid epsilon_m")
    inverse_epsilon = 1.0 / epsilon_m
    sigma_dps_no_mass = sigma_dps_cut * inverse_epsilon
    sigma_dps_no_mass_stat = sigma_dps_cut_stat * inverse_epsilon
    sigma_dps_no_mass_syst_up = sigma_dps_cut_syst_up * inverse_epsilon
    sigma_dps_no_mass_syst_down = sigma_dps_cut_syst_down * inverse_epsilon

    single_nb = float(single["visible_central_nb"])
    single_stat_nb = float(single["visible_stat_nb"])
    single_bin_syst_nb = float(single["visible_bin_syst_nb"])
    single_lumi_relative = float(single["luminosity_relative"])
    sigma_eff_mb = (
        0.001
        * 0.5
        * single_nb
        * single_nb
        / (branching * branching * sigma_dps_no_mass)
    )

    # Original sigma_eff propagation: linear statistical components and exact
    # nonlinear endpoints for systematic sources.  Pair BR cancels exactly.
    sigma_eff_stat_from_dps = (
        sigma_eff_mb * sigma_dps_cut_stat / sigma_dps_cut
    )
    sigma_eff_stat_from_single = (
        sigma_eff_mb * 2.0 * single_stat_nb / single_nb
    )
    sigma_eff_stat = quadrature(
        sigma_eff_stat_from_dps, sigma_eff_stat_from_single
    )
    sigma_eff_method_up = inverse_up(sigma_eff_mb, method_down_relative)
    sigma_eff_method_down = inverse_down(sigma_eff_mb, method_up_relative)
    sigma_eff_pair_up = inverse_up(sigma_eff_mb, pair_non_br_relative)
    sigma_eff_pair_down = inverse_down(sigma_eff_mb, pair_non_br_relative)
    sigma_eff_single_bin_up = (
        sigma_eff_mb * ((single_nb + single_bin_syst_nb) / single_nb) ** 2
        - sigma_eff_mb
    )
    sigma_eff_single_bin_down = (
        sigma_eff_mb
        - sigma_eff_mb * ((single_nb - single_bin_syst_nb) / single_nb) ** 2
    )
    sigma_eff_single_lumi_up = (
        sigma_eff_mb * (1.0 + single_lumi_relative) ** 2 - sigma_eff_mb
    )
    sigma_eff_single_lumi_down = (
        sigma_eff_mb - sigma_eff_mb * (1.0 - single_lumi_relative) ** 2
    )
    sigma_eff_syst_up = quadrature(
        sigma_eff_method_up,
        sigma_eff_pair_up,
        sigma_eff_single_bin_up,
        sigma_eff_single_lumi_up,
    )
    sigma_eff_syst_down = quadrature(
        sigma_eff_method_down,
        sigma_eff_pair_down,
        sigma_eff_single_bin_down,
        sigma_eff_single_lumi_down,
    )

    assumptions = {
        "f_dps_statistical_sources": (
            "CR_DPS_yield_and_total_cross_section_independent_quadrature"
        ),
        "sigma_dps_statistical_sources": "CR_DPS_yield_only",
        "half_vs_full_selfmix_shift": (
            "cross_check_not_formal_statistical_uncertainty"
        ),
        "f_dps_systematic_sources": (
            "CR_boundary_and_full_direct_DPS_difference_only"
        ),
        "f_dps_cross_section_systematics": "excluded_by_user_definition",
        "sigma_dps_systematic_sources": (
            "two_method_sources_and_known_sigma_total_sources_in_quadrature"
        ),
        "fitter_stability": "pending_not_included_not_zero",
        "epsilon_m_uncertainty": "pending_not_included_not_zero",
        "branching_fraction_sigma_eff": "cancelled_exactly",
    }
    inputs = {
        "status": "frozen_snapshot_for_yield_based_candidate",
        "source_paths": {
            name: str(path.resolve()) for name, path in input_paths.items()
        },
        "source_sha256": {
            name: sha256(path) for name, path in input_paths.items()
        },
        "control_region": {
            "definition": "abs(delta_y)>=1.8_and_abs(delta_phi)<pi/2",
            "dps_yield": data_control,
            "dps_yield_error": data_control_error,
            "relative_error": control_stat_relative,
        },
        "full_fiducial_dps_yields": {
            "nominal": dps_yield,
            **variation_dps_yields,
            "full_direct_dps_reference": full_direct_dps_yield,
        },
        "pair_total": {
            "yield": total_yield,
            "yield_error": total_yield_error,
            "central_pb": pair_pb,
            "stat_pb": pair_stat_pb,
            "known_systematic_sources_relative": pair_sources,
            "known_systematic_relative": pair_known_relative,
            "transferred_fitter_relative_not_used": transferred_fitter_relative,
            "fitter_stability": "pending_not_included_not_zero",
        },
        "mass_extrapolation": {
            "epsilon_m": epsilon_m,
            "uncertainty": None,
            "uncertainty_status": "pending_not_included_not_zero",
        },
        "single_jpsi_cross_section": single,
        "branching_fraction": branching_input,
        "assumptions": assumptions,
    }
    results = {
        "status": (
            "candidate_complete_included_sources_pending_fitter_stability_"
            "and_epsilon_m_uncertainty"
        ),
        "assumptions": assumptions,
        "f_dps": {
            "status": "nominal_updated_user_defined_yield_based_scheme",
            "central": f_dps,
            "stat": f_stat,
            "stat_percent": percent(f_stat, f_dps),
            "stat_components": {
                "cr_dps_yield": f_stat_control,
                "cr_dps_yield_percent": percent(f_stat_control, f_dps),
                "sigma_total": f_stat_total,
                "sigma_total_percent": percent(f_stat_total, f_dps),
            },
            "cr_syst_up": f_cr_up,
            "cr_syst_down": f_cr_down,
            "cr_syst_up_percent": percent(f_cr_up, f_dps),
            "cr_syst_down_percent": percent(f_cr_down, f_dps),
            "dps_method_syst": f_dps_method,
            "dps_method_syst_percent": percent(f_dps_method, f_dps),
            "combined_syst_up": f_syst_up,
            "combined_syst_down": f_syst_down,
            "combined_syst_up_percent": percent(f_syst_up, f_dps),
            "combined_syst_down_percent": percent(f_syst_down, f_dps),
            "sigma_total_systematics": "not_included_by_user_definition",
        },
        "sigma_dps_mjj_ge_7p5_pb": {
            "central": sigma_dps_cut,
            "stat": sigma_dps_cut_stat,
            "stat_percent": percent(sigma_dps_cut_stat, sigma_dps_cut),
            "included_syst_up": sigma_dps_cut_syst_up,
            "included_syst_down": sigma_dps_cut_syst_down,
            "included_syst_up_percent": percent(
                sigma_dps_cut_syst_up, sigma_dps_cut
            ),
            "included_syst_down_percent": percent(
                sigma_dps_cut_syst_down, sigma_dps_cut
            ),
            "total_syst": None,
            "total_syst_status": "pending_fitter_stability_not_zero",
        },
        "sigma_dps_no_mjj_cut_pb": {
            "central": sigma_dps_no_mass,
            "stat": sigma_dps_no_mass_stat,
            "stat_percent": percent(sigma_dps_no_mass_stat, sigma_dps_no_mass),
            "included_syst_up": sigma_dps_no_mass_syst_up,
            "included_syst_down": sigma_dps_no_mass_syst_down,
            "included_syst_up_percent": percent(
                sigma_dps_no_mass_syst_up, sigma_dps_no_mass
            ),
            "included_syst_down_percent": percent(
                sigma_dps_no_mass_syst_down, sigma_dps_no_mass
            ),
            "total_syst": None,
            "total_syst_status": (
                "pending_fitter_stability_and_epsilon_m_uncertainty_not_zero"
            ),
        },
        "sigma_eff_mb": {
            "central": sigma_eff_mb,
            "stat": sigma_eff_stat,
            "stat_percent": percent(sigma_eff_stat, sigma_eff_mb),
            "included_syst_up": sigma_eff_syst_up,
            "included_syst_down": sigma_eff_syst_down,
            "included_syst_up_percent": percent(sigma_eff_syst_up, sigma_eff_mb),
            "included_syst_down_percent": percent(
                sigma_eff_syst_down, sigma_eff_mb
            ),
            "total_syst": None,
            "total_syst_status": (
                "pending_fitter_stability_and_epsilon_m_uncertainty_not_zero"
            ),
        },
    }

    component_rows = [
        ["f_DPS", "CR_DPS_yield", "stat", f_stat_control, f_stat_control, "fraction", "included_independent_of_sigma_total"],
        ["f_DPS", "sigma_total", "stat", f_stat_total, f_stat_total, "fraction", "included_independent_of_CR_DPS_yield"],
        ["f_DPS", "CR_boundary", "syst", f_cr_up, f_cr_down, "fraction", "included_no_sigma_total_systematics"],
        ["f_DPS", "DPS_template", "syst", f_dps_method, f_dps_method, "fraction", "included_no_sigma_total_systematics"],
        ["sigma_DPS_mJJ_ge_7p5", "CR_DPS_yield", "stat", sigma_dps_cut_stat, sigma_dps_cut_stat, "pb", "only_statistical_source"],
        ["sigma_DPS_mJJ_ge_7p5", "CR_boundary", "syst", sigma_dps_cut * cr_up_relative, sigma_dps_cut * cr_down_relative, "pb", "included"],
        ["sigma_DPS_mJJ_ge_7p5", "DPS_template", "syst", sigma_dps_cut * dps_method_relative, sigma_dps_cut * dps_method_relative, "pb", "included"],
    ]
    for name, relative in pair_sources.items():
        component_rows.append([
            "sigma_DPS_mJJ_ge_7p5",
            f"sigma_total_{name}",
            "syst",
            sigma_dps_cut * relative,
            sigma_dps_cut * relative,
            "pb",
            "included",
        ])
    component_rows.extend([
        ["sigma_DPS_mJJ_ge_7p5", "sigma_total_fitter_stability", "syst", "", "", "pb", "pending_not_included_not_zero"],
        ["sigma_eff", "sigma_DPS", "stat", sigma_eff_stat_from_dps, sigma_eff_stat_from_dps, "mb", "included"],
        ["sigma_eff", "ATLAS_single_Jpsi", "stat", sigma_eff_stat_from_single, sigma_eff_stat_from_single, "mb", "included"],
        ["sigma_eff", "data_driven_method", "syst_endpoint", sigma_eff_method_up, sigma_eff_method_down, "mb", "included"],
        ["sigma_eff", "sigma_total_non_BR", "syst_endpoint", sigma_eff_pair_up, sigma_eff_pair_down, "mb", "included_fitter_pending"],
        ["sigma_eff", "pair_branching_fraction", "syst", 0.0, 0.0, "mb", "cancelled_exactly"],
        ["sigma_eff", "ATLAS_single_bin", "syst_endpoint", sigma_eff_single_bin_up, sigma_eff_single_bin_down, "mb", "included"],
        ["sigma_eff", "ATLAS_luminosity", "syst_endpoint", sigma_eff_single_lumi_up, sigma_eff_single_lumi_down, "mb", "included"],
        ["sigma_eff", "sigma_total_fitter_stability", "syst", "", "", "mb", "pending_not_included_not_zero"],
        ["sigma_eff", "epsilon_m", "stat_or_model", "", "", "mb", "pending_not_included_not_zero"],
    ])

    output_path.mkdir(parents=True, exist_ok=False)
    (output_path / "inputs.json").write_text(
        json.dumps(inputs, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output_path / "results.json").write_text(
        json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    with (output_path / "results.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "quantity", "phase_space", "central", "stat", "stat_percent",
            "included_syst_up", "included_syst_up_percent",
            "included_syst_down", "included_syst_down_percent", "unit",
            "final_systematic_status",
        ])
        writer.writerow([
            "f_DPS", "mJJ>=7.5", f_dps, f_stat, percent(f_stat, f_dps),
            f_syst_up, percent(f_syst_up, f_dps), f_syst_down,
            percent(f_syst_down, f_dps), "fraction",
            "complete_excludes_sigma_total_systematics_by_definition",
        ])
        writer.writerow([
            "sigma_DPS", "mJJ>=7.5", sigma_dps_cut, sigma_dps_cut_stat,
            percent(sigma_dps_cut_stat, sigma_dps_cut), sigma_dps_cut_syst_up,
            percent(sigma_dps_cut_syst_up, sigma_dps_cut),
            sigma_dps_cut_syst_down,
            percent(sigma_dps_cut_syst_down, sigma_dps_cut), "pb",
            "pending_fitter_stability_not_zero",
        ])
        writer.writerow([
            "sigma_DPS", "no_mJJ_cut", sigma_dps_no_mass,
            sigma_dps_no_mass_stat, percent(sigma_dps_no_mass_stat, sigma_dps_no_mass),
            sigma_dps_no_mass_syst_up,
            percent(sigma_dps_no_mass_syst_up, sigma_dps_no_mass),
            sigma_dps_no_mass_syst_down,
            percent(sigma_dps_no_mass_syst_down, sigma_dps_no_mass), "pb",
            "pending_fitter_stability_and_epsilon_m_uncertainty_not_zero",
        ])
        writer.writerow([
            "sigma_eff", "single_Jpsi_matched_no_pair_mass_cut", sigma_eff_mb,
            sigma_eff_stat, percent(sigma_eff_stat, sigma_eff_mb),
            sigma_eff_syst_up, percent(sigma_eff_syst_up, sigma_eff_mb),
            sigma_eff_syst_down, percent(sigma_eff_syst_down, sigma_eff_mb),
            "mb", "pending_fitter_stability_and_epsilon_m_uncertainty_not_zero",
        ])

    with (output_path / "uncertainty_components.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "quantity", "source", "type", "effect_up", "effect_down",
            "unit", "treatment", "effect_up_percent", "effect_down_percent",
        ])
        for row in component_rows:
            if row[0] == "f_DPS":
                central = f_dps
            elif row[0] == "sigma_eff":
                central = sigma_eff_mb
            else:
                central = sigma_dps_cut
            up_percent = percent(float(row[3]), central) if row[3] != "" else ""
            down_percent = percent(float(row[4]), central) if row[4] != "" else ""
            writer.writerow(row + [up_percent, down_percent])

    summary_lines = [
        "status=candidate_complete_included_sources_pending_fitter_stability_and_epsilon_m_uncertainty",
        "f_dps_status=nominal_updated_user_defined_yield_based_scheme",
        f"control_dps_yield={data_control:.12g}",
        f"control_dps_yield_error={data_control_error:.12g}",
        f"control_dps_yield_stat_percent={100.0 * control_stat_relative:.12g}",
        f"full_fiducial_dps_yield={dps_yield:.12g}",
        f"full_fiducial_dps_yield_stat={dps_yield_stat:.12g}",
        f"total_pair_yield={total_yield:.12g}",
        f"total_pair_yield_stat={total_yield_error:.12g}",
        f"pair_total_pb={pair_pb:.12g}",
        f"pair_total_stat_pb={pair_stat_pb:.12g}",
        f"pair_total_known_syst_relative={pair_known_relative:.12g}",
        f"pair_total_known_syst_pb={pair_pb * pair_known_relative:.12g}",
        f"pair_total_transferred_fitter_relative_not_used={transferred_fitter_relative:.12g}",
        "pair_total_fitter_stability_status=pending_not_included_not_zero",
        f"f_dps={f_dps:.12g}",
        f"f_dps_stat={f_stat:.12g}",
        f"f_dps_stat_percent={percent(f_stat, f_dps):.12g}",
        f"f_dps_stat_from_control_yield={f_stat_control:.12g}",
        f"f_dps_stat_from_control_yield_percent={percent(f_stat_control, f_dps):.12g}",
        f"f_dps_stat_from_sigma_total={f_stat_total:.12g}",
        f"f_dps_stat_from_sigma_total_percent={percent(f_stat_total, f_dps):.12g}",
        f"f_dps_cr_syst_up={f_cr_up:.12g}",
        f"f_dps_cr_syst_down={f_cr_down:.12g}",
        f"f_dps_cr_syst_up_percent={percent(f_cr_up, f_dps):.12g}",
        f"f_dps_cr_syst_down_percent={percent(f_cr_down, f_dps):.12g}",
        f"f_dps_dps_template_syst={f_dps_method:.12g}",
        f"f_dps_dps_template_syst_percent={percent(f_dps_method, f_dps):.12g}",
        f"f_dps_total_syst_up={f_syst_up:.12g}",
        f"f_dps_total_syst_down={f_syst_down:.12g}",
        f"f_dps_total_syst_up_percent={percent(f_syst_up, f_dps):.12g}",
        f"f_dps_total_syst_down_percent={percent(f_syst_down, f_dps):.12g}",
        "f_dps_sigma_total_systematics=not_included_by_user_definition",
        f"epsilon_m={epsilon_m:.15g}",
        "epsilon_m_uncertainty_status=pending_not_included_not_zero",
        f"sigma_dps_mjj_ge_7p5_pb={sigma_dps_cut:.12g}",
        f"sigma_dps_mjj_ge_7p5_stat_pb={sigma_dps_cut_stat:.12g}",
        f"sigma_dps_mjj_ge_7p5_stat_percent={percent(sigma_dps_cut_stat, sigma_dps_cut):.12g}",
        f"sigma_dps_mjj_ge_7p5_included_syst_up_pb={sigma_dps_cut_syst_up:.12g}",
        f"sigma_dps_mjj_ge_7p5_included_syst_down_pb={sigma_dps_cut_syst_down:.12g}",
        f"sigma_dps_mjj_ge_7p5_included_syst_up_percent={percent(sigma_dps_cut_syst_up, sigma_dps_cut):.12g}",
        f"sigma_dps_mjj_ge_7p5_included_syst_down_percent={percent(sigma_dps_cut_syst_down, sigma_dps_cut):.12g}",
        "sigma_dps_total_systematic_status=pending_fitter_stability_not_included_not_zero",
        f"sigma_dps_no_mjj_cut_pb={sigma_dps_no_mass:.12g}",
        f"sigma_dps_no_mjj_cut_stat_pb={sigma_dps_no_mass_stat:.12g}",
        f"sigma_dps_no_mjj_cut_included_syst_up_pb={sigma_dps_no_mass_syst_up:.12g}",
        f"sigma_dps_no_mjj_cut_included_syst_down_pb={sigma_dps_no_mass_syst_down:.12g}",
        f"sigma_eff_mb={sigma_eff_mb:.12g}",
        f"sigma_eff_stat_mb={sigma_eff_stat:.12g}",
        f"sigma_eff_stat_percent={percent(sigma_eff_stat, sigma_eff_mb):.12g}",
        f"sigma_eff_included_syst_up_mb={sigma_eff_syst_up:.12g}",
        f"sigma_eff_included_syst_down_mb={sigma_eff_syst_down:.12g}",
        f"sigma_eff_included_syst_up_percent={percent(sigma_eff_syst_up, sigma_eff_mb):.12g}",
        f"sigma_eff_included_syst_down_percent={percent(sigma_eff_syst_down, sigma_eff_mb):.12g}",
        "sigma_eff_total_systematic_status=pending_fitter_stability_and_epsilon_m_uncertainty_not_included_not_zero",
        "half_vs_full_selfmix_shift_role=cross_check_not_formal_statistical_uncertainty",
        "new_roofit_fits_run=0",
    ]
    (output_path / "summary.txt").write_text(
        "\n".join(summary_lines) + "\n", encoding="utf-8"
    )

    generated_at = datetime.now().astimezone().isoformat(timespec="seconds")
    metadata_lines = [
        "status=complete_candidate_pending_fitter_stability_and_epsilon_m_uncertainty",
        f"output_tag={output_path.name}",
        f"generated_at={generated_at}",
        f"root_version={root_version}",
        f"git_head={git_head}",
        f"script_sha256={sha256(Path(__file__).resolve())}",
        "record_type=deterministic_postprocess_no_ROOT_fits",
        "central_definition=full_fiducial_DPS_yield_divided_by_accepted_total_fit_yield",
        "f_dps_stat_definition=CR_DPS_yield_and_sigma_total_stat_independent_quadrature",
        "sigma_dps_stat_definition=CR_DPS_yield_only",
        "f_dps_syst_definition=CR_boundary_and_full_direct_DPS_difference_no_sigma_total_systematics",
        "sigma_dps_syst_definition=two_method_sources_plus_known_sigma_total_sources_fitter_pending",
        "sigma_eff_definition=original_inverse_sigma_DPS_propagation_with_ATLAS_single_Jpsi_inputs",
        f"transferred_fitter_relative_not_used={transferred_fitter_relative:.12g}",
        "updated_fitter_stability=pending_not_included_not_zero",
        "epsilon_m_uncertainty=pending_not_included_not_zero",
        "user_gate=2026-09-03_yield_based_uncertainty_update_and_leave_Data_fitter_stability_pending",
    ]
    (output_path / "run.metadata.txt").write_text(
        "\n".join(metadata_lines) + "\n", encoding="utf-8"
    )
    source_lines = [
        f"{sha256(path)}  {path.resolve()}" for path in input_paths.values()
    ]
    (output_path / "source_checksums.txt").write_text(
        "\n".join(source_lines) + "\n", encoding="utf-8"
    )
    (output_path / "complete.marker").write_text(
        "status=complete_candidate_pending_fitter_stability_and_epsilon_m_uncertainty\n",
        encoding="utf-8",
    )
    artifact_names = [
        "inputs.json",
        "results.json",
        "results.csv",
        "uncertainty_components.csv",
        "summary.txt",
        "run.metadata.txt",
        "source_checksums.txt",
        "complete.marker",
    ]
    artifact_lines = [
        f"{sha256(output_path / name)}  {name}" for name in artifact_names
    ]
    (output_path / "artifact_checksums.txt").write_text(
        "\n".join(artifact_lines) + "\n", encoding="utf-8"
    )

    print(
        f"f_DPS = {f_dps:.12g} +/- {f_stat:.12g} (stat.) "
        f"+{f_syst_up:.12g}/-{f_syst_down:.12g} (syst.; no sigma_tot syst.)"
    )
    print(
        f"sigma_DPS(mJJ>=7.5) = {sigma_dps_cut:.12g} +/- "
        f"{sigma_dps_cut_stat:.12g} (stat.) +{sigma_dps_cut_syst_up:.12g}/-"
        f"{sigma_dps_cut_syst_down:.12g} (included syst.; fitter pending) pb"
    )
    print(
        f"sigma_eff = {sigma_eff_mb:.12g} +/- {sigma_eff_stat:.12g} (stat.) "
        f"+{sigma_eff_syst_up:.12g}/-{sigma_eff_syst_down:.12g} "
        "(included syst.; fitter and epsilon_m pending) mb"
    )


if __name__ == "__main__":
    main()
