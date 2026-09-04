#!/usr/bin/env python3
"""Propagate an accepted Data fitter systematic into the Data-driven nominal.

This deterministic postprocess consumes the immutable yield-based candidate
whose statistical definition already includes the CR-yield and half-vs-full
self-mixing template terms.  It does not run ROOT/RooFit fits.

The accepted fitter variation is a total-cross-section systematic, so it:

* does not enter the f_DPS method systematic;
* enters sigma_DPS in quadrature with the other total-cross-section sources;
* enters sigma_eff through the non-BR pair-cross-section denominator term.

The epsilon_m uncertainty remains explicit and pending.  It is never treated
as zero by this script.
"""

from __future__ import annotations

import copy
import csv
import hashlib
import json
import math
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path


EXPECTED_ROOT_VERSION = "6.40.02"
EXPECTED_REQUEST = "H018_recalculate_fitter_stability_with_legacy_variations"
EXPECTED_RUN_TAG = (
    "h018_fitter_stability_accmix23_effmix19_sps0p85_dps0p15_"
    "root640_v8_20260903"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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


def quadrature(*values: float) -> float:
    return math.sqrt(sum(value * value for value in values))


def percent(value: float, central: float) -> float:
    return 100.0 * value / central


def inverse_up(central: float, relative_down: float) -> float:
    return central / (1.0 - relative_down) - central


def inverse_down(central: float, relative_up: float) -> float:
    return central - central / (1.0 + relative_up)


def require_close(name: str, actual: float, expected: float) -> None:
    if not math.isclose(actual, expected, rel_tol=5e-11, abs_tol=1e-12):
        raise RuntimeError(f"{name} mismatch: {actual} vs {expected}")


def verify_checksums(manifest: Path) -> None:
    for line in manifest.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        expected, name = line.split(None, 1)
        path = Path(name.lstrip(" *"))
        if not path.is_absolute():
            path = manifest.parent / path
        actual = sha256(path)
        if actual != expected:
            raise RuntimeError(
                f"checksum mismatch for {path}: {actual} != {expected}"
            )


def write_json(path: Path, value: object) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    if len(sys.argv) != 6:
        raise SystemExit(
            "usage: finalize_yield_based_nominal_with_fitter.py "
            "BASE_CANDIDATE_DIR FITTER_METADATA FITTER_RECOMBINATION_CSV "
            "FITTER_CHECKSUMS OUTPUT_DIR"
        )

    base_dir, fitter_metadata_path, fitter_csv_path, fitter_checksums, output_dir = (
        map(Path, sys.argv[1:])
    )
    output_dir = output_dir.resolve()
    if output_dir.exists():
        raise RuntimeError(f"refusing to overwrite existing output: {output_dir}")

    base_files = {
        "inputs": base_dir / "inputs.json",
        "results": base_dir / "results.json",
        "results_csv": base_dir / "results.csv",
        "uncertainty_components": base_dir / "uncertainty_components.csv",
        "summary": base_dir / "summary.txt",
        "metadata": base_dir / "run.metadata.txt",
        "artifact_checksums": base_dir / "artifact_checksums.txt",
    }
    source_paths = {
        **base_files,
        "fitter_metadata": fitter_metadata_path,
        "fitter_recombination": fitter_csv_path,
        "fitter_checksums": fitter_checksums,
        "finalization_script": Path(__file__).resolve(),
    }
    for label, path in source_paths.items():
        if not path.is_file():
            raise RuntimeError(f"missing {label}: {path}")

    verify_checksums(base_files["artifact_checksums"])
    verify_checksums(fitter_checksums)

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

    base_inputs = json.loads(base_files["inputs"].read_text(encoding="utf-8"))
    base_results = json.loads(base_files["results"].read_text(encoding="utf-8"))
    base_metadata = read_kv(base_files["metadata"])
    fitter_metadata = read_kv(fitter_metadata_path)

    expected_base_status = (
        "candidate_complete_included_sources_pending_fitter_stability_"
        "and_epsilon_m_uncertainty"
    )
    if base_results.get("status") != expected_base_status:
        raise RuntimeError("base candidate has unexpected result status")
    if base_metadata.get("updated_fitter_stability") != (
        "pending_not_included_not_zero"
    ):
        raise RuntimeError("base candidate does not leave fitter stability pending")
    if fitter_metadata.get("status") != "complete":
        raise RuntimeError("fitter campaign is incomplete")
    if fitter_metadata.get("request") != EXPECTED_REQUEST:
        raise RuntimeError("unexpected fitter request")
    if fitter_metadata.get("run_tag") != EXPECTED_RUN_TAG:
        raise RuntimeError("unexpected fitter run tag")
    if fitter_metadata.get("root_version") != EXPECTED_ROOT_VERSION:
        raise RuntimeError("fitter campaign ROOT version mismatch")
    if fitter_metadata.get("differential_reference_selection") != "user_confirmed":
        raise RuntimeError("fitter differential reference is not user confirmed")

    inherited_inputs = base_inputs["inherited_inputs"]
    inherited_pair = inherited_inputs["pair_total"]
    accepted_table_hash = inherited_inputs["source_sha256"]["data_cross_sections"]
    if fitter_metadata.get("nominal_table_sha256") != accepted_table_hash:
        raise RuntimeError("fitter nominal table does not match Data-driven input")

    with fitter_csv_path.open(newline="", encoding="utf-8") as handle:
        total_rows = [
            row for row in csv.DictReader(handle) if row["scope"] == "total"
        ]
    if len(total_rows) != 1:
        raise RuntimeError("fitter recombination table does not have one total row")
    fitter_total = total_rows[0]
    fitter_relative = float(fitter_total["fitter_relative"])
    if fitter_relative <= 0.0:
        raise RuntimeError("accepted fitter relative uncertainty is not positive")

    require_close(
        "pair central",
        float(fitter_total["cross_section_pb_per_unit"]),
        float(inherited_pair["central_pb"]),
    )
    require_close(
        "pair stat",
        float(fitter_total["stat_pb_per_unit"]),
        float(inherited_pair["stat_pb"]),
    )
    pair_sources = dict(inherited_pair["known_systematic_sources_relative"])
    source_columns = {
        "acceptance_efficiency": "correction_relative",
        "branching_fraction": "br_relative",
        "cms_luminosity": "luminosity_relative",
        "lifetime_variable": "lifetime_relative",
    }
    for name, column in source_columns.items():
        require_close(name, float(fitter_total[column]), float(pair_sources[name]))

    pair_total_relative = quadrature(*pair_sources.values(), fitter_relative)
    require_close(
        "accepted total systematic",
        pair_total_relative,
        float(fitter_total["total_systematic_relative"]),
    )
    pair_non_br_relative = quadrature(
        pair_sources["acceptance_efficiency"],
        pair_sources["cms_luminosity"],
        pair_sources["lifetime_variable"],
        fitter_relative,
    )
    pair_non_br_without_fitter = quadrature(
        pair_sources["acceptance_efficiency"],
        pair_sources["cms_luminosity"],
        pair_sources["lifetime_variable"],
    )

    revised = copy.deepcopy(base_results)
    revised["status"] = (
        "confirmed_nominal_user_accepted_fitter_v8_"
        "epsilon_m_uncertainty_pending_not_zero"
    )
    revised["assumptions"]["fitter_stability"] = (
        "v8_user_accepted_included_in_sigma_DPS_and_sigma_eff"
    )
    revised["assumptions"]["epsilon_m_uncertainty"] = (
        "pending_not_included_not_zero"
    )
    revised["assumptions"]["nominal_promotion"] = (
        "2026-09-03_user_accepted_fitter_v8_and_closed_current_method_as_nominal"
    )

    f_result = revised["f_dps"]
    f_result["status"] = "confirmed_nominal_complete"
    method_up_relative = float(f_result["combined_syst_up_percent"]) / 100.0
    method_down_relative = float(f_result["combined_syst_down_percent"]) / 100.0

    sigma_cut = revised["sigma_dps_mjj_ge_7p5_pb"]
    sigma_cut_central = float(sigma_cut["central"])
    sigma_cut_up_relative = quadrature(method_up_relative, pair_total_relative)
    sigma_cut_down_relative = quadrature(method_down_relative, pair_total_relative)
    sigma_cut_up = sigma_cut_central * sigma_cut_up_relative
    sigma_cut_down = sigma_cut_central * sigma_cut_down_relative
    sigma_cut.update(
        {
            "included_syst_up": sigma_cut_up,
            "included_syst_down": sigma_cut_down,
            "included_syst_up_percent": 100.0 * sigma_cut_up_relative,
            "included_syst_down_percent": 100.0 * sigma_cut_down_relative,
            "total_syst": {"up": sigma_cut_up, "down": sigma_cut_down},
            "total_syst_status": "complete_all_defined_sources_included",
        }
    )

    sigma_no_mass = revised["sigma_dps_no_mjj_cut_pb"]
    sigma_no_mass_central = float(sigma_no_mass["central"])
    sigma_no_mass_up = sigma_no_mass_central * sigma_cut_up_relative
    sigma_no_mass_down = sigma_no_mass_central * sigma_cut_down_relative
    sigma_no_mass.update(
        {
            "included_syst_up": sigma_no_mass_up,
            "included_syst_down": sigma_no_mass_down,
            "included_syst_up_percent": 100.0 * sigma_cut_up_relative,
            "included_syst_down_percent": 100.0 * sigma_cut_down_relative,
            "total_syst": None,
            "total_syst_status": "pending_epsilon_m_uncertainty_not_zero",
        }
    )

    with base_files["uncertainty_components"].open(
        newline="", encoding="utf-8"
    ) as handle:
        component_reader = csv.DictReader(handle)
        component_fields = list(component_reader.fieldnames or [])
        component_rows = list(component_reader)

    def find_component(quantity: str, source: str) -> dict[str, str]:
        matches = [
            row
            for row in component_rows
            if row["quantity"] == quantity and row["source"] == source
        ]
        if len(matches) != 1:
            raise RuntimeError(f"expected one component {quantity}/{source}")
        return matches[0]

    sigma_cut_fitter = find_component(
        "sigma_DPS_mJJ_ge_7p5", "sigma_total_fitter_stability"
    )
    sigma_cut_fitter_effect = sigma_cut_central * fitter_relative
    sigma_cut_fitter.update(
        {
            "effect_up": str(sigma_cut_fitter_effect),
            "effect_down": str(sigma_cut_fitter_effect),
            "treatment": "included_user_accepted_v8",
            "effect_up_percent": str(100.0 * fitter_relative),
            "effect_down_percent": str(100.0 * fitter_relative),
        }
    )

    sigma_eff = revised["sigma_eff_mb"]
    sigma_eff_central = float(sigma_eff["central"])
    pair_component = find_component("sigma_eff", "sigma_total_non_BR")
    require_close(
        "base sigma_eff pair up",
        float(pair_component["effect_up"]),
        inverse_up(sigma_eff_central, pair_non_br_without_fitter),
    )
    require_close(
        "base sigma_eff pair down",
        float(pair_component["effect_down"]),
        inverse_down(sigma_eff_central, pair_non_br_without_fitter),
    )
    pair_effect_up = inverse_up(sigma_eff_central, pair_non_br_relative)
    pair_effect_down = inverse_down(sigma_eff_central, pair_non_br_relative)
    pair_component.update(
        {
            "effect_up": str(pair_effect_up),
            "effect_down": str(pair_effect_down),
            "treatment": "included_with_user_accepted_fitter_v8",
            "effect_up_percent": str(percent(pair_effect_up, sigma_eff_central)),
            "effect_down_percent": str(percent(pair_effect_down, sigma_eff_central)),
        }
    )
    fitter_eff_component = find_component(
        "sigma_eff", "sigma_total_fitter_stability"
    )
    fitter_eff_up = inverse_up(sigma_eff_central, fitter_relative)
    fitter_eff_down = inverse_down(sigma_eff_central, fitter_relative)
    fitter_eff_component.update(
        {
            "effect_up": str(fitter_eff_up),
            "effect_down": str(fitter_eff_down),
            "treatment": "included_in_sigma_total_non_BR_quadrature_before_endpoint",
            "effect_up_percent": str(percent(fitter_eff_up, sigma_eff_central)),
            "effect_down_percent": str(percent(fitter_eff_down, sigma_eff_central)),
        }
    )

    method_component = find_component("sigma_eff", "data_driven_method")
    atlas_bin_component = find_component("sigma_eff", "ATLAS_single_bin")
    atlas_lumi_component = find_component("sigma_eff", "ATLAS_luminosity")
    sigma_eff_up = quadrature(
        float(method_component["effect_up"]),
        pair_effect_up,
        float(atlas_bin_component["effect_up"]),
        float(atlas_lumi_component["effect_up"]),
    )
    sigma_eff_down = quadrature(
        float(method_component["effect_down"]),
        pair_effect_down,
        float(atlas_bin_component["effect_down"]),
        float(atlas_lumi_component["effect_down"]),
    )
    sigma_eff.update(
        {
            "included_syst_up": sigma_eff_up,
            "included_syst_down": sigma_eff_down,
            "included_syst_up_percent": percent(sigma_eff_up, sigma_eff_central),
            "included_syst_down_percent": percent(
                sigma_eff_down, sigma_eff_central
            ),
            "total_syst": None,
            "total_syst_status": "pending_epsilon_m_uncertainty_not_zero",
        }
    )

    component_rows.extend(
        [
            {
                "quantity": "sigma_DPS_mJJ_ge_7p5",
                "source": "combined",
                "type": "syst_total",
                "effect_up": sigma_cut_up,
                "effect_down": sigma_cut_down,
                "unit": "pb",
                "treatment": "complete_all_defined_sources_in_quadrature",
                "effect_up_percent": 100.0 * sigma_cut_up_relative,
                "effect_down_percent": 100.0 * sigma_cut_down_relative,
            },
            {
                "quantity": "sigma_eff",
                "source": "combined_quantified",
                "type": "syst_subtotal",
                "effect_up": sigma_eff_up,
                "effect_down": sigma_eff_down,
                "unit": "mb",
                "treatment": "epsilon_m_uncertainty_pending_not_included_not_zero",
                "effect_up_percent": percent(sigma_eff_up, sigma_eff_central),
                "effect_down_percent": percent(sigma_eff_down, sigma_eff_central),
            },
        ]
    )

    with base_files["results_csv"].open(newline="", encoding="utf-8") as handle:
        result_reader = csv.DictReader(handle)
        result_fields = list(result_reader.fieldnames or [])
        result_rows = list(result_reader)
    result_updates = {
        ("f_DPS", "mJJ>=7.5"): (
            float(f_result["combined_syst_up"]),
            float(f_result["combined_syst_down"]),
            "confirmed_nominal_complete_excludes_sigma_total_systematics_by_definition",
        ),
        ("sigma_DPS", "mJJ>=7.5"): (
            sigma_cut_up,
            sigma_cut_down,
            "confirmed_nominal_complete_all_defined_sources_included",
        ),
        ("sigma_DPS", "no_mJJ_cut"): (
            sigma_no_mass_up,
            sigma_no_mass_down,
            "confirmed_nominal_quantified_sources_complete_epsilon_m_uncertainty_pending_not_zero",
        ),
        ("sigma_eff", "single_Jpsi_matched_no_pair_mass_cut"): (
            sigma_eff_up,
            sigma_eff_down,
            "confirmed_nominal_quantified_sources_complete_epsilon_m_uncertainty_pending_not_zero",
        ),
    }
    for row in result_rows:
        key = (row["quantity"], row["phase_space"])
        if key not in result_updates:
            raise RuntimeError(f"unexpected result row {key}")
        up, down, status = result_updates[key]
        central = float(row["central"])
        row["included_syst_up"] = up
        row["included_syst_up_percent"] = percent(up, central)
        row["included_syst_down"] = down
        row["included_syst_down_percent"] = percent(down, central)
        row["final_systematic_status"] = status

    source_hashes = {name: sha256(path) for name, path in source_paths.items()}
    revised_inputs = {
        "status": "frozen_snapshot_for_confirmed_nominal_fitter_v8_propagation",
        "base_candidate": {
            "path": str(base_dir.resolve()),
            "artifact_checksums_verified": True,
        },
        "fitter_stability": {
            "status": "accepted_by_user_and_included",
            "run_tag": fitter_metadata["run_tag"],
            "relative": fitter_relative,
            "percent": 100.0 * fitter_relative,
            "envelope": fitter_metadata["differential_reference"],
            "official_checksums_verified": True,
            "total_systematic_relative_after_recombination": pair_total_relative,
        },
        "epsilon_m_uncertainty": {
            "value": None,
            "status": "pending_not_included_not_zero",
        },
        "inherited_inputs": base_inputs,
        "source_paths": {
            name: str(path.resolve()) for name, path in source_paths.items()
        },
        "source_sha256": source_hashes,
    }

    generated_at = datetime.now().astimezone().isoformat(timespec="seconds")
    summary_lines = [
        f"status={revised['status']}",
        "promotion_status=confirmed_nominal_user_accepted",
        f"fitter_stability_source={fitter_metadata['run_tag']}",
        f"fitter_stability_relative={fitter_relative:.12g}",
        f"fitter_stability_percent={100.0 * fitter_relative:.12g}",
        f"sigma_total_systematic_relative={pair_total_relative:.12g}",
        f"sigma_total_systematic_percent={100.0 * pair_total_relative:.12g}",
        f"f_dps={float(f_result['central']):.12g}",
        f"f_dps_stat={float(f_result['stat']):.12g}",
        f"f_dps_stat_percent={float(f_result['stat_percent']):.12g}",
        f"f_dps_total_syst_up={float(f_result['combined_syst_up']):.12g}",
        f"f_dps_total_syst_down={float(f_result['combined_syst_down']):.12g}",
        f"f_dps_total_syst_up_percent={float(f_result['combined_syst_up_percent']):.12g}",
        f"f_dps_total_syst_down_percent={float(f_result['combined_syst_down_percent']):.12g}",
        "f_dps_sigma_total_systematics=not_included_by_user_definition",
        f"sigma_dps_mjj_ge_7p5_pb={sigma_cut_central:.12g}",
        f"sigma_dps_mjj_ge_7p5_stat_pb={float(sigma_cut['stat']):.12g}",
        f"sigma_dps_mjj_ge_7p5_stat_percent={float(sigma_cut['stat_percent']):.12g}",
        f"sigma_dps_mjj_ge_7p5_syst_up_pb={sigma_cut_up:.12g}",
        f"sigma_dps_mjj_ge_7p5_syst_down_pb={sigma_cut_down:.12g}",
        f"sigma_dps_mjj_ge_7p5_syst_up_percent={100.0 * sigma_cut_up_relative:.12g}",
        f"sigma_dps_mjj_ge_7p5_syst_down_percent={100.0 * sigma_cut_down_relative:.12g}",
        "sigma_dps_mjj_ge_7p5_systematic_status=complete_all_defined_sources_included",
        f"sigma_dps_no_mjj_cut_pb={sigma_no_mass_central:.12g}",
        f"sigma_dps_no_mjj_cut_stat_pb={float(sigma_no_mass['stat']):.12g}",
        f"sigma_dps_no_mjj_cut_syst_up_pb={sigma_no_mass_up:.12g}",
        f"sigma_dps_no_mjj_cut_syst_down_pb={sigma_no_mass_down:.12g}",
        f"sigma_dps_no_mjj_cut_syst_up_percent={100.0 * sigma_cut_up_relative:.12g}",
        f"sigma_dps_no_mjj_cut_syst_down_percent={100.0 * sigma_cut_down_relative:.12g}",
        "sigma_dps_no_mjj_cut_systematic_status=epsilon_m_uncertainty_pending_not_included_not_zero",
        f"sigma_eff_mb={sigma_eff_central:.12g}",
        f"sigma_eff_stat_mb={float(sigma_eff['stat']):.12g}",
        f"sigma_eff_stat_percent={float(sigma_eff['stat_percent']):.12g}",
        f"sigma_eff_quantified_syst_up_mb={sigma_eff_up:.12g}",
        f"sigma_eff_quantified_syst_down_mb={sigma_eff_down:.12g}",
        f"sigma_eff_quantified_syst_up_percent={percent(sigma_eff_up, sigma_eff_central):.12g}",
        f"sigma_eff_quantified_syst_down_percent={percent(sigma_eff_down, sigma_eff_central):.12g}",
        "sigma_eff_systematic_status=epsilon_m_uncertainty_pending_not_included_not_zero",
        "epsilon_m_uncertainty_status=pending_not_included_not_zero",
        "new_roofit_fits_run=0",
    ]
    metadata_lines = [
        "status=complete_confirmed_nominal_epsilon_m_uncertainty_pending_not_zero",
        "promotion_status=accepted_by_user",
        f"output_tag={output_dir.name}",
        f"generated_at={generated_at}",
        f"root_version={root_version}",
        f"git_head={git_head}",
        f"script_sha256={source_hashes['finalization_script']}",
        "record_type=deterministic_fitter_v8_propagation_no_ROOT_fits",
        "base_candidate_artifacts_verified=true",
        "fitter_official_checksums_verified=true",
        f"fitter_run_tag={fitter_metadata['run_tag']}",
        f"fitter_relative={fitter_relative:.12g}",
        "fitter_treatment=excluded_from_f_DPS_included_in_sigma_DPS_and_sigma_eff",
        "sigma_dps_stat_definition=CR_DPS_yield_and_half_vs_full_selfmix_template_independent_quadrature",
        "f_dps_stat_definition=sigma_DPS_stat_and_sigma_total_stat_independent_quadrature",
        "epsilon_m_uncertainty=pending_not_included_not_zero",
        "user_gate=2026-09-03_accept_fitter_v8_and_close_current_Data_driven_method_as_nominal",
    ]

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(prefix=f".{output_dir.name}.tmp.", dir=output_dir.parent)
    )
    try:
        write_json(staging / "inputs.json", revised_inputs)
        write_json(staging / "results.json", revised)
        write_csv(staging / "results.csv", result_fields, result_rows)
        write_csv(
            staging / "uncertainty_components.csv", component_fields, component_rows
        )
        (staging / "summary.txt").write_text(
            "\n".join(summary_lines) + "\n", encoding="utf-8"
        )
        (staging / "run.metadata.txt").write_text(
            "\n".join(metadata_lines) + "\n", encoding="utf-8"
        )
        (staging / "source_checksums.txt").write_text(
            "".join(
                f"{source_hashes[name]}  {source_paths[name].resolve()}\n"
                for name in sorted(source_paths)
            ),
            encoding="utf-8",
        )
        (staging / "complete.marker").write_text(
            "status=complete_confirmed_nominal_epsilon_m_uncertainty_pending_not_zero\n",
            encoding="utf-8",
        )
        artifacts = [
            "inputs.json",
            "results.json",
            "results.csv",
            "uncertainty_components.csv",
            "summary.txt",
            "run.metadata.txt",
            "source_checksums.txt",
            "complete.marker",
        ]
        (staging / "artifact_checksums.txt").write_text(
            "".join(f"{sha256(staging / name)}  {name}\n" for name in artifacts),
            encoding="utf-8",
        )
        if output_dir.exists():
            raise RuntimeError(f"refusing to overwrite existing output: {output_dir}")
        staging.rename(output_dir)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise

    print(f"created {output_dir}")
    print(
        f"f_DPS={float(f_result['central']):.12g} +/- "
        f"{float(f_result['stat']):.12g} stat "
        f"+{float(f_result['combined_syst_up']):.12g}/"
        f"-{float(f_result['combined_syst_down']):.12g} method syst"
    )
    print(
        f"sigma_DPS(mJJ>=7.5)={sigma_cut_central:.12g} +/- "
        f"{float(sigma_cut['stat']):.12g} stat +{sigma_cut_up:.12g}/"
        f"-{sigma_cut_down:.12g} syst pb"
    )
    print(
        f"sigma_eff={sigma_eff_central:.12g} +/- "
        f"{float(sigma_eff['stat']):.12g} stat +{sigma_eff_up:.12g}/"
        f"-{sigma_eff_down:.12g} quantified syst mb; epsilon_m pending"
    )


if __name__ == "__main__":
    main()
