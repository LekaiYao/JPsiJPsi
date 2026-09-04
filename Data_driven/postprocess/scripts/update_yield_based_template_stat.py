#!/usr/bin/env python3
"""Add the accepted half-vs-full self-mixing fluctuation to yield-based stats.

This is a deterministic revision of an immutable yield-based candidate.  It
does not rerun ROOT/RooFit.  The revised definitions are:

* sigma_DPS stat = CR DPS-yield stat (+) half-vs-full template stat;
* f_DPS stat = revised sigma_DPS stat (+) total-cross-section stat;
* sigma_eff stat keeps the original propagation from sigma_DPS and the
  external single-J/psi statistical uncertainty.

Here (+) denotes independent quadrature.  Systematic uncertainties and all
central values are copied unchanged.  Data fitter stability and epsilon_m
uncertainty remain explicitly pending.
"""

from __future__ import annotations

import copy
import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path


EXPECTED_ROOT_VERSION = "6.40.02"


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


def require_close(name: str, actual: float, expected: float) -> None:
    if not math.isclose(actual, expected, rel_tol=5e-11, abs_tol=1e-12):
        raise RuntimeError(f"{name} mismatch: {actual} vs {expected}")


def verify_checksums(manifest: Path) -> None:
    for line in manifest.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        expected, name = line.split(None, 1)
        name = name.lstrip(" *")
        path = manifest.parent / name
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


def effect_row(
    quantity: str,
    source: str,
    kind: str,
    effect: float,
    unit: str,
    treatment: str,
    central: float,
) -> dict[str, object]:
    return {
        "quantity": quantity,
        "source": source,
        "type": kind,
        "effect_up": effect,
        "effect_down": effect,
        "unit": unit,
        "treatment": treatment,
        "effect_up_percent": percent(effect, central),
        "effect_down_percent": percent(effect, central),
    }


def main() -> None:
    if len(sys.argv) != 4:
        raise SystemExit(
            "usage: update_yield_based_template_stat.py "
            "BASE_CANDIDATE_DIR TEMPLATE_STAT_SUMMARY OUTPUT_DIR"
        )

    base_dir, template_summary_path, output_dir = map(Path, sys.argv[1:])
    output_dir = output_dir.resolve()
    if output_dir.exists():
        raise RuntimeError(f"refusing to overwrite existing output: {output_dir}")
    if not base_dir.is_dir():
        raise RuntimeError(f"missing base candidate: {base_dir}")
    if not template_summary_path.is_file():
        raise RuntimeError(f"missing template-stat summary: {template_summary_path}")

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

    base_files = {
        "inputs": base_dir / "inputs.json",
        "results": base_dir / "results.json",
        "results_csv": base_dir / "results.csv",
        "uncertainty_components": base_dir / "uncertainty_components.csv",
        "summary": base_dir / "summary.txt",
        "metadata": base_dir / "run.metadata.txt",
        "artifact_checksums": base_dir / "artifact_checksums.txt",
    }
    for label, path in base_files.items():
        if not path.is_file():
            raise RuntimeError(f"missing base {label}: {path}")
    verify_checksums(base_files["artifact_checksums"])

    base_inputs = json.loads(base_files["inputs"].read_text(encoding="utf-8"))
    base_results = json.loads(base_files["results"].read_text(encoding="utf-8"))
    base_summary = read_kv(base_files["summary"])
    base_metadata = read_kv(base_files["metadata"])
    template_summary = read_kv(template_summary_path)

    expected_status = (
        "candidate_complete_included_sources_pending_fitter_stability_"
        "and_epsilon_m_uncertainty"
    )
    if base_results.get("status") != expected_status:
        raise RuntimeError("base candidate has unexpected status")
    if base_metadata.get("updated_fitter_stability") != (
        "pending_not_included_not_zero"
    ):
        raise RuntimeError("base candidate does not leave fitter stability pending")
    if template_summary.get("status") != "nominal_complete_user_confirmed":
        raise RuntimeError("template statistical source is not user confirmed")
    if template_summary.get("stat_definition") != (
        "absolute_shift_half_selfmix_vs_full_selfmix"
    ):
        raise RuntimeError("unexpected template statistical definition")

    template_stat_absolute_old_fraction = float(
        template_summary["stat_uncertainty"]
    )
    template_old_central = float(template_summary["f_dps"])
    template_relative = template_stat_absolute_old_fraction / template_old_central
    require_close(
        "template relative statistical fluctuation",
        100.0 * template_relative,
        float(template_summary["stat_uncertainty_percent"]),
    )

    revised = copy.deepcopy(base_results)
    assumptions = revised["assumptions"]
    assumptions["half_vs_full_selfmix_shift"] = (
        "formal_template_statistical_source"
    )
    assumptions["sigma_dps_statistical_sources"] = (
        "CR_DPS_yield_and_half_vs_full_selfmix_template_independent_quadrature"
    )
    assumptions["f_dps_statistical_sources"] = (
        "sigma_DPS_stat_and_total_cross_section_stat_independent_quadrature"
    )
    assumptions["statistical_definition_revision"] = (
        "2026-09-03_user_clarification"
    )

    f_result = revised["f_dps"]
    f_central = float(f_result["central"])
    f_cr_stat = float(f_result["stat_components"]["cr_dps_yield"])
    f_total_stat = float(f_result["stat_components"]["sigma_total"])
    f_template_stat = f_central * template_relative
    f_sigma_dps_stat = quadrature(f_cr_stat, f_template_stat)
    f_stat = quadrature(f_sigma_dps_stat, f_total_stat)
    f_result["stat"] = f_stat
    f_result["stat_percent"] = percent(f_stat, f_central)
    f_result["stat_components"] = {
        "cr_dps_yield": f_cr_stat,
        "cr_dps_yield_percent": percent(f_cr_stat, f_central),
        "half_vs_full_selfmix_template": f_template_stat,
        "half_vs_full_selfmix_template_percent": 100.0 * template_relative,
        "sigma_dps_combined": f_sigma_dps_stat,
        "sigma_dps_combined_percent": percent(f_sigma_dps_stat, f_central),
        "sigma_total": f_total_stat,
        "sigma_total_percent": percent(f_total_stat, f_central),
    }

    sigma_cut = revised["sigma_dps_mjj_ge_7p5_pb"]
    sigma_cut_central = float(sigma_cut["central"])
    cr_relative = float(base_inputs["control_region"]["relative_error"])
    sigma_cut_cr_stat = sigma_cut_central * cr_relative
    sigma_cut_template_stat = sigma_cut_central * template_relative
    sigma_cut_stat = quadrature(sigma_cut_cr_stat, sigma_cut_template_stat)
    sigma_cut["stat"] = sigma_cut_stat
    sigma_cut["stat_percent"] = percent(sigma_cut_stat, sigma_cut_central)
    sigma_cut["stat_components"] = {
        "cr_dps_yield": sigma_cut_cr_stat,
        "cr_dps_yield_percent": 100.0 * cr_relative,
        "half_vs_full_selfmix_template": sigma_cut_template_stat,
        "half_vs_full_selfmix_template_percent": 100.0 * template_relative,
    }

    sigma_no_mass = revised["sigma_dps_no_mjj_cut_pb"]
    sigma_no_mass_central = float(sigma_no_mass["central"])
    sigma_no_mass_cr_stat = sigma_no_mass_central * cr_relative
    sigma_no_mass_template_stat = sigma_no_mass_central * template_relative
    sigma_no_mass_stat = quadrature(
        sigma_no_mass_cr_stat, sigma_no_mass_template_stat
    )
    sigma_no_mass["stat"] = sigma_no_mass_stat
    sigma_no_mass["stat_percent"] = percent(
        sigma_no_mass_stat, sigma_no_mass_central
    )
    sigma_no_mass["stat_components"] = {
        "cr_dps_yield": sigma_no_mass_cr_stat,
        "cr_dps_yield_percent": 100.0 * cr_relative,
        "half_vs_full_selfmix_template": sigma_no_mass_template_stat,
        "half_vs_full_selfmix_template_percent": 100.0 * template_relative,
    }

    sigma_eff = revised["sigma_eff_mb"]
    sigma_eff_central = float(sigma_eff["central"])
    sigma_dps_relative = sigma_cut_stat / sigma_cut_central
    sigma_eff_from_dps = sigma_eff_central * sigma_dps_relative
    single = base_inputs["single_jpsi_cross_section"]
    single_relative = (
        2.0 * float(single["visible_stat_nb"]) / float(single["visible_central_nb"])
    )
    sigma_eff_from_single = sigma_eff_central * single_relative
    sigma_eff_stat = quadrature(sigma_eff_from_dps, sigma_eff_from_single)
    sigma_eff["stat"] = sigma_eff_stat
    sigma_eff["stat_percent"] = percent(sigma_eff_stat, sigma_eff_central)
    sigma_eff["stat_components"] = {
        "sigma_dps": sigma_eff_from_dps,
        "sigma_dps_percent": 100.0 * sigma_dps_relative,
        "atlas_single_jpsi": sigma_eff_from_single,
        "atlas_single_jpsi_percent": 100.0 * single_relative,
    }

    # Central values and every systematic field must remain byte-for-byte
    # numerically identical to the base candidate.
    for quantity in (
        "f_dps",
        "sigma_dps_mjj_ge_7p5_pb",
        "sigma_dps_no_mjj_cut_pb",
        "sigma_eff_mb",
    ):
        require_close(
            f"{quantity} central",
            float(revised[quantity]["central"]),
            float(base_results[quantity]["central"]),
        )
        for key in (
            "included_syst_up",
            "included_syst_down",
            "combined_syst_up",
            "combined_syst_down",
        ):
            if key in base_results[quantity]:
                require_close(
                    f"{quantity} {key}",
                    float(revised[quantity][key]),
                    float(base_results[quantity][key]),
                )

    with base_files["results_csv"].open(newline="", encoding="utf-8") as handle:
        result_rows = list(csv.DictReader(handle))
    with base_files["results_csv"].open(newline="", encoding="utf-8") as handle:
        result_fields = list(csv.DictReader(handle).fieldnames or [])
    revised_stats = {
        ("f_DPS", "mJJ>=7.5"): f_stat,
        ("sigma_DPS", "mJJ>=7.5"): sigma_cut_stat,
        ("sigma_DPS", "no_mJJ_cut"): sigma_no_mass_stat,
        ("sigma_eff", "single_Jpsi_matched_no_pair_mass_cut"): sigma_eff_stat,
    }
    for row in result_rows:
        key = (row["quantity"], row["phase_space"])
        stat = revised_stats[key]
        row["stat"] = stat
        row["stat_percent"] = 100.0 * stat / float(row["central"])

    with base_files["uncertainty_components"].open(
        newline="", encoding="utf-8"
    ) as handle:
        component_reader = csv.DictReader(handle)
        component_fields = list(component_reader.fieldnames or [])
        base_components = list(component_reader)
    system_components = [row for row in base_components if row["type"] != "stat"]
    f_system = [row for row in system_components if row["quantity"] == "f_DPS"]
    cut_system = [
        row
        for row in system_components
        if row["quantity"] == "sigma_DPS_mJJ_ge_7p5"
    ]
    eff_system = [
        row for row in system_components if row["quantity"] == "sigma_eff"
    ]
    atlas_stat = next(
        row
        for row in base_components
        if row["quantity"] == "sigma_eff"
        and row["source"] == "ATLAS_single_Jpsi"
        and row["type"] == "stat"
    )
    atlas_stat["effect_up"] = sigma_eff_from_single
    atlas_stat["effect_down"] = sigma_eff_from_single
    atlas_stat["effect_up_percent"] = percent(
        sigma_eff_from_single, sigma_eff_central
    )
    atlas_stat["effect_down_percent"] = percent(
        sigma_eff_from_single, sigma_eff_central
    )
    component_rows = [
        effect_row(
            "f_DPS", "CR_DPS_yield", "stat", f_cr_stat, "fraction",
            "included_in_sigma_DPS_stat", f_central,
        ),
        effect_row(
            "f_DPS", "half_vs_full_selfmix_template", "stat",
            f_template_stat, "fraction", "included_in_sigma_DPS_stat",
            f_central,
        ),
        effect_row(
            "f_DPS", "sigma_DPS_combined", "stat_subtotal",
            f_sigma_dps_stat, "fraction",
            "combined_independent_of_sigma_total_stat", f_central,
        ),
        effect_row(
            "f_DPS", "sigma_total", "stat", f_total_stat, "fraction",
            "included_independent_of_sigma_DPS_stat", f_central,
        ),
        *f_system,
        effect_row(
            "sigma_DPS_mJJ_ge_7p5", "CR_DPS_yield", "stat",
            sigma_cut_cr_stat, "pb", "included", sigma_cut_central,
        ),
        effect_row(
            "sigma_DPS_mJJ_ge_7p5", "half_vs_full_selfmix_template", "stat",
            sigma_cut_template_stat, "pb", "included", sigma_cut_central,
        ),
        effect_row(
            "sigma_DPS_mJJ_ge_7p5", "combined", "stat_subtotal",
            sigma_cut_stat, "pb", "independent_quadrature", sigma_cut_central,
        ),
        *cut_system,
        effect_row(
            "sigma_eff", "sigma_DPS", "stat", sigma_eff_from_dps, "mb",
            "included", sigma_eff_central,
        ),
        atlas_stat,
        *eff_system,
    ]

    source_paths = {
        **base_files,
        "template_stat_summary": template_summary_path,
        "revision_script": Path(__file__).resolve(),
    }
    source_hashes = {name: sha256(path) for name, path in source_paths.items()}
    revised_inputs = {
        "status": "frozen_snapshot_for_statistical_definition_revision",
        "base_candidate": {
            "path": str(base_dir.resolve()),
            "artifact_checksums_verified": True,
        },
        "template_statistical_source": {
            "path": str(template_summary_path.resolve()),
            "absolute_shift_in_old_fraction": template_stat_absolute_old_fraction,
            "old_fraction_central": template_old_central,
            "relative": template_relative,
            "percent": 100.0 * template_relative,
        },
        "inherited_inputs": base_inputs,
        "source_paths": {name: str(path.resolve()) for name, path in source_paths.items()},
        "source_sha256": source_hashes,
    }

    generated_at = datetime.now().astimezone().isoformat(timespec="seconds")
    summary_lines = [
        f"status={revised['status']}",
        "f_dps_status=nominal_updated_user_clarified_statistical_scheme",
        f"f_dps={f_central:.12g}",
        f"f_dps_stat={f_stat:.12g}",
        f"f_dps_stat_percent={percent(f_stat, f_central):.12g}",
        f"f_dps_stat_from_cr_dps_yield={f_cr_stat:.12g}",
        f"f_dps_stat_from_cr_dps_yield_percent={percent(f_cr_stat, f_central):.12g}",
        f"f_dps_stat_from_template={f_template_stat:.12g}",
        f"f_dps_stat_from_template_percent={100.0 * template_relative:.12g}",
        f"f_dps_stat_from_sigma_dps_combined={f_sigma_dps_stat:.12g}",
        f"f_dps_stat_from_sigma_dps_combined_percent={percent(f_sigma_dps_stat, f_central):.12g}",
        f"f_dps_stat_from_sigma_total={f_total_stat:.12g}",
        f"f_dps_stat_from_sigma_total_percent={percent(f_total_stat, f_central):.12g}",
        f"f_dps_total_syst_up={f_result['combined_syst_up']:.12g}",
        f"f_dps_total_syst_down={f_result['combined_syst_down']:.12g}",
        f"f_dps_total_syst_up_percent={f_result['combined_syst_up_percent']:.12g}",
        f"f_dps_total_syst_down_percent={f_result['combined_syst_down_percent']:.12g}",
        "f_dps_sigma_total_systematics=not_included_by_user_definition",
        f"sigma_dps_mjj_ge_7p5_pb={sigma_cut_central:.12g}",
        f"sigma_dps_mjj_ge_7p5_stat_pb={sigma_cut_stat:.12g}",
        f"sigma_dps_mjj_ge_7p5_stat_percent={percent(sigma_cut_stat, sigma_cut_central):.12g}",
        f"sigma_dps_stat_from_cr_yield_pb={sigma_cut_cr_stat:.12g}",
        f"sigma_dps_stat_from_cr_yield_percent={100.0 * cr_relative:.12g}",
        f"sigma_dps_stat_from_template_pb={sigma_cut_template_stat:.12g}",
        f"sigma_dps_stat_from_template_percent={100.0 * template_relative:.12g}",
        f"sigma_dps_mjj_ge_7p5_included_syst_up_pb={sigma_cut['included_syst_up']:.12g}",
        f"sigma_dps_mjj_ge_7p5_included_syst_down_pb={sigma_cut['included_syst_down']:.12g}",
        f"sigma_dps_mjj_ge_7p5_included_syst_up_percent={sigma_cut['included_syst_up_percent']:.12g}",
        f"sigma_dps_mjj_ge_7p5_included_syst_down_percent={sigma_cut['included_syst_down_percent']:.12g}",
        f"sigma_dps_total_systematic_status={sigma_cut['total_syst_status']}",
        f"sigma_dps_no_mjj_cut_pb={sigma_no_mass_central:.12g}",
        f"sigma_dps_no_mjj_cut_stat_pb={sigma_no_mass_stat:.12g}",
        f"sigma_dps_no_mjj_cut_included_syst_up_pb={sigma_no_mass['included_syst_up']:.12g}",
        f"sigma_dps_no_mjj_cut_included_syst_down_pb={sigma_no_mass['included_syst_down']:.12g}",
        f"sigma_eff_mb={sigma_eff_central:.12g}",
        f"sigma_eff_stat_mb={sigma_eff_stat:.12g}",
        f"sigma_eff_stat_percent={percent(sigma_eff_stat, sigma_eff_central):.12g}",
        f"sigma_eff_stat_from_sigma_dps_mb={sigma_eff_from_dps:.12g}",
        f"sigma_eff_stat_from_sigma_dps_percent={100.0 * sigma_dps_relative:.12g}",
        f"sigma_eff_stat_from_atlas_single_jpsi_mb={sigma_eff_from_single:.12g}",
        f"sigma_eff_stat_from_atlas_single_jpsi_percent={100.0 * single_relative:.12g}",
        f"sigma_eff_included_syst_up_mb={sigma_eff['included_syst_up']:.12g}",
        f"sigma_eff_included_syst_down_mb={sigma_eff['included_syst_down']:.12g}",
        f"sigma_eff_included_syst_up_percent={sigma_eff['included_syst_up_percent']:.12g}",
        f"sigma_eff_included_syst_down_percent={sigma_eff['included_syst_down_percent']:.12g}",
        f"sigma_eff_total_systematic_status={sigma_eff['total_syst_status']}",
        "pair_total_fitter_stability_status=pending_not_included_not_zero",
        "epsilon_m_uncertainty_status=pending_not_included_not_zero",
        "new_roofit_fits_run=0",
    ]
    metadata_lines = [
        "status=complete_candidate_pending_fitter_stability_and_epsilon_m_uncertainty",
        f"output_tag={output_dir.name}",
        f"generated_at={generated_at}",
        f"root_version={root_version}",
        f"git_head={git_head}",
        f"script_sha256={source_hashes['revision_script']}",
        "record_type=deterministic_statistical_revision_no_ROOT_fits",
        "base_candidate_artifacts_verified=true",
        "sigma_dps_stat_definition=CR_DPS_yield_and_half_vs_full_selfmix_template_independent_quadrature",
        "f_dps_stat_definition=sigma_DPS_stat_and_sigma_total_stat_independent_quadrature",
        "systematic_values=unchanged_from_base_candidate",
        "updated_fitter_stability=pending_not_included_not_zero",
        "epsilon_m_uncertainty=pending_not_included_not_zero",
        "user_gate=2026-09-03_statistical_definition_clarification",
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
            staging / "uncertainty_components.csv",
            component_fields,
            component_rows,
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
        (staging / "complete.marker").write_text("complete\n", encoding="utf-8")
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
    print(f"f_DPS={f_central:.12g} +/- {f_stat:.12g} stat")
    print(
        f"sigma_DPS(mJJ>=7.5)={sigma_cut_central:.12g} "
        f"+/- {sigma_cut_stat:.12g} stat pb"
    )
    print(f"sigma_eff={sigma_eff_central:.12g} +/- {sigma_eff_stat:.12g} stat mb")


if __name__ == "__main__":
    main()
