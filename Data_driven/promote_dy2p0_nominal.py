#!/usr/bin/env python3
"""Consolidate the accepted dy2p0 CR variation into a new nominal tag."""

from __future__ import annotations

import csv
import hashlib
import math
import re
import subprocess
import sys
from pathlib import Path


REPO = Path("/eos/home-l/leyao/26JJ/JPsiJPsi")
HERE = REPO / "Data_driven"
RESULTS = HERE / "results"
SOURCE_TAG = (
    "newdata167_accmix23_effmix19_fixedhalf_seed20260902_"
    "combcomb_sigcombfix0_root640_v2_20260902"
)
DY2P0_TAG = "cr_dy2p0_fixedhalf_newdata167_root640_v1_20260902"


def read_kv(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        if key in result:
            raise RuntimeError(f"duplicate key {key} in {path}")
        result[key] = value
    return result


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_kv(path: Path, rows: list[tuple[str, object]]) -> None:
    path.write_text(
        "".join(f"{key}={value}\n" for key, value in rows), encoding="utf-8"
    )


def validate_checksums(directory: Path) -> None:
    subprocess.run(
        ["sha256sum", "-c", "artifact_checksums.txt"],
        cwd=directory,
        check=True,
        stdout=subprocess.DEVNULL,
    )


def read_fit_qa(paths: list[Path]) -> dict[str, object]:
    rows: list[dict[str, str]] = []
    for path in paths:
        with path.open(newline="", encoding="utf-8") as handle:
            rows.extend(csv.DictReader(handle))
    if not rows:
        raise RuntimeError("no fit QA rows")
    accepted = sum(int(row["accepted"]) for row in rows)
    return {
        "fits": len(rows),
        "accepted": accepted,
        "failed": len(rows) - accepted,
        "min_covQual": min(int(row["covQual"]) for row in rows),
        "max_edm": max(float(row["edm"]) for row in rows),
        "comb_fixed": sum(int(row["comb_comb_fixed_zero"]) for row in rows),
        "sig_fixed": sum(int(row["sig_comb_fixed_zero"]) for row in rows),
    }


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit(f"usage: {Path(sys.argv[0]).name} UNIQUE_OUTPUT_TAG")
    output_tag = sys.argv[1]
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", output_tag):
        raise RuntimeError("invalid output tag")

    source = RESULTS / SOURCE_TAG
    dy2p0 = RESULTS / DY2P0_TAG
    output = RESULTS / output_tag
    if output.exists():
        raise RuntimeError(f"refusing to overwrite existing tag: {output}")
    for directory in (source, dy2p0):
        validate_checksums(directory)

    nominal = read_kv(source / "summary.txt")
    nominal_meta = read_kv(source / "run.metadata.txt")
    dy = read_kv(dy2p0 / "summary.txt")
    dy_meta = read_kv(dy2p0 / "run.metadata.txt")
    if nominal.get("status") != "nominal_complete_user_confirmed":
        raise RuntimeError("source tag is not a confirmed nominal")
    if dy.get("status") != "candidate_complete_pending_user_confirmation":
        raise RuntimeError("dy2p0 tag is not the reviewed candidate")
    if dy.get("definition") != "abs_delta_y_ge_2p0_and_abs_delta_phi_lt_pi_over_2":
        raise RuntimeError("unexpected dy2p0 definition")
    if dy_meta.get("nominal_tag") != SOURCE_TAG:
        raise RuntimeError("dy2p0 source nominal mismatch")
    root_version = subprocess.check_output(["root-config", "--version"], text=True).strip()
    if root_version != "6.40.02":
        raise RuntimeError(f"ROOT version mismatch: {root_version}")

    f_dps = float(nominal["f_dps"])
    stat = float(nominal["stat_uncertainty"])
    usage = float(nominal["dps_mc_usage_systematic"])
    phi70 = float(nominal["cr_phi70_signed_shift"])
    phi110 = float(nominal["cr_phi110_signed_shift"])
    dy_shift = float(dy["signed_shift"])
    if not math.isclose(float(dy["nominal_f_dps"]), f_dps, abs_tol=5e-13):
        raise RuntimeError("dy2p0 central value mismatch")
    cr_up = max(0.0, phi70, phi110, dy_shift)
    cr_down = max(0.0, -phi70, -phi110, -dy_shift)
    total_up = math.hypot(cr_up, usage)
    total_down = math.hypot(cr_down, usage)

    qa = read_fit_qa(
        [
            source / "pp_data_2d/fit_results.csv",
            source / "cr/phi70_pp/fit_results.csv",
            source / "cr/phi110_pp/fit_results.csv",
            dy2p0 / "pp_data_dy2p0/fit_results.csv",
        ]
    )
    if qa["failed"] != 0 or qa["min_covQual"] < 3:
        raise RuntimeError(f"fit QA failed: {qa}")

    output.mkdir(parents=True)
    summary_rows: list[tuple[str, object]] = [
        ("status", "nominal_complete_user_confirmed"),
        ("artifact", "newdata_fixed_half_dpsmc_nominal_fdps_dy2p0_cr"),
        ("source_nominal_tag", SOURCE_TAG),
        ("dy2p0_source_tag", DY2P0_TAG),
        ("user_confirmation", "dy2p0_replaces_dy2p4_accepted_2026-09-02"),
        ("method", nominal["method"]),
        ("f_dps", nominal["f_dps"]),
        ("f_sps", nominal["f_sps"]),
        ("stat_definition", nominal["stat_definition"]),
        ("stat_uncertainty", nominal["stat_uncertainty"]),
        ("stat_uncertainty_percent", nominal["stat_uncertainty_percent"]),
        ("stat_signed_half_minus_full", nominal["stat_signed_half_minus_full"]),
        ("stat_full_reference_f_dps", nominal["stat_full_reference_f_dps"]),
        ("stat_seed", nominal["stat_seed"]),
        ("stat_sampling", nominal["stat_sampling"]),
        ("stat_input_entries", nominal["stat_input_entries"]),
        ("stat_selected_entries", nominal["stat_selected_entries"]),
        ("cr_systematic_definition", "zero_leakage_phi70_phi110_dy2p0_envelope_same_fixed_half"),
        ("cr_systematic_up", f"{cr_up:.12f}"),
        ("cr_systematic_down", f"{cr_down:.12f}"),
        ("cr_systematic_symmetric_max", f"{max(cr_up, cr_down):.12f}"),
        ("cr_systematic_up_percent", f"{100.0 * cr_up / f_dps:.12f}"),
        ("cr_systematic_down_percent", f"{100.0 * cr_down / f_dps:.12f}"),
        ("cr_phi70_f_dps", nominal["cr_phi70_f_dps"]),
        ("cr_phi70_signed_shift", nominal["cr_phi70_signed_shift"]),
        ("cr_phi70_negative_sps_bins", nominal["cr_phi70_negative_sps_bins"]),
        ("cr_phi110_f_dps", nominal["cr_phi110_f_dps"]),
        ("cr_phi110_signed_shift", nominal["cr_phi110_signed_shift"]),
        ("cr_phi110_negative_sps_bins", nominal["cr_phi110_negative_sps_bins"]),
        ("cr_dy2p0_f_dps", dy["dy2p0_f_dps"]),
        ("cr_dy2p0_signed_shift", dy["signed_shift"]),
        ("cr_dy2p0_negative_sps_bins", dy["negative_sps_bins"]),
        ("cr_dy2p0_data_control", dy["data_control"]),
        ("cr_dy2p0_data_control_error", dy["data_control_error"]),
        ("cr_dy2p0_dps_mc_control", dy["dps_mc_control"]),
        ("cr_dy2p0_alpha", dy["alpha"]),
        ("cr_dy2p4_role", "superseded_CR_crosscheck_not_in_nominal_systematic"),
        ("cr_dy2p4_f_dps", nominal["cr_dy2p4_f_dps"]),
        ("cr_dy2p4_signed_shift", nominal["cr_dy2p4_signed_shift"]),
        ("dps_mc_usage_systematic_definition", nominal["dps_mc_usage_systematic_definition"]),
        ("dps_mc_self_mixed_f_dps", nominal["dps_mc_self_mixed_f_dps"]),
        ("dps_mc_usage_signed_shift", nominal["dps_mc_usage_signed_shift"]),
        ("dps_mc_usage_systematic", nominal["dps_mc_usage_systematic"]),
        ("dps_mc_usage_systematic_percent", nominal["dps_mc_usage_systematic_percent"]),
        ("total_systematic_definition", "quadrature_CR_and_DPS_MC_event_mixing_usage"),
        ("total_systematic_up", f"{total_up:.12f}"),
        ("total_systematic_down", f"{total_down:.12f}"),
        ("total_systematic_up_percent", f"{100.0 * total_up / f_dps:.12f}"),
        ("total_systematic_down_percent", f"{100.0 * total_down / f_dps:.12f}"),
        ("data_event_mixing_role", nominal["data_event_mixing_role"]),
        ("event_mixing_f_dps", nominal["event_mixing_f_dps"]),
        ("event_mixing_signed_shift_from_nominal", nominal["event_mixing_signed_shift_from_nominal"]),
        ("sps_leakage_fraction", "0_for_central_and_all_CR_variations"),
        ("weightdata_selected_entries", nominal["weightdata_selected_entries"]),
        ("weightdata_fit_range_entries", nominal["weightdata_fit_range_entries"]),
        ("dps_reco_selected_entries", nominal["dps_reco_selected_entries"]),
        ("pp_fits_total", qa["fits"]),
        ("pp_fit_failed", qa["failed"]),
        ("pp_fit_min_covQual", qa["min_covQual"]),
        ("pp_fit_max_edm", f"{qa['max_edm']:.15g}"),
        ("pp_fit_comb_comb_fixed_zero", qa["comb_fixed"]),
        ("pp_fit_sig_comb_fixed_zero", qa["sig_fixed"]),
        ("pp_fit_error_convention", nominal["pp_fit_error_convention"]),
        ("pp_comb_comb_boundary_strategy", nominal["pp_comb_comb_boundary_strategy"]),
        ("pp_sig_comb_boundary_strategy", nominal["pp_sig_comb_boundary_strategy"]),
        ("sigma_dps_status", "requires_repropagation_with_dy2p0_CR_systematic"),
        ("sigma_eff_status", "requires_repropagation_with_dy2p0_CR_systematic"),
        ("bootstrap_status", nominal["bootstrap_status"]),
        ("category_systematic_status", nominal["category_systematic_status"]),
        ("pure_dps_closure_status", nominal["pure_dps_closure_status"]),
    ]
    write_kv(output / "summary.txt", summary_rows)

    with (source / "cr_variations.csv").open(newline="", encoding="utf-8") as handle:
        old_rows = list(csv.DictReader(handle))
        fieldnames = handle.seek(0) or next(csv.reader(handle))
    selected = [row for row in old_rows if row["variation"] in {"nominal", "phi70", "phi110"}]
    selected.append(
        {
            "variation": "dy2p0",
            "dy_min": "2.0",
            "phi_max": "1.5707963267948966",
            "cells": dy["pp_fits"],
            "f_dps": dy["dy2p0_f_dps"],
            "signed_shift": dy["signed_shift"],
            "negative_sps_bins": dy["negative_sps_bins"],
            "pp_accepted": dy["pp_accepted_fits"],
        }
    )
    with (output / "cr_variations.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(selected)

    metadata_rows: list[tuple[str, object]] = [
        ("status", "nominal_complete_user_confirmed"),
        ("artifact", "newdata_fixed_half_dpsmc_nominal_fdps_dy2p0_cr_consolidation"),
        ("output_tag", output_tag),
        ("source_nominal_tag", SOURCE_TAG),
        ("source_nominal_summary_sha256", sha256(source / "summary.txt")),
        ("source_nominal_metadata_sha256", sha256(source / "run.metadata.txt")),
        ("dy2p0_source_tag", DY2P0_TAG),
        ("dy2p0_summary_sha256", sha256(dy2p0 / "summary.txt")),
        ("dy2p0_metadata_sha256", sha256(dy2p0 / "run.metadata.txt")),
        ("user_confirmation", "dy2p0_replaces_dy2p4_accepted_2026-09-02"),
        ("root_version", root_version),
        ("git_head", subprocess.check_output(["git", "-C", str(REPO), "rev-parse", "HEAD"], text=True).strip()),
        ("model_sha256", nominal_meta["model_sha256"]),
        ("weightdata_sha256", nominal_meta["weightdata_sha256"]),
        ("half_dps_mc_sha256", dy_meta["half_dps_mc_sha256"]),
        ("acceptance_sha256", nominal_meta["acceptance_sha256"]),
        ("efficiency_sha256", nominal_meta["efficiency_sha256"]),
        ("seed_policy", nominal_meta["seed_policy"]),
        ("fit_policy", nominal_meta["fit_policy"]),
        ("cr_variations", "phi70_phi110_dy2p0"),
        ("dy2p4_role", "superseded_CR_crosscheck_not_in_nominal_systematic"),
        ("cr_systematic", f"+{cr_up:.12f}/-{cr_down:.12f}"),
        ("dps_mc_usage_systematic", f"{usage:.12f}"),
        ("total_systematic", f"+{total_up:.12f}/-{total_down:.12f}"),
        ("pp_fits_referenced", qa["fits"]),
        ("pp_fit_failed", qa["failed"]),
        ("pp_fit_min_covQual", qa["min_covQual"]),
        ("pp_fit_max_edm", f"{qa['max_edm']:.15g}"),
        ("pp_fit_comb_comb_fixed_zero", qa["comb_fixed"]),
        ("pp_fit_sig_comb_fixed_zero", qa["sig_fixed"]),
        ("new_roofit_fits_run", 0),
        ("real_newlines", "true"),
    ]
    write_kv(output / "run.metadata.txt", metadata_rows)
    (output / "complete.marker").write_text(
        "nominal_complete_user_confirmed\n", encoding="utf-8"
    )
    artifact_paths = [
        output / "summary.txt",
        output / "cr_variations.csv",
        output / "run.metadata.txt",
        output / "complete.marker",
    ]
    (output / "artifact_checksums.txt").write_text(
        "".join(f"{sha256(path)}  {path.name}\n" for path in artifact_paths),
        encoding="utf-8",
    )
    validate_checksums(output)
    print(f"DY2P0_NOMINAL_PROMOTED tag={output_tag}")


if __name__ == "__main__":
    main()
