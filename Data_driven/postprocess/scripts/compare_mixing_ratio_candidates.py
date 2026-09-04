#!/usr/bin/env python3
"""Compare two complete f_DPS and sigma_eff candidate generations."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path


def read_kv(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key] = value
    return values


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def compare(old: float, new: float) -> tuple[float, float]:
    delta = new - old
    return delta, 100.0 * delta / old


def main() -> None:
    if len(sys.argv) != 6:
        raise SystemExit(
            "usage: compare_mixing_ratio_candidates.py "
            "MIX0P8_F MIX0P8_SIGMA MIX0P85_F MIX0P85_SIGMA OUTPUT_DIR"
        )
    old_f_path, old_sigma_path, new_f_path, new_sigma_path, output = map(
        Path, sys.argv[1:]
    )
    if output.exists():
        raise RuntimeError(f"refusing to overwrite {output}")
    for path in (old_f_path, old_sigma_path, new_f_path, new_sigma_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"missing input {path}")

    old_f = read_kv(old_f_path)
    old_sigma = read_kv(old_sigma_path)
    new_f = read_kv(new_f_path)
    new_sigma = read_kv(new_sigma_path)
    for label, values in (
        ("mix0p8_f", old_f),
        ("mix0p85_f", new_f),
    ):
        if values.get("status") != "candidate_complete_pending_user_confirmation":
            raise RuntimeError(f"unexpected {label} status")
    if float(old_f["f_dps"]) != float(old_sigma["f_dps"]):
        raise RuntimeError("mix0p8 f_DPS mismatch between f and sigma summaries")
    if float(new_f["f_dps"]) != float(new_sigma["f_dps"]):
        raise RuntimeError("mix0p85 f_DPS mismatch between f and sigma summaries")

    definitions = [
        ("f_DPS", "central", old_f, new_f, "f_dps", "fraction"),
        ("f_DPS", "stat", old_f, new_f, "stat_uncertainty", "fraction"),
        ("f_DPS", "CR_syst_up", old_f, new_f, "cr_systematic_up", "fraction"),
        ("f_DPS", "CR_syst_down", old_f, new_f, "cr_systematic_down", "fraction"),
        ("f_DPS", "DPS_MC_usage_syst", old_f, new_f,
         "dps_mc_usage_systematic", "fraction"),
        ("f_DPS", "total_syst_up", old_f, new_f,
         "total_systematic_up", "fraction"),
        ("f_DPS", "total_syst_down", old_f, new_f,
         "total_systematic_down", "fraction"),
        ("pair_total", "central", old_sigma, new_sigma, "pair_total_pb", "pb"),
        ("sigma_DPS_mJJ_ge_7p5", "central", old_sigma, new_sigma,
         "sigma_dps_mjj_ge_7p5_pb", "pb"),
        ("sigma_DPS_mJJ_ge_7p5", "stat", old_sigma, new_sigma,
         "sigma_dps_mjj_ge_7p5_stat_pb", "pb"),
        ("sigma_DPS_mJJ_ge_7p5", "syst_up", old_sigma, new_sigma,
         "sigma_dps_mjj_ge_7p5_syst_up_pb", "pb"),
        ("sigma_DPS_mJJ_ge_7p5", "syst_down", old_sigma, new_sigma,
         "sigma_dps_mjj_ge_7p5_syst_down_pb", "pb"),
        ("sigma_DPS_no_mJJ_cut", "central", old_sigma, new_sigma,
         "sigma_dps_no_mjj_cut_pb", "pb"),
        ("sigma_eff", "central", old_sigma, new_sigma, "sigma_eff_mb", "mb"),
        ("sigma_eff", "stat", old_sigma, new_sigma, "sigma_eff_stat_mb", "mb"),
        ("sigma_eff", "syst_up", old_sigma, new_sigma,
         "sigma_eff_syst_up_mb", "mb"),
        ("sigma_eff", "syst_down", old_sigma, new_sigma,
         "sigma_eff_syst_down_mb", "mb"),
    ]
    rows = []
    for quantity, component, old_values, new_values, key, unit in definitions:
        old = float(old_values[key])
        new = float(new_values[key])
        delta, relative = compare(old, new)
        rows.append({
            "quantity": quantity,
            "component": component,
            "mix0p8_dps0p2": old,
            "mix0p85_dps0p15": new,
            "new_minus_old": delta,
            "relative_to_mix0p8_percent": relative,
            "unit": unit,
        })

    output.mkdir(parents=True)
    comparison_path = output / "comparison.csv"
    with comparison_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    sources = {
        "mix0p8_f_summary": str(old_f_path.resolve()),
        "mix0p8_sigma_summary": str(old_sigma_path.resolve()),
        "mix0p85_f_summary": str(new_f_path.resolve()),
        "mix0p85_sigma_summary": str(new_sigma_path.resolve()),
    }
    inputs = {
        "status": "complete_frozen_comparison_inputs",
        "sources": sources,
        "sha256": {
            label: sha256(Path(path)) for label, path in sources.items()
        },
        "comparison_orientation": "mix0p85_minus_mix0p8",
        "shared_dps_inventory": "deduplicated_60_files_excluding_11_to_15",
        "shared_sampling": "floor_3835_over_2_equals_1917_seed_20260902",
    }
    (output / "inputs.json").write_text(
        json.dumps(inputs, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    def row(quantity: str, component: str) -> dict[str, object]:
        return next(
            item for item in rows
            if item["quantity"] == quantity and item["component"] == component
        )

    f_row = row("f_DPS", "central")
    dps_row = row("sigma_DPS_mJJ_ge_7p5", "central")
    eff_row = row("sigma_eff", "central")
    summary = [
        "status=candidate_comparison_complete_pending_user_confirmation",
        "orientation=mix0p85_dps0p15_minus_mix0p8_dps0p2",
        "shared_dps_files=60_excluding_11_12_13_14_15",
        "shared_dps_selected_entries=3835",
        "shared_nominal_selected_entries=1917",
        "shared_seed=20260902",
        f"mix0p8_f_dps={float(old_f['f_dps']):.12g}",
        f"mix0p85_f_dps={float(new_f['f_dps']):.12g}",
        f"f_dps_delta={float(f_row['new_minus_old']):.12g}",
        f"f_dps_relative_change_percent={float(f_row['relative_to_mix0p8_percent']):.12g}",
        f"sigma_dps_mjj_ge_7p5_delta_pb={float(dps_row['new_minus_old']):.12g}",
        f"sigma_dps_mjj_ge_7p5_relative_change_percent={float(dps_row['relative_to_mix0p8_percent']):.12g}",
        f"sigma_eff_delta_mb={float(eff_row['new_minus_old']):.12g}",
        f"sigma_eff_relative_change_percent={float(eff_row['relative_to_mix0p8_percent']):.12g}",
        "mix0p8_pair_status=accepted_nominal_user_confirmed",
        "mix0p85_pair_status=temporary_total_only_pending_fitter_stability",
        "epsilon_m_uncertainty_status=pending_not_included_not_zero_for_both",
    ]
    (output / "summary.txt").write_text(
        "\n".join(summary) + "\n", encoding="utf-8"
    )
    repo = Path("/eos/home-l/leyao/26JJ/JPsiJPsi")
    metadata = [
        "status=candidate_comparison_complete_pending_user_confirmation",
        "artifact=sps0p8_dps0p2_vs_sps0p85_dps0p15_dedup60_comparison",
        f"git_head={subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'HEAD'], text=True).strip()}",
        "new_roofit_fits_run=0",
        f"script_sha256={sha256(Path(__file__).resolve())}",
    ]
    (output / "run.metadata.txt").write_text(
        "\n".join(metadata) + "\n", encoding="utf-8"
    )
    artifact_paths = sorted(
        path for path in output.iterdir()
        if path.name not in {"artifact_checksums.txt", "complete.marker"}
    )
    with (output / "artifact_checksums.txt").open("w", encoding="utf-8") as handle:
        for path in artifact_paths:
            handle.write(f"{sha256(path)}  {path.resolve()}\n")
    (output / "complete.marker").write_text(
        "candidate_comparison_complete_pending_user_confirmation\n",
        encoding="utf-8",
    )
    print(
        "MIXING_RATIO_COMPARISON_COMPLETE "
        f"f_DPS_delta={float(f_row['new_minus_old']):.12g} "
        f"sigma_eff_delta_mb={float(eff_row['new_minus_old']):.12g}"
    )


if __name__ == "__main__":
    main()
