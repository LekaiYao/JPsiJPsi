#!/usr/bin/env python3
"""Replace the correction systematic in a confirmed data-driven artifact.

The central values, statistical uncertainties, data-driven method systematic,
fitter systematic, and epsilon_m treatment are inherited unchanged.  Only the
pair acceptance/efficiency source and quantities that depend on it are
recomputed.  No ROOT/RooFit fit is run.
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
EXPECTED_RULE = (
    "For SPS and DPS separately, preserve the signed closure residuals and "
    "calculate (1 + acceptance_residual) * (1 + efficiency_residual) - 1; "
    "sys3 is the larger absolute combined residual."
)


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
    return central / (1.0 - relative_down) - central


def inverse_down(central: float, relative_up: float) -> float:
    return central - central / (1.0 + relative_up)


def close(actual: float, expected: float) -> bool:
    return math.isclose(actual, expected, rel_tol=5e-11, abs_tol=1e-12)


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
            raise RuntimeError(f"checksum mismatch for {path}: {actual} != {expected}")


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_csv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def find_row(
    rows: list[dict[str, str]], quantity: str, source: str, kind: str | None = None
) -> dict[str, str]:
    matches = [
        row
        for row in rows
        if row["quantity"] == quantity
        and row["source"] == source
        and (kind is None or row["type"] == kind)
    ]
    if len(matches) != 1:
        raise RuntimeError(f"expected one component {quantity}/{source}/{kind}")
    return matches[0]


def main() -> None:
    if len(sys.argv) != 5:
        raise SystemExit(
            "usage: update_confirmed_nominal_correction.py "
            "BASE_CONFIRMED_DIR RECOMBINATION_CSV CORRECTION_DIR OUTPUT_DIR"
        )
    base_dir, recombination_path, correction_dir, output_dir = map(Path, sys.argv[1:])
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
    correction_files = {
        "metadata": correction_dir / "metadata.json",
        "sys3": correction_dir / "sys3.csv",
        "checksums": correction_dir / "checksums.sha256",
    }
    sources = {
        **{f"base_{key}": value for key, value in base_files.items()},
        **{f"correction_{key}": value for key, value in correction_files.items()},
        "recombination": recombination_path,
        "revision_script": Path(__file__).resolve(),
    }
    for label, path in sources.items():
        if not path.is_file():
            raise FileNotFoundError(f"missing {label}: {path}")
    verify_checksums(base_files["artifact_checksums"])
    verify_checksums(correction_files["checksums"])

    root_version = subprocess.check_output(["root-config", "--version"], text=True).strip()
    if root_version != EXPECTED_ROOT_VERSION:
        raise RuntimeError(f"ROOT version mismatch: {root_version}")
    repo = Path(__file__).resolve().parents[3]
    git_head = subprocess.check_output(
        ["git", "-C", str(repo), "rev-parse", "HEAD"], text=True
    ).strip()

    base_inputs = json.loads(base_files["inputs"].read_text(encoding="utf-8"))
    revised = copy.deepcopy(json.loads(base_files["results"].read_text(encoding="utf-8")))
    if not revised["status"].startswith("confirmed_nominal_user_accepted_fitter_v8"):
        raise RuntimeError("base data-driven artifact is not the confirmed fitter-v8 nominal")

    correction_metadata = json.loads(correction_files["metadata"].read_text(encoding="utf-8"))
    if correction_metadata.get("combination_rule") != EXPECTED_RULE:
        raise RuntimeError("unexpected correction combination rule")
    with correction_files["sys3"].open(newline="", encoding="utf-8") as handle:
        correction_total_rows = [
            row for row in csv.DictReader(handle) if row["scope"] == "total"
        ]
    if len(correction_total_rows) != 1:
        raise RuntimeError("expected one total correction row")
    correction_total = correction_total_rows[0]
    correction_relative = float(correction_total["sys3"])
    if not close(
        correction_relative,
        max(
            abs(float(correction_total["sps_combined"])),
            abs(float(correction_total["dps_combined"])),
        ),
    ):
        raise RuntimeError("total correction is not the signed-product envelope")

    with recombination_path.open(newline="", encoding="utf-8") as handle:
        total_rows = [row for row in csv.DictReader(handle) if row["scope"] == "total"]
    if len(total_rows) != 1:
        raise RuntimeError("expected one total recombination row")
    total = total_rows[0]
    if not close(float(total["correction_relative"]), correction_relative):
        raise RuntimeError("recombination and correction artifact disagree")

    pair_total_relative = float(total["total_systematic_relative"])
    fitter_relative = float(total["fitter_relative"])
    pair_non_br_relative = quadrature(
        correction_relative,
        float(total["luminosity_relative"]),
        float(total["lifetime_relative"]),
        fitter_relative,
    )

    f_result = revised["f_dps"]
    method_up_relative = float(f_result["combined_syst_up_percent"]) / 100.0
    method_down_relative = float(f_result["combined_syst_down_percent"]) / 100.0
    sigma_up_relative = quadrature(method_up_relative, pair_total_relative)
    sigma_down_relative = quadrature(method_down_relative, pair_total_relative)

    sigma_cut = revised["sigma_dps_mjj_ge_7p5_pb"]
    sigma_cut_central = float(sigma_cut["central"])
    sigma_cut_up = sigma_cut_central * sigma_up_relative
    sigma_cut_down = sigma_cut_central * sigma_down_relative
    sigma_cut.update(
        {
            "included_syst_up": sigma_cut_up,
            "included_syst_down": sigma_cut_down,
            "included_syst_up_percent": 100.0 * sigma_up_relative,
            "included_syst_down_percent": 100.0 * sigma_down_relative,
            "total_syst": {"up": sigma_cut_up, "down": sigma_cut_down},
        }
    )

    sigma_no_mass = revised["sigma_dps_no_mjj_cut_pb"]
    sigma_no_mass_central = float(sigma_no_mass["central"])
    sigma_no_mass_up = sigma_no_mass_central * sigma_up_relative
    sigma_no_mass_down = sigma_no_mass_central * sigma_down_relative
    sigma_no_mass.update(
        {
            "included_syst_up": sigma_no_mass_up,
            "included_syst_down": sigma_no_mass_down,
            "included_syst_up_percent": 100.0 * sigma_up_relative,
            "included_syst_down_percent": 100.0 * sigma_down_relative,
        }
    )

    with base_files["uncertainty_components"].open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        component_fields = list(reader.fieldnames or [])
        components = list(reader)
    correction_component = find_row(
        components, "sigma_DPS_mJJ_ge_7p5", "sigma_total_acceptance_efficiency"
    )
    correction_effect = sigma_cut_central * correction_relative
    correction_component.update(
        {
            "effect_up": str(correction_effect),
            "effect_down": str(correction_effect),
            "treatment": "included_signed_multiplicative_SPS_DPS_envelope",
            "effect_up_percent": str(100.0 * correction_relative),
            "effect_down_percent": str(100.0 * correction_relative),
        }
    )
    sigma_cut_combined = find_row(
        components, "sigma_DPS_mJJ_ge_7p5", "combined", "syst_total"
    )
    sigma_cut_combined.update(
        {
            "effect_up": str(sigma_cut_up),
            "effect_down": str(sigma_cut_down),
            "effect_up_percent": str(100.0 * sigma_up_relative),
            "effect_down_percent": str(100.0 * sigma_down_relative),
        }
    )

    sigma_eff = revised["sigma_eff_mb"]
    sigma_eff_central = float(sigma_eff["central"])
    pair_component = find_row(components, "sigma_eff", "sigma_total_non_BR")
    pair_effect_up = inverse_up(sigma_eff_central, pair_non_br_relative)
    pair_effect_down = inverse_down(sigma_eff_central, pair_non_br_relative)
    pair_component.update(
        {
            "effect_up": str(pair_effect_up),
            "effect_down": str(pair_effect_down),
            "treatment": "included_with_multiplicative_correction_and_fitter_v8",
            "effect_up_percent": str(percent(pair_effect_up, sigma_eff_central)),
            "effect_down_percent": str(percent(pair_effect_down, sigma_eff_central)),
        }
    )
    method_component = find_row(components, "sigma_eff", "data_driven_method")
    atlas_bin = find_row(components, "sigma_eff", "ATLAS_single_bin")
    atlas_lumi = find_row(components, "sigma_eff", "ATLAS_luminosity")
    sigma_eff_up = quadrature(
        float(method_component["effect_up"]),
        pair_effect_up,
        float(atlas_bin["effect_up"]),
        float(atlas_lumi["effect_up"]),
    )
    sigma_eff_down = quadrature(
        float(method_component["effect_down"]),
        pair_effect_down,
        float(atlas_bin["effect_down"]),
        float(atlas_lumi["effect_down"]),
    )
    sigma_eff.update(
        {
            "included_syst_up": sigma_eff_up,
            "included_syst_down": sigma_eff_down,
            "included_syst_up_percent": percent(sigma_eff_up, sigma_eff_central),
            "included_syst_down_percent": percent(sigma_eff_down, sigma_eff_central),
        }
    )
    sigma_eff_combined = find_row(
        components, "sigma_eff", "combined_quantified", "syst_subtotal"
    )
    sigma_eff_combined.update(
        {
            "effect_up": str(sigma_eff_up),
            "effect_down": str(sigma_eff_down),
            "effect_up_percent": str(percent(sigma_eff_up, sigma_eff_central)),
            "effect_down_percent": str(percent(sigma_eff_down, sigma_eff_central)),
        }
    )

    revised["status"] = (
        "confirmed_nominal_user_accepted_fitter_v8_multiplicative_correction_"
        "epsilon_m_uncertainty_pending_not_zero"
    )
    revised["assumptions"]["correction_systematic"] = (
        "signed_multiplicative_acceptance_efficiency_closure_envelope"
    )

    with base_files["results_csv"].open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        result_fields = list(reader.fieldnames or [])
        result_rows = list(reader)
    result_updates = {
        ("sigma_DPS", "mJJ>=7.5"): (sigma_cut_up, sigma_cut_down),
        ("sigma_DPS", "no_mJJ_cut"): (sigma_no_mass_up, sigma_no_mass_down),
        ("sigma_eff", "single_Jpsi_matched_no_pair_mass_cut"): (
            sigma_eff_up,
            sigma_eff_down,
        ),
    }
    for row in result_rows:
        key = (row["quantity"], row["phase_space"])
        if key not in result_updates:
            continue
        up, down = result_updates[key]
        central = float(row["central"])
        row["included_syst_up"] = up
        row["included_syst_down"] = down
        row["included_syst_up_percent"] = percent(up, central)
        row["included_syst_down_percent"] = percent(down, central)

    source_hashes = {name: sha256(path) for name, path in sources.items()}
    revised_inputs = {
        "status": "frozen_snapshot_for_confirmed_multiplicative_correction_propagation",
        "base_confirmed_artifact": {
            "path": str(base_dir.resolve()),
            "artifact_checksums_verified": True,
        },
        "correction_systematic": {
            "status": "accepted_by_user_and_included",
            "rule": EXPECTED_RULE,
            "relative": correction_relative,
            "percent": 100.0 * correction_relative,
            "sps_combined_signed": float(correction_total["sps_combined"]),
            "dps_combined_signed": float(correction_total["dps_combined"]),
            "envelope_model": correction_total["envelope_model"],
            "pair_total_systematic_relative": pair_total_relative,
        },
        "fitter_stability": {
            **base_inputs["fitter_stability"],
            "total_systematic_relative_after_recombination": pair_total_relative,
        },
        "epsilon_m_uncertainty": {
            "value": None,
            "status": "pending_not_included_not_zero",
        },
        "inherited_inputs": base_inputs,
        "source_paths": {name: str(path.resolve()) for name, path in sources.items()},
        "source_sha256": source_hashes,
    }

    summary = [
        f"status={revised['status']}",
        "promotion_status=confirmed_nominal_user_accepted",
        "correction_combination=signed_multiplicative_acceptance_efficiency_closure_envelope",
        f"correction_sps_combined_signed={float(correction_total['sps_combined']):.12g}",
        f"correction_dps_combined_signed={float(correction_total['dps_combined']):.12g}",
        f"correction_envelope_model={correction_total['envelope_model']}",
        f"correction_systematic_relative={correction_relative:.12g}",
        f"correction_systematic_percent={100.0 * correction_relative:.12g}",
        f"sigma_total_systematic_relative={pair_total_relative:.12g}",
        f"sigma_total_systematic_percent={100.0 * pair_total_relative:.12g}",
        f"f_dps={float(f_result['central']):.12g}",
        f"f_dps_stat={float(f_result['stat']):.12g}",
        f"f_dps_method_syst_up={float(f_result['combined_syst_up']):.12g}",
        f"f_dps_method_syst_down={float(f_result['combined_syst_down']):.12g}",
        f"sigma_dps_mjj_ge_7p5_pb={sigma_cut_central:.12g}",
        f"sigma_dps_mjj_ge_7p5_stat_pb={float(sigma_cut['stat']):.12g}",
        f"sigma_dps_mjj_ge_7p5_syst_up_pb={sigma_cut_up:.12g}",
        f"sigma_dps_mjj_ge_7p5_syst_down_pb={sigma_cut_down:.12g}",
        f"sigma_dps_no_mjj_cut_pb={sigma_no_mass_central:.12g}",
        f"sigma_dps_no_mjj_cut_stat_pb={float(sigma_no_mass['stat']):.12g}",
        f"sigma_dps_no_mjj_cut_syst_up_pb={sigma_no_mass_up:.12g}",
        f"sigma_dps_no_mjj_cut_syst_down_pb={sigma_no_mass_down:.12g}",
        f"sigma_eff_mb={sigma_eff_central:.12g}",
        f"sigma_eff_stat_mb={float(sigma_eff['stat']):.12g}",
        f"sigma_eff_quantified_syst_up_mb={sigma_eff_up:.12g}",
        f"sigma_eff_quantified_syst_down_mb={sigma_eff_down:.12g}",
        "epsilon_m_uncertainty_status=pending_not_included_not_zero",
        "new_roofit_fits_run=0",
    ]
    metadata = [
        "status=complete_confirmed_nominal_epsilon_m_uncertainty_pending_not_zero",
        f"output_tag={output_dir.name}",
        f"generated_at={datetime.now().astimezone().isoformat(timespec='seconds')}",
        f"root_version={root_version}",
        f"git_head={git_head}",
        "record_type=deterministic_correction_systematic_revision_no_ROOT_fits",
        "base_artifact_checksums_verified=true",
        "correction_artifact_checksums_verified=true",
        "user_gate=2026-09-03_signed_multiplicative_correction_nominal",
        "epsilon_m_uncertainty=pending_not_included_not_zero",
    ]

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.tmp.", dir=output_dir.parent))
    try:
        write_json(staging / "inputs.json", revised_inputs)
        write_json(staging / "results.json", revised)
        write_csv(staging / "results.csv", result_fields, result_rows)
        write_csv(staging / "uncertainty_components.csv", component_fields, components)
        (staging / "summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
        (staging / "run.metadata.txt").write_text("\n".join(metadata) + "\n", encoding="utf-8")
        (staging / "source_checksums.txt").write_text(
            "".join(
                f"{source_hashes[name]}  {sources[name].resolve()}\n"
                for name in sorted(sources)
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
        staging.rename(output_dir)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise

    print(f"created={output_dir}")
    print(f"correction_systematic={100.0 * correction_relative:.12g}%")
    print(f"sigma_total_systematic={100.0 * pair_total_relative:.12g}%")
    print(f"sigma_DPS_syst=+{sigma_cut_up:.12g}/-{sigma_cut_down:.12g} pb")
    print(f"sigma_eff_syst=+{sigma_eff_up:.12g}/-{sigma_eff_down:.12g} mb")


if __name__ == "__main__":
    main()
