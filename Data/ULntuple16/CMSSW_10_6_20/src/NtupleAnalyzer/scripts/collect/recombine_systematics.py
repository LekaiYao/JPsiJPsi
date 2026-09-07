#!/usr/bin/env python3
"""Build the formal Data systematic table from frozen fit and closure inputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import re
import shutil
import sys
from datetime import datetime
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
ANALYSIS_DIR = SCRIPT_DIR.parents[1]
REPO = ANALYSIS_DIR.parents[4]
FIT_RESULTS = ANALYSIS_DIR / "fit_results"
SYSTEMATICS_RUNS = FIT_RESULTS / "runs/systematics"
DEFAULT_BASE = (
    FIT_RESULTS / "_CURRENT/differential_cross_sections/cross_sections.csv"
)
DEFAULT_CORRECTION_DIR = (
    REPO
    / "GEN_nofilter/DPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc/closure_results"
    / "nominal_mixed23_sps0p85_dps0p15_multiplicative_v3_20260903"
)
DEFAULT_FITTER_DIR = (
    REPO
    / "docs/answers/H018_recalculate_fitter_stability_with_legacy_variations/machine"
)
DEFAULT_FITTER_ARTIFACT = FIT_RESULTS / "_CURRENT/fitter_stability"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Recombine the confirmed Data systematic components."
    )
    parser.add_argument("tag", help="new tag in the form systematics-YYYYMMDD-rNN")
    parser.add_argument("--base", type=Path, default=DEFAULT_BASE)
    parser.add_argument(
        "--correction-dir", type=Path, default=DEFAULT_CORRECTION_DIR
    )
    parser.add_argument("--fitter-dir", type=Path, default=DEFAULT_FITTER_DIR)
    parser.add_argument(
        "--fitter-artifact", type=Path, default=DEFAULT_FITTER_ARTIFACT
    )
    parser.add_argument("--accepted-at", default=datetime.now().date().isoformat())
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"refusing to write an empty table: {path}")
    with path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def bin_key(variable: str, low: str, high: str) -> tuple[str, float, float]:
    return variable, round(float(low), 8), round(float(high), 8)


def optional_float(value: str) -> float | None:
    return None if value == "" else float(value)


def close(left: float, right: float, tolerance: float = 5.0e-13) -> bool:
    return math.isclose(left, right, rel_tol=tolerance, abs_tol=tolerance)


def validate_correction(rows: list[dict[str, str]]) -> None:
    if len(rows) != 36:
        raise RuntimeError(f"expected 36 correction rows, found {len(rows)}")
    for row in rows:
        candidates: dict[str, float] = {}
        for model in ("sps", "dps"):
            acceptance = optional_float(row[f"{model}_acceptance"])
            efficiency = optional_float(row[f"{model}_efficiency"])
            combined = optional_float(row[f"{model}_combined"])
            if acceptance is None:
                if combined is not None:
                    raise RuntimeError(f"unexpected {model} combined value: {row}")
                continue
            if efficiency is None or combined is None:
                raise RuntimeError(f"incomplete {model} correction inputs: {row}")
            expected = (1.0 + acceptance) * (1.0 + efficiency) - 1.0
            if not close(combined, expected):
                raise RuntimeError(
                    f"signed multiplicative closure mismatch for {model}: {row}"
                )
            candidates[model] = combined
        if not candidates:
            raise RuntimeError(f"no correction candidate: {row}")
        envelope_model = max(candidates, key=lambda model: abs(candidates[model]))
        envelope = abs(candidates[envelope_model])
        if row["envelope_model"] != envelope_model or not close(
            float(row["sys3"]), envelope
        ):
            raise RuntimeError(f"incorrect SPS/DPS correction envelope: {row}")


def snapshot(source: Path, target: Path) -> None:
    if not source.is_file():
        raise FileNotFoundError(source)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def main() -> None:
    args = parse_args()
    if len(args.tag) > 48 or not re.fullmatch(
        r"systematics-[0-9]{8}-r[0-9]{2}", args.tag
    ):
        raise SystemExit(
            "tag must match systematics-YYYYMMDD-rNN and be at most 48 characters"
        )
    output = SYSTEMATICS_RUNS / args.tag
    if output.exists() or output.is_symlink():
        raise SystemExit(f"refusing to overwrite output: {output}")

    correction_csv = args.correction_dir / "sys3.csv"
    total_fitter_csv = args.fitter_dir / "total_fitter_variations.csv"
    differential_fitter_csv = (
        args.fitter_dir / "differential_fitter_variations.csv"
    )
    sources = {
        "base_cross_sections.csv": args.base,
        "correction_sys3.csv": correction_csv,
        "correction_closure_components.csv": (
            args.correction_dir / "closure_components.csv"
        ),
        "correction_metadata.json": args.correction_dir / "metadata.json",
        "correction_checksums.sha256": (
            args.correction_dir / "checksums.sha256"
        ),
        "fitter_total_variations.csv": total_fitter_csv,
        "fitter_differential_variations.csv": differential_fitter_csv,
        "fitter_stability_run.metadata.txt": (
            args.fitter_artifact / "run.metadata.txt"
        ),
        "fitter_stability_checksums.sha256": (
            args.fitter_artifact / "checksums.sha256"
        ),
    }
    for source in sources.values():
        if not source.is_file():
            raise FileNotFoundError(source)

    base_rows = read_rows(args.base)
    correction_rows = read_rows(correction_csv)
    total_fitter_rows = read_rows(total_fitter_csv)
    differential_fitter_rows = read_rows(differential_fitter_csv)
    if len(base_rows) != 36:
        raise RuntimeError(f"expected 36 base rows, found {len(base_rows)}")
    validate_correction(correction_rows)

    total_envelopes = [
        row for row in total_fitter_rows if int(row["is_envelope"]) == 1
    ]
    if len(total_envelopes) != 1:
        raise RuntimeError("expected exactly one total fitter envelope")
    total_envelope = total_envelopes[0]
    if total_envelope["variation"] != "f3_single_gaussian":
        raise RuntimeError("the confirmed total fitter envelope is not f3_single_gaussian")
    total_fitter = float(total_envelope["relative_shift"])

    differential_fitter: dict[tuple[str, float, float], float] = {}
    for row in differential_fitter_rows:
        if row["variation"] != "f3_single_gaussian":
            raise RuntimeError(f"unexpected differential fitter model: {row}")
        key = bin_key(row["variable"], row["bin_min"], row["bin_max"])
        if key in differential_fitter:
            raise RuntimeError(f"duplicate differential fitter bin: {key}")
        differential_fitter[key] = float(row["fitter_relative"])
    if len(differential_fitter) != 35:
        raise RuntimeError(
            f"expected 35 differential fitter bins, found {len(differential_fitter)}"
        )

    total_correction_rows = [
        row for row in correction_rows if row["scope"] == "total"
    ]
    if len(total_correction_rows) != 1:
        raise RuntimeError("expected exactly one total correction row")
    total_correction = float(total_correction_rows[0]["sys3"])
    differential_correction = {
        bin_key(row["variable"], row["low"], row["high"]): float(row["sys3"])
        for row in correction_rows
        if row["scope"] == "bin"
    }
    if len(differential_correction) != 35:
        raise RuntimeError(
            f"expected 35 differential correction bins, found {len(differential_correction)}"
        )

    result_rows: list[dict[str, object]] = []
    for base in base_rows:
        if base["scope"] == "total":
            correction = total_correction
            fitter = total_fitter
        elif base["scope"] == "differential":
            key = bin_key(base["variable"], base["bin_min"], base["bin_max"])
            if key not in differential_correction or key not in differential_fitter:
                raise RuntimeError(f"missing systematic input for bin {key}")
            correction = differential_correction[key]
            fitter = differential_fitter[key]
        else:
            raise RuntimeError(f"unexpected base scope: {base['scope']}")
        branching = float(base["br_relative"])
        luminosity = float(base["luminosity_relative"])
        lifetime = float(base["lifetime_relative"])
        total_systematic = math.sqrt(
            branching**2
            + luminosity**2
            + correction**2
            + lifetime**2
            + fitter**2
        )
        cross_section = float(base["sigma_pb_per_unit"])
        result_rows.append(
            {
                "scope": base["scope"],
                "variable": base["variable"],
                "bin_min": base["bin_min"],
                "bin_max": base["bin_max"],
                "sigma_pb_per_unit": f"{cross_section:.15g}",
                "stat_pb_per_unit": base["stat_pb_per_unit"],
                "fitter_relative": f"{fitter:.15g}",
                "br_relative": f"{branching:.15g}",
                "luminosity_relative": f"{luminosity:.15g}",
                "correction_relative": f"{correction:.15g}",
                "lifetime_relative": f"{lifetime:.15g}",
                "total_systematic_relative": f"{total_systematic:.15g}",
                "total_systematic_pb_per_unit": (
                    f"{cross_section * total_systematic:.15g}"
                ),
            }
        )

    output.mkdir(parents=True)
    for name, source in sources.items():
        snapshot(source, output / "inputs" / name)
    snapshot(Path(__file__), output / "code" / Path(__file__).name)
    result_csv = output / "cross_sections.csv"
    write_rows(result_csv, result_rows)

    total = result_rows[0]
    metadata = [
        "status=complete",
        "promotion_status=user_confirmed_formal_systematic_source",
        f"accepted_at={args.accepted_at}",
        f"run_tag={args.tag}",
        "rows=36 (one total plus 35 differential bins)",
        "central_and_statistical_source=base corrected-error 4D fit table",
        "correction_component_source=project closure residuals in correction_sys3.csv",
        "correction_candidate=(1+signed_acceptance_residual)*(1+signed_efficiency_residual)-1 separately for SPS and DPS",
        "correction_envelope=max(abs(SPS_candidate),abs(DPS_candidate)) over available models",
        "sparse_bin_policy=when SPS acceptance closure is unavailable, use the complete DPS candidate; never substitute zero",
        "fitter_total=max absolute n_P_P shift among eight approved one-at-a-time variations",
        "fitter_total_envelope=f3_single_gaussian",
        "fitter_differential=per-bin absolute n_P_P shift for f3_single_gaussian",
        "combination=quadrature of branching,luminosity,correction,lifetime,fitter relative components",
        f"total_sigma_pb={total['sigma_pb_per_unit']}",
        f"total_stat_pb={total['stat_pb_per_unit']}",
        f"total_correction_relative={total['correction_relative']}",
        f"total_fitter_relative={total['fitter_relative']}",
        f"total_systematic_relative={total['total_systematic_relative']}",
        f"total_systematic_pb={total['total_systematic_pb_per_unit']}",
        f"python_version={sys.version.split()[0]}",
        f"generated_at={datetime.now().astimezone().isoformat(timespec='seconds')}",
    ]
    for name, source in sources.items():
        metadata.append(f"source_{name}_sha256={sha256(source)}")
    metadata.append(f"sha256 cross_sections.csv={sha256(result_csv)}")
    metadata.append(f"sha256 recombine_systematics.py={sha256(Path(__file__))}")
    (output / "run.metadata.txt").write_text(
        "\n".join(metadata) + "\n", encoding="utf-8"
    )
    (output / "complete.marker").write_text(
        "status=complete_user_confirmed_formal_systematic_source\n",
        encoding="utf-8",
    )

    checksum_targets = sorted(
        path for path in output.rglob("*") if path.is_file()
    )
    with (output / "checksums.sha256").open("x", encoding="utf-8") as handle:
        for path in checksum_targets:
            handle.write(f"{sha256(path)}  {path.relative_to(output)}\n")
    print(output)


if __name__ == "__main__":
    main()
