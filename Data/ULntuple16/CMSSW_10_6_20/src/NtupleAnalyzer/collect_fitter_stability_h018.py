#!/usr/bin/env python3
"""Collect the H018 legacy-variation fitter-stability results."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


VARIATIONS = [
    ("f1_double_gaussian", "f1->f1^(1)", "mass signal: double Gaussian"),
    ("f1_double_crystal_ball", "f1->f1^(2)", "mass signal: double Crystal Ball"),
    ("f1_triple_gaussian", "f1->f1^(3)", "mass signal: triple Gaussian"),
    ("f2_linear", "f2->f2^(1)", "mass combinatorial background: linear"),
    ("f3_single_gaussian", "f3->f3^(1)", "prompt lifetime: single Gaussian"),
    ("f3_triple_gaussian", "f3->f3^(2)", "prompt lifetime: triple Gaussian"),
    (
        "f4_exp_double_gaussian",
        "f4->f4^(1)",
        "nonprompt lifetime: exponential convolved with double Gaussian",
    ),
    (
        "f4_double_exp_gaussian",
        "f4->f4^(2)",
        "nonprompt lifetime: double exponential convolved with Gaussian",
    ),
]


def read_single_row(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1:
        raise RuntimeError(f"expected exactly one row in {path}, found {len(rows)}")
    return rows[0]


def write_rows(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def bin_key(row: dict[str, str]) -> tuple[str, float, float]:
    return (row["variable"], float(row["bin_min"]), float(row["bin_max"]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--nominal-csv", required=True, type=Path)
    parser.add_argument("--answer-dir", required=True, type=Path)
    parser.add_argument("--differential-reference", required=True)
    args = parser.parse_args()

    with args.nominal_csv.open(newline="") as handle:
        nominal_rows = list(csv.DictReader(handle))
    total_nominal_rows = [row for row in nominal_rows if row["scope"] == "total"]
    differential_nominal_rows = [
        row for row in nominal_rows if row["scope"] == "differential"
    ]
    if len(total_nominal_rows) != 1 or len(differential_nominal_rows) != 35:
        raise RuntimeError("nominal H017 table must contain one total and 35 differential rows")
    total_nominal = total_nominal_rows[0]
    nominal_yield = float(total_nominal["yield"])
    if not math.isclose(nominal_yield, 5687.8051561631491, rel_tol=0.0, abs_tol=1e-9):
        raise RuntimeError(f"unexpected frozen nominal total yield: {nominal_yield}")
    if (
        int(total_nominal["status"]) != 0
        or int(total_nominal["covQual"]) != 3
        or float(total_nominal["edm"]) >= 0.01
    ):
        raise RuntimeError("frozen nominal total fit fails the strict H018 gate")

    total_rows: list[dict[str, object]] = []
    maximum_shift = -1.0
    maximum_variation = ""
    for variation, an_label, description in VARIATIONS:
        row = read_single_row(args.run_dir / "total" / variation / "summary.csv")
        status = int(row["status"])
        covariance_quality = int(row["covQual"])
        optimizer_edm = float(row["optimizer_edm"])
        post_asymptotic_edm = float(row["post_asymptotic_edm"])
        post_covariance_exception = "post_covariance_edm" in row["recovery"]
        if (
            status != 0
            or covariance_quality != 3
            or optimizer_edm >= 0.01
            or (post_asymptotic_edm >= 0.01 and not post_covariance_exception)
            or (post_asymptotic_edm < 0.01 and post_covariance_exception)
        ):
            raise RuntimeError(
                f"total variation {variation} fails the approved H018 gate"
            )
        variation_yield = float(row["yield"])
        absolute_shift = abs(variation_yield - nominal_yield)
        relative_shift = absolute_shift / nominal_yield
        if relative_shift > maximum_shift:
            maximum_shift = relative_shift
            maximum_variation = variation
        total_rows.append(
            {
                "variation": variation,
                "an_label": an_label,
                "description": description,
                "nominal_yield": f"{nominal_yield:.15g}",
                "variation_yield": f"{variation_yield:.15g}",
                "variation_stat_error": f"{float(row['yield_error']):.15g}",
                "absolute_shift": f"{absolute_shift:.15g}",
                "relative_shift": f"{relative_shift:.15g}",
                "relative_shift_percent": f"{100.0 * relative_shift:.12g}",
                "status": status,
                "covQual": covariance_quality,
                "optimizer_edm": f"{optimizer_edm:.15g}",
                "post_asymptotic_edm": f"{post_asymptotic_edm:.15g}",
                "post_covariance_edm_exception": int(post_covariance_exception),
                "recovery": row["recovery"],
                "is_envelope": 0,
            }
        )
    for row in total_rows:
        row["is_envelope"] = int(row["variation"] == maximum_variation)
    if args.differential_reference != maximum_variation:
        raise RuntimeError(
            "differential reference does not match the total fitter envelope: "
            f"{args.differential_reference} != {maximum_variation}"
        )

    diff_files = sorted((args.run_dir / "differential").glob("*.root.csv"))
    if len(diff_files) != 35:
        raise RuntimeError(f"expected 35 H018 differential result rows, found {len(diff_files)}")
    variation_by_bin = {bin_key(read_single_row(path)): read_single_row(path) for path in diff_files}
    if len(variation_by_bin) != 35:
        raise RuntimeError("duplicate H018 differential bin keys")

    differential_rows: list[dict[str, object]] = []
    recombination_rows: list[dict[str, object]] = []
    for nominal in differential_nominal_rows:
        key = bin_key(nominal)
        if key not in variation_by_bin:
            raise RuntimeError(f"missing H018 differential variation for bin {key}")
        variation = variation_by_bin[key]
        status = int(variation["status"])
        covariance_quality = int(variation["covQual"])
        edm = float(variation["edm"])
        if status != 0 or covariance_quality != 3 or edm >= 0.01:
            raise RuntimeError(f"differential variation {key} fails the strict H018 gate")
        nominal_bin_yield = float(nominal["yield"])
        variation_yield = float(variation["yield"])
        absolute_shift = abs(variation_yield - nominal_bin_yield)
        relative_shift = absolute_shift / nominal_bin_yield
        differential_rows.append(
            {
                "variable": key[0],
                "bin_min": f"{key[1]:.12g}",
                "bin_max": f"{key[2]:.12g}",
                "nominal_yield": f"{nominal_bin_yield:.15g}",
                "nominal_stat_error": f"{float(nominal['stat_yield']):.15g}",
                "variation": args.differential_reference,
                "variation_yield": f"{variation_yield:.15g}",
                "variation_stat_error": f"{float(variation['yield_error']):.15g}",
                "absolute_shift": f"{absolute_shift:.15g}",
                "fitter_relative": f"{relative_shift:.15g}",
                "fitter_percent": f"{100.0 * relative_shift:.12g}",
                "nominal_status": nominal["status"],
                "nominal_covQual": nominal["covQual"],
                "nominal_edm": nominal["edm"],
                "variation_status": status,
                "variation_covQual": covariance_quality,
                "variation_edm": f"{edm:.15g}",
                "nominal_n_Sig_Comb_fallback": nominal[
                    "nominal_n_Sig_Comb_n_Comb_Sig_shared_fallback"
                ],
                "nominal_n_Comb_Comb_fallback": nominal[
                    "nominal_n_Comb_Comb_fallback"
                ],
                "variation_n_Sig_Comb_fallback": variation[
                    "n_Sig_Comb_fallback"
                ],
                "variation_n_Comb_Comb_fallback": variation[
                    "n_Comb_Comb_fallback"
                ],
                "variation_strategy": variation["strategy"],
                "variation_offset": variation["offset"],
            }
        )
        br = float(nominal["br_relative"])
        lumi = float(nominal["luminosity_relative"])
        correction = float(nominal["correction_relative"])
        lifetime = float(nominal["lifetime_relative"])
        combined = math.sqrt(
            br * br
            + lumi * lumi
            + correction * correction
            + lifetime * lifetime
            + relative_shift * relative_shift
        )
        cross_section = float(nominal["sigma_pb_per_unit"])
        recombination_rows.append(
            {
                "scope": "differential",
                "variable": key[0],
                "bin_min": f"{key[1]:.12g}",
                "bin_max": f"{key[2]:.12g}",
                "cross_section_pb_per_unit": f"{cross_section:.15g}",
                "stat_pb_per_unit": nominal["stat_pb_per_unit"],
                "br_relative": f"{br:.15g}",
                "luminosity_relative": f"{lumi:.15g}",
                "correction_relative": f"{correction:.15g}",
                "lifetime_relative": f"{lifetime:.15g}",
                "fitter_relative": f"{relative_shift:.15g}",
                "total_systematic_relative": f"{combined:.15g}",
                "total_systematic_pb_per_unit": f"{cross_section * combined:.15g}",
            }
        )

    total_br = float(total_nominal["br_relative"])
    total_lumi = float(total_nominal["luminosity_relative"])
    total_correction = float(total_nominal["correction_relative"])
    total_lifetime = float(total_nominal["lifetime_relative"])
    total_combined = math.sqrt(
        total_br * total_br
        + total_lumi * total_lumi
        + total_correction * total_correction
        + total_lifetime * total_lifetime
        + maximum_shift * maximum_shift
    )
    total_cross_section = float(total_nominal["sigma_pb_per_unit"])
    recombination_rows.insert(
        0,
        {
            "scope": "total",
            "variable": "total",
            "bin_min": "0",
            "bin_max": "0",
            "cross_section_pb_per_unit": f"{total_cross_section:.15g}",
            "stat_pb_per_unit": total_nominal["stat_pb_per_unit"],
            "br_relative": f"{total_br:.15g}",
            "luminosity_relative": f"{total_lumi:.15g}",
            "correction_relative": f"{total_correction:.15g}",
            "lifetime_relative": f"{total_lifetime:.15g}",
            "fitter_relative": f"{maximum_shift:.15g}",
            "total_systematic_relative": f"{total_combined:.15g}",
            "total_systematic_pb_per_unit": f"{total_cross_section * total_combined:.15g}",
        },
    )

    machine = args.answer_dir / "machine"
    write_rows(
        machine / "total_fitter_variations.csv",
        list(total_rows[0].keys()),
        total_rows,
    )
    write_rows(
        machine / "differential_fitter_variations.csv",
        list(differential_rows[0].keys()),
        differential_rows,
    )
    write_rows(
        machine / "systematic_recombination.csv",
        list(recombination_rows[0].keys()),
        recombination_rows,
    )
    print(
        "H018_COLLECTION_COMPLETE,"
        f"envelope={maximum_variation},"
        f"fitter_relative={maximum_shift:.17g},"
        f"total_systematic_relative={total_combined:.17g}"
    )


if __name__ == "__main__":
    main()
