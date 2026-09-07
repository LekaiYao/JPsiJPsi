#!/usr/bin/env python3
import csv
import hashlib
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path


BINS = [
    ("delta_y", 0, 0.5), ("delta_y", 0.5, 1), ("delta_y", 1, 1.5),
    ("delta_y", 1.5, 2), ("delta_y", 2, 2.5), ("delta_y", 2.5, 4),
    ("delta_phi", 0, 0.3927), ("delta_phi", 0.3927, 0.7854),
    ("delta_phi", 0.7854, 1.1781), ("delta_phi", 1.1781, 1.5708),
    ("delta_phi", 1.5708, 1.9635), ("delta_phi", 1.9635, 2.3562),
    ("delta_phi", 2.3562, 2.7489), ("delta_phi", 2.7489, 3.1416),
    ("evt_mass", 7.5, 17.5), ("evt_mass", 17.5, 27.5),
    ("evt_mass", 27.5, 37.5), ("evt_mass", 37.5, 47.5),
    ("evt_mass", 47.5, 57.5), ("evt_mass", 57.5, 67.5),
    ("evt_mass", 67.5, 107.5), ("evt_y", 0, 0.4),
    ("evt_y", 0.4, 0.8), ("evt_y", 0.8, 1.2),
    ("evt_y", 1.2, 1.6), ("evt_y", 1.6, 2),
    ("evt_pt", 0, 5), ("evt_pt", 5, 10), ("evt_pt", 10, 15),
    ("evt_pt", 15, 20), ("evt_pt", 20, 25), ("evt_pt", 25, 30),
    ("evt_pt", 30, 35), ("evt_pt", 35, 40), ("evt_pt", 40, 80),
]

HEADER = [
    "scope", "variable", "bin_min", "bin_max", "yield", "stat_yield",
    "reference_yield", "fitter_relative", "sigma_pb_per_unit",
    "stat_pb_per_unit", "fitter_pb_per_unit", "br_relative",
    "luminosity_relative", "correction_relative", "lifetime_relative",
    "total_systematic_relative", "total_systematic_pb_per_unit",
    "nominal_n_Sig_Comb_n_Comb_Sig_shared_fallback",
    "nominal_n_Comb_Comb_fallback",
    "reference_n_Sig_Comb_n_Comb_Sig_shared_fallback",
    "reference_n_Comb_Comb_fallback", "status", "covQual", "edm",
    "reference_status", "reference_covQual", "reference_edm",
]


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def key_values(path):
    values = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key] = value
    return values


def read_single_row(path):
    rows = list(csv.reader(path.open(newline="")))
    if len(rows) != 1 or len(rows[0]) != len(HEADER):
        raise RuntimeError(f"Malformed row file: {path}")
    return rows[0]


def validate_fit_row(row, scope, variable, low, high):
    if row[0] != scope or row[1] != variable:
        raise RuntimeError(f"Wrong row identity: {row[:4]}")
    if abs(float(row[2]) - low) > 1e-9 or abs(float(row[3]) - high) > 1e-9:
        raise RuntimeError(f"Wrong bin range: {row[:4]}")
    if float(row[4]) <= 0 or float(row[5]) <= 0:
        raise RuntimeError(f"Non-positive yield or error: {row[:6]}")
    if int(row[21]) != 0 or int(row[22]) != 3 or float(row[23]) >= 0.01:
        raise RuntimeError(f"Nominal fit gate failed: {row[:4]}")
    if int(row[24]) != 0 or int(row[25]) != 3 or float(row[26]) >= 0.01:
        raise RuntimeError(f"Reference fit gate failed: {row[:4]}")


def main():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: collect_differential_newdata167.py CAMPAIGN_TAG")
    campaign_tag = sys.argv[1]
    analysis_dir = Path(__file__).resolve().parent
    output_dir = analysis_dir / "fit_results" / campaign_tag
    output_csv = output_dir / "cross_sections.csv"
    metadata_path = output_dir / "run.metadata.txt"
    if output_csv.exists() or metadata_path.exists():
        raise SystemExit(f"Refusing to overwrite completed output: {output_dir}")

    total = read_single_row(output_dir / "total_row.csv")
    validate_fit_row(total, "total", "total", 0.0, 0.0)
    rows = []
    for index, (variable, low, high) in enumerate(BINS):
        metadata = key_values(output_dir / "jobs" / f"{index:02d}.metadata.txt")
        if metadata.get("status") != "complete":
            raise RuntimeError(f"Incomplete job metadata for index {index}")
        row = read_single_row(output_dir / "rows" / f"{index:02d}.csv")
        validate_fit_row(row, "differential", variable, low, high)
        rows.append(row)

    with output_csv.open("x", newline="") as output:
        writer = csv.writer(output)
        writer.writerow(HEADER)
        writer.writerow(total)
        writer.writerows(rows)

    integrals = defaultdict(float)
    fallback_counts = defaultdict(int)
    for row in rows:
        variable = row[1]
        integrals[variable] += float(row[8]) * (float(row[3]) - float(row[2]))
        fallback_counts["nominal_sigcomb"] += int(row[17])
        fallback_counts["nominal_combcomb"] += int(row[18])
        fallback_counts["reference_sigcomb"] += int(row[19])
        fallback_counts["reference_combcomb"] += int(row[20])

    source_sys3 = Path(
        "/eos/home-l/leyao/26JJ/JPsiJPsi/GEN_nofilter/DPS/"
        "ULPythia2016/CMSSW_10_2_5/src/4mu_acc/closure_results/"
        "nominal_mixed23_sps0p85_dps0p15_avgacc_v2_20260902/sys3_cpp.txt"
    )
    lines = [
        "status=complete",
        f"campaign_tag={campaign_tag}",
        "completed_bins=35",
        "fits=70 nominal/reference differential fits plus frozen total nominal/reference",
        "error_convention=ROOT native AsymptoticError(true), no external seed",
        "fit_gate=status=0,covQual=3,EDM<0.01",
        "luminosity_fb=36.684",
        "systematic_status=base_fit_table_only; formal result is fit_results/_CURRENT/systematics",
        "base_total_correction_systematic_relative=0.08805293089840492",
        f"base_sys3_source={source_sys3}",
        f"base_sys3_sha256={sha256(source_sys3)}",
        "base_correction_fallback_note=four sparse bins use the accepted DPS-only fallback; undefined SPS closure is not treated as zero",
        f"nominal_sigcomb_fallback_bins={fallback_counts['nominal_sigcomb']}",
        f"nominal_combcomb_fallback_bins={fallback_counts['nominal_combcomb']}",
        f"reference_sigcomb_fallback_bins={fallback_counts['reference_sigcomb']}",
        f"reference_combcomb_fallback_bins={fallback_counts['reference_combcomb']}",
    ]
    for variable in ("delta_y", "delta_phi", "evt_mass", "evt_y", "evt_pt"):
        lines.append(f"{variable}_integrated_sigma_pb={integrals[variable]:.17g}")
    lines.extend([
        f"generated_at={datetime.now().astimezone().isoformat(timespec='seconds')}",
        f"sha256 cross_sections.csv={sha256(output_csv)}",
        f"sha256 Fit_4D_diff.cpp={sha256(analysis_dir / 'Fit_4D_diff.cpp')}",
        f"sha256 Fit_valid.cpp={sha256(analysis_dir / 'Fit_valid.cpp')}",
        f"sha256 Plot_4D.hpp={sha256(analysis_dir / 'Plot_4D.hpp')}",
        f"sha256 worker={sha256(analysis_dir / 'run_differential_newdata167_condor.sh')}",
        f"sha256 collector={sha256(Path(__file__))}",
        "",
    ])
    metadata_path.write_text("\n".join(lines))
    print(output_dir)


if __name__ == "__main__":
    main()
