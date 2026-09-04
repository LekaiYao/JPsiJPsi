#!/usr/bin/env python3
import argparse
import csv
import hashlib
import re
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


PROJECTIONS = ["Jpsi_mass1", "Jpsi_mass2", "Jpsi_ctau1", "Jpsi_ctau2"]
COLORS = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd"]
MARKERS = ["o", "s", "^", "D"]
BIN_COUNTS = list(range(10, 61, 5))


def parse_key_values(path):
    values = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key] = value
    return values


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(
        description="Validate, aggregate, and plot the integral-based chi2 binning scan."
    )
    parser.add_argument("campaign_tag")
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9_-]+", args.campaign_tag):
        raise SystemExit("Invalid campaign tag")

    analysis_dir = Path(__file__).resolve().parent
    fit_results = analysis_dir / "fit_results"
    output_dir = fit_results / args.campaign_tag
    if output_dir.exists():
        raise SystemExit(f"Refusing to overwrite existing summary: {output_dir}")

    rows = []
    symmetrization_modes = set()
    sumw2_modes = set()
    for bins in BIN_COUNTS:
        result_tag = f"{args.campaign_tag}_bins{bins}"
        result_dir = fit_results / result_tag
        metadata = parse_key_values(result_dir / "job.metadata.txt")
        if metadata.get("status") != "complete":
            raise RuntimeError(f"Incomplete worker metadata: {result_dir}")
        if int(metadata.get("projection_bins", "-1")) != bins:
            raise RuntimeError(f"Worker bin count mismatch: {result_dir}")
        symmetrization_modes.add(
            metadata.get("projection_pair_symmetrization", "disabled")
        )
        sumw2_modes.add(metadata.get("projection_sumw2", "standard event-level"))

        configuration = parse_key_values(result_dir / "configuration.txt")
        expected_total = float(configuration["nominal_expected_events"])

        with (result_dir / "observed_chi2.csv").open(newline="") as source:
            observed = {row["projection"]: row for row in csv.DictReader(source)}
        if set(observed) != set(PROJECTIONS):
            raise RuntimeError(f"Projection set mismatch: {result_dir}")

        expected_by_projection = defaultdict(list)
        with (result_dir / "projection_bins.csv").open(newline="") as source:
            for row in csv.DictReader(source):
                expected_by_projection[row["projection"]].append(float(row["expected"]))

        for projection in PROJECTIONS:
            expected = expected_by_projection[projection]
            if len(expected) != bins:
                raise RuntimeError(
                    f"{result_tag}: {projection} has {len(expected)} bins, expected {bins}"
                )
            expected_sum = sum(expected)
            closure_relative = abs(expected_sum - expected_total) / expected_total
            target = expected_total / bins
            maximum_equal_probability_deviation = max(
                abs(value - target) / target for value in expected
            )
            if closure_relative > 1e-7:
                raise RuntimeError(
                    f"{result_tag}: {projection} closure={closure_relative:.3g}"
                )

            chi2 = float(observed[projection]["chi2"])
            minimum_effective_entries = float(
                observed[projection]["minimum_effective_entries"]
            )
            rows.append(
                {
                    "projection": projection,
                    "bins": bins,
                    "chi2": chi2,
                    "ndf": bins - 1,
                    "chi2_ndf": chi2 / (bins - 1),
                    "minimum_effective_entries": minimum_effective_entries,
                    "expected_sum": expected_sum,
                    "expected_total": expected_total,
                    "closure_relative": closure_relative,
                    "maximum_equal_probability_deviation_relative":
                        maximum_equal_probability_deviation,
                    "result_tag": result_tag,
                }
            )

        for projection in PROJECTIONS:
            for suffix in ("pdf", "png"):
                plot = result_dir / "plots" / f"projection_{projection}.{suffix}"
                if not plot.is_file() or plot.stat().st_size == 0:
                    raise RuntimeError(f"Missing projection plot: {plot}")

    if len(symmetrization_modes) != 1 or len(sumw2_modes) != 1:
        raise RuntimeError("Worker projection conventions are inconsistent")
    symmetrization_mode = next(iter(symmetrization_modes))
    sumw2_mode = next(iter(sumw2_modes))

    output_dir.mkdir(parents=False)
    csv_path = output_dir / "chi2_binning_scan.csv"
    fieldnames = list(rows[0])
    with csv_path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    figure, axes = plt.subplots(2, 1, figsize=(7.2, 8.0), sharex=True)
    for projection, color, marker in zip(PROJECTIONS, COLORS, MARKERS):
        selected = [row for row in rows if row["projection"] == projection]
        axes[0].plot(
            [row["bins"] for row in selected],
            [row["chi2_ndf"] for row in selected],
            marker=marker,
            color=color,
            linewidth=1.5,
            markersize=5,
            label=projection,
        )
        axes[1].plot(
            [row["bins"] for row in selected],
            [row["minimum_effective_entries"] for row in selected],
            marker=marker,
            color=color,
            linewidth=1.5,
            markersize=5,
            label=projection,
        )

    axes[0].axhline(1.0, color="black", linestyle="--", linewidth=1.0)
    axes[0].set_ylabel(r"$\chi^2$/ndf")
    title_prefix = (
        "Pair-symmetrized " if symmetrization_mode == "enabled" else ""
    )
    axes[0].set_title(
        f"{title_prefix}adaptive-integral projection chi2 binning scan"
    )
    axes[0].legend(frameon=False, ncol=2)
    axes[1].set_xlabel("Number of nominal-PDF quantile bins")
    axes[1].set_ylabel("Minimum effective entries per bin")
    for axis in axes:
        axis.grid(alpha=0.25)
        axis.set_xlim(8, 62)
        axis.set_xticks(BIN_COUNTS)

    figure.tight_layout()
    figure.savefig(output_dir / "chi2_binning_scan.pdf")
    figure.savefig(output_dir / "chi2_binning_scan.png", dpi=180)
    plt.close(figure)

    metadata_path = output_dir / "run.metadata.txt"
    metadata_path.write_text(
        "\n".join(
            [
                "status=complete",
                f"campaign_tag={args.campaign_tag}",
                "projection_bins=10,15,20,25,30,35,40,45,50,55,60",
                "worker_results=11",
                "projections=Jpsi_mass1,Jpsi_mass2,Jpsi_ctau1,Jpsi_ctau2",
                f"projection_pair_symmetrization={symmetrization_mode}",
                f"projection_sumw2={sumw2_mode}",
                "projection_integration=adaptive Gauss-Kronrod interval integral normalized by the corresponding full-range integral",
                "integration_tolerances=absolute 1e-10; relative 1e-9; maximum 100000 subintervals",
                "validation=all worker metadata complete; four numeric projections and eight plot files per binning; expected-yield closure relative <= 1e-7",
                f"generated_at={datetime.now().astimezone().isoformat(timespec='seconds')}",
                f"sha256 Fit_Check.cpp={sha256(analysis_dir / 'Fit_Check.cpp')}",
                f"sha256 chi2_binning_scan.csv={sha256(csv_path)}",
                f"sha256 chi2_binning_scan.pdf={sha256(output_dir / 'chi2_binning_scan.pdf')}",
                f"sha256 chi2_binning_scan.png={sha256(output_dir / 'chi2_binning_scan.png')}",
                "",
            ]
        )
    )
    print(output_dir)


if __name__ == "__main__":
    main()
