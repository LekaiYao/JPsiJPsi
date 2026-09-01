#!/usr/bin/env python3
"""Validate and summarize the fixed-500 full-workflow bootstrap."""

from __future__ import annotations

import csv
import math
import pathlib
import statistics
import sys


def read_kv(path: pathlib.Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key] = value
    return values


def quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    low = math.floor(position)
    high = math.ceil(position)
    if low == high:
        return ordered[low]
    fraction = position - low
    return ordered[low] * (1.0 - fraction) + ordered[high] * fraction


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("usage: summarize_bootstrap.py CENTRAL_TAG BOOTSTRAP_TAG")
    repo = pathlib.Path(__file__).resolve().parents[1]
    central_tag, bootstrap_tag = sys.argv[1:]
    central = repo / "Data_driven/results" / central_tag
    output = repo / "Data_driven/results" / bootstrap_tag
    summary_path = output / "summary.txt"
    if summary_path.exists():
        raise RuntimeError(f"refusing to overwrite {summary_path}")

    nominal = float(read_kv(central / "templates_2d/summary.txt")["f_dps"])
    submission = read_kv(output / "submission_metadata.txt")
    if submission.get("central_tag") != central_tag:
        raise RuntimeError("bootstrap submission central tag does not match")
    if submission.get("root_version") != "6.40.02":
        raise RuntimeError("bootstrap submission did not use ROOT 6.40.02")
    base_seed = int(submission["base_seed"])

    job_summaries = sorted(output.glob("job_*/job_summary.txt"))
    if len(job_summaries) != 100:
        raise RuntimeError(f"expected 100 job summaries, found {len(job_summaries)}")
    for path in job_summaries:
        values = read_kv(path)
        if values.get("state") != "complete" or values.get("successes") != "5":
            raise RuntimeError(f"incomplete bootstrap job: {path}")

    rows: list[dict[str, object]] = []
    failed: list[pathlib.Path] = []
    for status_path in sorted(output.glob("job_*/attempt_*/attempt_status.txt")):
        status = read_kv(status_path)
        if status.get("status") != "success":
            failed.append(status_path)
            continue
        metadata_path = status_path.parent / "metadata.txt"
        metadata = read_kv(metadata_path)
        if metadata.get("status") != "complete":
            raise RuntimeError(f"successful attempt has incomplete metadata: {metadata_path}")
        replica_id = int(status["attempt_id"])
        if int(metadata["replica_id"]) != replica_id:
            raise RuntimeError(f"replica ID mismatch: {metadata_path}")
        expected_seed = base_seed + 1_000_003 * replica_id
        if int(status["seed"]) != expected_seed:
            raise RuntimeError(f"seed mismatch: {status_path}")
        rows.append(
            {
                "replica_id": replica_id,
                "seed": expected_seed,
                "f_dps": float(metadata["primary_f_dps"]),
                "source_metadata": str(metadata_path.relative_to(repo)),
            }
        )

    rows.sort(key=lambda row: int(row["replica_id"]))
    expected_ids = {job * 10 + local for job in range(100) for local in range(5)}
    ids = [int(row["replica_id"]) for row in rows]
    if len(rows) != 500 or len(set(ids)) != 500 or set(ids) != expected_ids:
        raise RuntimeError("fixed-500 replica IDs are incomplete, duplicated, or replaced")
    if failed:
        raise RuntimeError(f"bootstrap contains {len(failed)} failed attempts")

    values = [float(row["f_dps"]) for row in rows]
    mean = statistics.mean(values)
    sd = statistics.stdev(values)
    mc_precision = sd / math.sqrt(2.0 * (len(values) - 1))
    with (output / "replicas.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    summary = [
        "status=complete",
        f"central_tag={central_tag}",
        f"bootstrap_tag={bootstrap_tag}",
        "successes=500",
        "failures=0",
        "attempts=500",
        "replacement_seeds=false",
        f"base_seed={base_seed}",
        f"nominal_f_dps={nominal:.12g}",
        f"mean={mean:.12g}",
        f"mean_shift_from_nominal={mean - nominal:.12g}",
        f"statistical_uncertainty={sd:.12g}",
        f"median={statistics.median(values):.12g}",
        f"q16={quantile(values, 0.16):.12g}",
        f"q84={quantile(values, 0.84):.12g}",
        f"minimum={min(values):.12g}",
        f"maximum={max(values):.12g}",
        f"sd_monte_carlo_precision_approx={mc_precision:.12g}",
        "root_version=6.40.02",
        "status_for_physics=pending_user_confirmation",
    ]
    summary_path.write_text("\n".join(summary) + "\n", encoding="utf-8")
    (output / "complete.marker").write_text("complete_pending_user_confirmation\n", encoding="utf-8")
    print("\n".join(summary))


if __name__ == "__main__":
    main()
