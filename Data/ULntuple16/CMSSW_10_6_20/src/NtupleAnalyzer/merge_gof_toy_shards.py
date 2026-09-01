#!/usr/bin/env python3

import argparse
import csv
import hashlib
import os
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
import re
import shutil
import tempfile


TOY_HEADER = [
    "toy_id", "seed", "n_raw", "sumw", "sumw2", "attempts", "accepted",
    "status", "covQual", "edm", "n_P_P_true", "n_P_P_fit",
    "n_P_P_error", "pull", "minNll", "chi2_mass1", "chi2_mass2",
    "chi2_ctau1", "chi2_ctau2",
]
ATTEMPT_HEADER = [
    "toy_id", "seed", "attempt", "mode", "status", "covQual", "edm",
    "minNll", "n_P_P_fit", "n_P_P_error", "near_boundary_count",
    "near_boundary_parameters", "exception",
]


def read_key_values(path):
    values = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key] = value
    return values


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def expected_seed(base_seed, global_toy_id):
    mask = (1 << 64) - 1
    mixed = (base_seed + 104729 * (global_toy_id + 1)) & mask
    mixed ^= mixed >> 16
    mixed = (mixed * 0x7FEB352D) & mask
    mixed ^= mixed >> 15
    seed = mixed % 900000000
    return seed if seed else 1


def normalized_configuration(path):
    lines = []
    for line in path.read_text().splitlines():
        lines.append("base_seed=<job-specific>" if line.startswith("base_seed=") else line)
    return "\n".join(lines)


def read_csv(path, expected_header):
    with path.open(newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        if header != expected_header:
            raise RuntimeError(f"Unexpected CSV header in {path}")
        return [row for row in reader if row]


def write_csv(path, header, rows):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Merge validated GOF toy shards into one global-ID campaign")
    parser.add_argument("campaign_tag")
    parser.add_argument("--jobs", type=int, default=100)
    parser.add_argument("--toys-per-job", type=int, default=50)
    parser.add_argument("--global-base-seed", type=int, default=20260901)
    args = parser.parse_args()

    if not re.fullmatch(r"[A-Za-z0-9_-]+", args.campaign_tag):
        raise RuntimeError("Invalid campaign tag")
    if args.jobs <= 0 or args.toys_per_job <= 0:
        raise RuntimeError("Job and toy counts must be positive")

    analysis_dir = Path(__file__).resolve().parent
    results_root = analysis_dir / "fit_results"
    campaign_dir = results_root / args.campaign_tag
    if not campaign_dir.is_dir():
        raise RuntimeError(f"Missing campaign directory: {campaign_dir}")

    output_names = [
        "toy_results.csv", "fit_attempts.csv", "configuration.txt",
        "projection_bins.csv", "observed_chi2.csv", "shard_manifest.csv",
        "merge_metadata.txt",
    ]
    existing = [name for name in output_names if (campaign_dir / name).exists()]
    if existing:
        raise RuntimeError("Refusing to overwrite merged outputs: " + ",".join(existing))

    expected_dirs = [
        results_root / f"{args.campaign_tag}_job{job_index:03d}"
        for job_index in range(args.jobs)
    ]
    actual_dirs = sorted(results_root.glob(f"{args.campaign_tag}_job*"))
    if set(actual_dirs) != set(expected_dirs):
        raise RuntimeError("Shard directory set does not match the expected job indices")

    merged_toys = []
    merged_attempts = []
    manifest_rows = []
    seen_global_ids = set()
    seen_seeds = set()
    attempts_by_toy = defaultdict(list)
    truth_values = set()
    result_attempt_counts = Counter()
    fit_modes = Counter()
    reference_observed = None
    reference_bins = None
    reference_configuration = None

    for job_index, shard_dir in enumerate(expected_dirs):
        metadata_path = shard_dir / "job.metadata.txt"
        metadata = read_key_values(metadata_path)
        global_start = job_index * args.toys_per_job
        global_end = global_start + args.toys_per_job - 1
        job_base_seed = args.global_base_seed + 104729 * global_start
        required_metadata = {
            "status": "complete",
            "campaign_tag": args.campaign_tag,
            "job_index": str(job_index),
            "global_toy_start": str(global_start),
            "global_toy_end": str(global_end),
            "local_toy_start": "0",
            "local_toy_end": str(args.toys_per_job - 1),
            "toys_requested": str(args.toys_per_job),
            "toys_accepted": str(args.toys_per_job),
            "global_base_seed": str(args.global_base_seed),
            "job_base_seed": str(job_base_seed),
            "toy0_diagnostics": "disabled",
            "root_version": "6.40.02",
        }
        for key, expected in required_metadata.items():
            if metadata.get(key) != expected:
                raise RuntimeError(
                    f"Metadata mismatch for {shard_dir.name}: {key}")
        if (shard_dir / "failure.txt").exists():
            raise RuntimeError(f"Failure artifact exists in {shard_dir.name}")
        if (shard_dir / "toy_0_fit_plots").exists() or \
           (shard_dir / "toy_0_fit.root").exists():
            raise RuntimeError(f"Toy-0 diagnostic artifact exists in {shard_dir.name}")

        observed_bytes = (shard_dir / "observed_chi2.csv").read_bytes()
        bins_bytes = (shard_dir / "projection_bins.csv").read_bytes()
        configuration = normalized_configuration(shard_dir / "configuration.txt")
        if job_index == 0:
            reference_observed = observed_bytes
            reference_bins = bins_bytes
            reference_configuration = configuration
        elif observed_bytes != reference_observed or bins_bytes != reference_bins:
            raise RuntimeError(f"Observed chi2/binning mismatch in {shard_dir.name}")
        elif configuration != reference_configuration:
            raise RuntimeError(f"Fit configuration mismatch in {shard_dir.name}")

        toy_rows = read_csv(shard_dir / "toy_results.csv", TOY_HEADER)
        if len(toy_rows) != args.toys_per_job:
            raise RuntimeError(f"Wrong toy row count in {shard_dir.name}")
        local_ids = set()
        for row in toy_rows:
            if len(row) != len(TOY_HEADER):
                raise RuntimeError(f"Malformed toy row in {shard_dir.name}")
            local_id = int(row[0])
            global_id = global_start + local_id
            if local_id < 0 or local_id >= args.toys_per_job or \
               local_id in local_ids or global_id in seen_global_ids:
                raise RuntimeError(f"Invalid or duplicate toy ID in {shard_dir.name}")
            local_ids.add(local_id)
            seen_global_ids.add(global_id)
            seed = int(row[1])
            if seed != expected_seed(args.global_base_seed, global_id):
                raise RuntimeError(f"Seed mismatch for global toy {global_id}")
            if seed in seen_seeds:
                raise RuntimeError(f"Duplicate seed for global toy {global_id}")
            seen_seeds.add(seed)
            if int(row[6]) != 1 or int(row[7]) != 0 or int(row[8]) != 3 or \
               float(row[9]) >= 0.01:
                raise RuntimeError(f"Fit-quality failure for global toy {global_id}")
            attempts = int(row[5])
            if attempts < 1 or attempts > 3:
                raise RuntimeError(f"Invalid attempt count for global toy {global_id}")
            result_attempt_counts[attempts] += 1
            truth_values.add(row[10])
            merged_row = row.copy()
            merged_row[0] = str(global_id)
            merged_toys.append(merged_row)
        if local_ids != set(range(args.toys_per_job)):
            raise RuntimeError(f"Local toy IDs are not contiguous in {shard_dir.name}")

        attempt_rows = read_csv(shard_dir / "fit_attempts.csv", ATTEMPT_HEADER)
        for row in attempt_rows:
            if len(row) != len(ATTEMPT_HEADER):
                raise RuntimeError(f"Malformed fit-attempt row in {shard_dir.name}")
            local_id = int(row[0])
            global_id = global_start + local_id
            if local_id < 0 or local_id >= args.toys_per_job:
                raise RuntimeError(f"Invalid attempt toy ID in {shard_dir.name}")
            if int(row[1]) != expected_seed(args.global_base_seed, global_id):
                raise RuntimeError(f"Attempt seed mismatch for global toy {global_id}")
            attempt_number = int(row[2])
            attempts_by_toy[global_id].append(attempt_number)
            fit_modes[row[3]] += 1
            merged_row = row.copy()
            merged_row[0] = str(global_id)
            merged_attempts.append(merged_row)

        manifest_rows.append([
            job_index, shard_dir.name, global_start, global_end, job_base_seed,
            sha256(metadata_path), sha256(shard_dir / "toy_results.csv"),
            sha256(shard_dir / "fit_attempts.csv"),
        ])

    expected_toys = args.jobs * args.toys_per_job
    if seen_global_ids != set(range(expected_toys)):
        raise RuntimeError("Merged global toy IDs are not contiguous")
    if len(seen_seeds) != expected_toys:
        raise RuntimeError("Merged toy seeds are not unique")
    if len(truth_values) != 1:
        raise RuntimeError("n_P_P truth differs between shards")
    for row in merged_toys:
        global_id = int(row[0])
        expected_attempts = int(row[5])
        if attempts_by_toy[global_id] != list(range(1, expected_attempts + 1)):
            raise RuntimeError(f"Attempt sequence mismatch for global toy {global_id}")

    merged_toys.sort(key=lambda row: int(row[0]))
    merged_attempts.sort(key=lambda row: (int(row[0]), int(row[2])))

    temporary_dir = Path(tempfile.mkdtemp(prefix=".merge-", dir=campaign_dir))
    try:
        write_csv(temporary_dir / "toy_results.csv", TOY_HEADER, merged_toys)
        write_csv(temporary_dir / "fit_attempts.csv", ATTEMPT_HEADER, merged_attempts)
        write_csv(
            temporary_dir / "shard_manifest.csv",
            ["job_index", "job_tag", "global_toy_start", "global_toy_end",
             "job_base_seed", "job_metadata_sha256", "toy_results_sha256",
             "fit_attempts_sha256"],
            manifest_rows)
        shutil.copyfile(expected_dirs[0] / "configuration.txt",
                        temporary_dir / "configuration.txt")
        shutil.copyfile(expected_dirs[0] / "projection_bins.csv",
                        temporary_dir / "projection_bins.csv")
        shutil.copyfile(expected_dirs[0] / "observed_chi2.csv",
                        temporary_dir / "observed_chi2.csv")

        metadata_lines = [
            "status=complete",
            f"generated_at={datetime.now().astimezone().isoformat()}",
            f"campaign_tag={args.campaign_tag}",
            f"jobs={args.jobs}",
            f"toys_per_job={args.toys_per_job}",
            f"toys_merged={expected_toys}",
            f"fit_attempt_rows={len(merged_attempts)}",
            f"global_toy_range=0--{expected_toys - 1}",
            f"unique_seeds={len(seen_seeds)}",
            f"global_base_seed={args.global_base_seed}",
            "seed_mapping=job_base_seed=global_base_seed+104729*global_toy_start",
            f"n_P_P_truth={next(iter(truth_values))}",
            "fit_gate=status=0,covQual=3,EDM<0.01,n_P_P error finite and positive",
            "toy0_diagnostics=disabled in all shards",
            f"attempts_1={result_attempt_counts[1]}",
            f"attempts_2={result_attempt_counts[2]}",
            f"attempts_3={result_attempt_counts[3]}",
        ]
        for mode, count in sorted(fit_modes.items()):
            metadata_lines.append(f"fit_mode_{mode}={count}")
        for name in output_names[:-1]:
            metadata_lines.append(f"sha256 {name}={sha256(temporary_dir / name)}")
        (temporary_dir / "merge_metadata.txt").write_text(
            "\n".join(metadata_lines) + "\n")

        for name in output_names:
            os.replace(temporary_dir / name, campaign_dir / name)
        temporary_dir.rmdir()
    except Exception:
        shutil.rmtree(temporary_dir, ignore_errors=True)
        raise

    print(
        f"Merged {expected_toys} toys and {len(merged_attempts)} fit attempts "
        f"into {campaign_dir}")


if __name__ == "__main__":
    main()
