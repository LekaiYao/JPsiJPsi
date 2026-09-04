#!/usr/bin/env python3

import argparse
import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[5]
MACRO = HERE / "single_jpsi_acceptance_closure.cpp"
SPS_PREFIX = Path("/eos/user/c/chensh/JPsiJPsi/GEN_nofilter/SPS/CMSSW_10_2_5/src/4mu_acc/direct")
DPS_PREFIX = Path("/eos/user/c/chensh/JPsiJPsi/GEN_nofilter/DPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc/direct")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def collect_inputs(prefix, stem):
    files = [(prefix / f"{stem}_{index}.root").resolve() for index in range(1, 11)]
    if len(set(files)) != len(files):
        raise RuntimeError(f"duplicate input paths in {stem}")
    metadata = []
    seen_hashes = {}
    for path in files:
        if not path.is_file():
            raise RuntimeError(f"missing input: {path}")
        digest = sha256(path)
        if digest in seen_hashes:
            raise RuntimeError(f"duplicate input content: {seen_hashes[digest]} and {path}")
        seen_hashes[digest] = path
        stat = path.stat()
        metadata.append({"path": str(path), "size": stat.st_size, "sha256": digest})
    return files, metadata


def root_string(value):
    return str(value).replace("\\", "\\\\").replace('"', '\\"')


def read_summary(path):
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 1:
        raise RuntimeError(f"expected one summary row in {path}")
    integers = {
        "sps_entries", "sps_jpsi_objects", "sps_single_fiducial", "sps_single_accepted",
        "dps_entries", "dps_jpsi_objects", "dps_single_fiducial", "dps_single_accepted",
        "unsupported_bins", "sps_empty_jpsi_event", "dps_empty_jpsi_event",
        "sps_inconsistent_jpsi_vectors", "dps_inconsistent_jpsi_vectors",
        "sps_malformed_muon", "dps_malformed_muon",
    }
    return {key: (int(value) if key in integers else float(value)) for key, value in rows[0].items()}


def command_output(command):
    return subprocess.run(command, text=True, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, check=True).stdout.strip()


def main():
    parser = argparse.ArgumentParser(
        description="Build single-J/psi SPS/DPS acceptance maps and close DPS with the SPS map.")
    parser.add_argument("--tag", required=True, help="Unique output tag; existing output is never overwritten.")
    parser.add_argument("--root", default="root")
    parser.add_argument("--output-root", type=Path, default=HERE / "closure_results")
    args = parser.parse_args()

    output = args.output_root.resolve() / args.tag
    if output.exists():
        raise RuntimeError(f"output already exists, refusing to overwrite: {output}")
    output.mkdir(parents=True)

    sps_files, sps_metadata = collect_inputs(SPS_PREFIX, "SPS_2016_JJ")
    dps_files, dps_metadata = collect_inputs(DPS_PREFIX, "DPS_2016_JJ")
    sps_list = output / "sps_gen_inputs.list"
    dps_list = output / "dps_gen_inputs.list"
    sps_list.write_text("".join(f"{path}\n" for path in sps_files))
    dps_list.write_text("".join(f"{path}\n" for path in dps_files))

    call = (f'{MACRO}("{root_string(output)}","{root_string(sps_list)}",'
            f'"{root_string(dps_list)}")')
    process = subprocess.run([args.root, "-l", "-b", "-q", call], text=True,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (output / "root.log").write_text(process.stdout)
    if process.returncode != 0:
        raise RuntimeError(f"ROOT closure failed with exit {process.returncode}; see {output / 'root.log'}")

    result = read_summary(output / "summary.csv")
    metadata = {
        "status": "complete",
        "tag": args.tag,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "repo_head": command_output(["git", "-C", str(REPO), "rev-parse", "HEAD"]),
        "root_version": command_output(["root-config", "--version"]),
        "source": str(MACRO),
        "source_sha256": sha256(MACRO),
        "runner_sha256": sha256(Path(__file__).resolve()),
        "selection": {
            "unit": "every individual object in the GENjpsi vectors; GEN_pair_id is not read",
            "fiducial": "10 <= pT < 40 GeV and -2 <= y < 2; no partner-J/psi or m(JJ) requirement",
            "accepted": "both decay muons have pT > 3.5 GeV and abs(eta) < 2.4",
            "dps_event_mixing": False,
        },
        "method": {
            "nominal_map": "SPS single-J/psi acceptance in 19 pT x 10 y bins",
            "closure_samples": "SPS self-closure and DPS closure both use individual J/psi from the original events",
            "closure_residual": "(N_fiducial - sum_accepted(1/A_SPS_bin)) / sum_accepted(1/A_SPS_bin)",
            "sps_self_closure_note": "The SPS map and SPS closure use the same bin counts, so this is a construction check rather than an independent validation.",
            "statistical_precision": "delta-method binomial propagation from DPS accepted counts and SPS map counts",
        },
        "inputs": {"sps": sps_metadata, "dps": dps_metadata},
        "result": result,
    }
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    artifacts = sorted(path for path in output.iterdir()
                       if path.is_file() and path.name != "checksums.sha256")
    (output / "checksums.sha256").write_text(
        "".join(f"{sha256(path)}  {path.name}\n" for path in artifacts))

    print(f"output={output}")
    print(f"sps_single_fiducial={result['sps_single_fiducial']}")
    print(f"dps_single_fiducial={result['dps_single_fiducial']}")
    print(f"sps_raw_acceptance={result['sps_raw_acceptance']:.10g}")
    print(f"dps_raw_acceptance={result['dps_raw_acceptance']:.10g}")
    print(f"sps_corrected_single={result['sps_corrected_single']:.10g}")
    print(f"sps_closure={result['sps_closure']:.10g}")
    print(f"dps_closure={result['dps_closure']:.10g}")
    print(f"dps_closure_stat_total={result['dps_closure_stat_total']:.10g}")


if __name__ == "__main__":
    main()
