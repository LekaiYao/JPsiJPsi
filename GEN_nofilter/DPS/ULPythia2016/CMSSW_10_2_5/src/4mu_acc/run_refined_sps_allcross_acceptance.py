#!/usr/bin/env python3

import argparse
import csv
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from run_single_jpsi_acceptance_closure import (
    DPS_PREFIX,
    REPO,
    SPS_PREFIX,
    collect_inputs,
    command_output,
    root_string,
    sha256,
)


HERE = Path(__file__).resolve().parent
SPS_HERE = REPO / "GEN_nofilter/SPS/CMSSW_10_2_5/src/4mu_acc"
BUILDER = SPS_HERE / "build_pairconditioned_acceptance_maps.cpp"
SPS_CLOSURE = SPS_HERE / "count.cpp"
DPS_CLOSURE = HERE / "count.cpp"
BASELINE_MAP = REPO / "Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/acceptance_sps_full10_v1.txt"
BASELINE_DIR = HERE / "closure_results/allcross_frozen_sps0p8_dedup60_v2_20260902"


def run_root(root, macro, arguments, log_path):
    call = f"{macro}(" + ",".join(f'"{root_string(value)}"' for value in arguments) + ")"
    process = subprocess.run([root, "-l", "-b", "-q", call], text=True,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    log_path.write_text(process.stdout)
    if process.returncode != 0:
        raise RuntimeError(f"ROOT failed with exit {process.returncode}; see {log_path}")
    return process.stdout


def parse_scalar(text, key):
    match = re.search(rf"^{re.escape(key)}=([^\s]+)$", text, re.MULTILINE)
    if not match:
        raise RuntimeError(f"missing {key} in ROOT output")
    return float(match.group(1))


def read_table(path):
    lines = path.read_text().splitlines()
    header = [int(value) for value in lines[0].split()]
    n_pt, n_y = header[:2]
    pt_edges = [float(value) for value in lines[1].split()]
    y_edges = [float(value) for value in lines[2].split()]
    rows = []
    for line in lines[3:3 + n_y]:
        values = line.split()
        if len(values) != 2 * n_pt:
            raise RuntimeError(f"bad map row in {path}")
        rows.append([int(float(value)) for value in values])
    if len(pt_edges) != n_pt + 1 or len(y_edges) != n_y + 1 or len(rows) != n_y:
        raise RuntimeError(f"bad map dimensions in {path}")
    return {"header": header, "pt_edges": pt_edges, "y_edges": y_edges, "rows": rows}


def baseline_components(path):
    result = {}
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            if row["scope"] == "total" and row["component"] in ("sps_acceptance", "dps_acceptance"):
                result[row["component"]] = float(row["residual"])
    if len(result) != 2:
        raise RuntimeError("missing baseline acceptance totals")
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Refine SPS pT bins from 10--14 GeV and rerun deterministic all-cross acceptance closure.")
    parser.add_argument("--tag", required=True)
    parser.add_argument("--root", default="root")
    parser.add_argument("--output-root", type=Path, default=HERE / "closure_results")
    args = parser.parse_args()

    output = args.output_root.resolve() / args.tag
    if output.exists():
        raise RuntimeError(f"output already exists, refusing to overwrite: {output}")
    output.mkdir(parents=True)

    baseline_metadata = json.loads((BASELINE_DIR / "metadata.json").read_text())
    current_sources = {
        "sps_acceptance": sha256(SPS_CLOSURE),
        "dps_acceptance": sha256(DPS_CLOSURE),
    }
    for key, digest in current_sources.items():
        if digest != baseline_metadata["source_sha256"][key]:
            raise RuntimeError(f"current {key} source differs from baseline artifact")

    sps_files, sps_metadata = collect_inputs(SPS_PREFIX, "SPS_2016_JJ")
    dps_files, dps_metadata = collect_inputs(DPS_PREFIX, "DPS_2016_JJ")
    sps_list = output / "sps_gen_inputs.list"
    dps_list = output / "dps_gen_inputs.list"
    sps_list.write_text("".join(f"{path}\n" for path in sps_files))
    dps_list.write_text("".join(f"{path}\n" for path in dps_files))

    coarse_map = output / "acceptance_sps_pairconditioned_coarse19x10.txt"
    refined_map = output / "acceptance_sps_pairconditioned_refined23x10.txt"
    builder_text = run_root(args.root, BUILDER, [coarse_map, refined_map, sps_list],
                            output / "build_maps.log")
    coarse = read_table(coarse_map)
    baseline_map = read_table(BASELINE_MAP)
    if coarse != baseline_map:
        raise RuntimeError("coarse map count regression does not reproduce frozen baseline map")
    refined = read_table(refined_map)
    expected_pt = [10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14,
                   15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 28, 30, 35, 40]
    if refined["pt_edges"] != expected_pt or refined["y_edges"] != coarse["y_edges"]:
        raise RuntimeError("refined map edges do not match the predefined 23x10 binning")

    sps_text = run_root(args.root, SPS_CLOSURE, [refined_map], output / "sps_acceptance.log")
    dps_text = run_root(args.root, DPS_CLOSURE, [refined_map, dps_list], output / "dps_acceptance.log")
    refined_result = {
        "sps_acceptance_residual": parse_scalar(sps_text, "nEvent") - parse_scalar(sps_text, "totWeight"),
        "dps_acceptance_residual": parse_scalar(dps_text, "closure"),
        "dps_jackknife_precision": parse_scalar(dps_text, "closure_jackknife_stat"),
        "sps_n_count": parse_scalar(sps_text, "nEvent"),
        "sps_n_weight": parse_scalar(sps_text, "totWeight"),
        "dps_n_count": parse_scalar(dps_text, "nEvent"),
        "dps_n_weight": parse_scalar(dps_text, "totWeight"),
    }
    refined_result["sps_acceptance_residual"] /= refined_result["sps_n_weight"]
    refined_result["acceptance_systematic"] = max(
        abs(refined_result["sps_acceptance_residual"]),
        abs(refined_result["dps_acceptance_residual"]))

    baseline = baseline_components(BASELINE_DIR / "closure_components.csv")
    baseline_jackknife = baseline_metadata["components"]["dps_acceptance"]["diagnostics"]["closure_jackknife_stat"]
    baseline_result = {
        "sps_acceptance_residual": baseline["sps_acceptance"],
        "dps_acceptance_residual": baseline["dps_acceptance"],
        "dps_jackknife_precision": baseline_jackknife,
        "acceptance_systematic": max(abs(baseline["sps_acceptance"]), abs(baseline["dps_acceptance"])),
    }

    fields = ["map", "n_pt", "n_y", "sps_acceptance_residual", "dps_acceptance_residual",
              "dps_jackknife_precision", "acceptance_systematic"]
    with (output / "comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerow({"map": "baseline_coarse", "n_pt": 19, "n_y": 10, **baseline_result})
        writer.writerow({"map": "refined_pt10to14", "n_pt": 23, "n_y": 10,
                         **{key: refined_result[key] for key in fields[3:]}})

    metadata = {
        "status": "complete",
        "tag": args.tag,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "root_version": command_output(["root-config", "--version"]),
        "repo_head": command_output(["git", "-C", str(REPO), "rev-parse", "HEAD"]),
        "method": "Frozen SPS pair-conditioned map plus deterministic DPS all-cross-slot closure; only pT binning changes.",
        "selection": {
            "sps_map": "m(JJ)>7.5 GeV; both Jpsi fiducial; each map numerator applies its own muon acceptance",
            "closure": "SPS original pairs and DPS different-source all-cross-slot pairs; both Jpsi fiducial; m(JJ)>7.5 GeV; pair product weight",
        },
        "binning": {"coarse_pt_edges": coarse["pt_edges"], "refined_pt_edges": refined["pt_edges"],
                    "y_edges": refined["y_edges"]},
        "coarse_map_regression": {
            "status": "exact_count_match",
            "frozen_map": str(BASELINE_MAP),
            "frozen_map_sha256": sha256(BASELINE_MAP),
            "generated_map_sha256": sha256(coarse_map),
            "builder_n_event": parse_scalar(builder_text, "nEvent"),
            "builder_n_accepted_event": parse_scalar(builder_text, "nAcceptedEvent"),
        },
        "source_sha256": {"builder": sha256(BUILDER), **current_sources,
                          "runner": sha256(Path(__file__).resolve())},
        "inputs": {"sps": sps_metadata, "dps": dps_metadata},
        "baseline_artifact": str(BASELINE_DIR),
        "baseline": baseline_result,
        "refined": refined_result,
        "difference": {
            "dps_residual_absolute": refined_result["dps_acceptance_residual"] - baseline_result["dps_acceptance_residual"],
            "dps_residual_relative": refined_result["dps_acceptance_residual"] / baseline_result["dps_acceptance_residual"] - 1,
            "systematic_absolute": refined_result["acceptance_systematic"] - baseline_result["acceptance_systematic"],
        },
    }
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    artifacts = sorted(path for path in output.iterdir()
                       if path.is_file() and path.name != "checksums.sha256")
    (output / "checksums.sha256").write_text(
        "".join(f"{sha256(path)}  {path.name}\n" for path in artifacts))

    print(f"output={output}")
    print(f"baseline_sps_closure={baseline_result['sps_acceptance_residual']:.10g}")
    print(f"baseline_dps_closure={baseline_result['dps_acceptance_residual']:.10g}")
    print(f"refined_sps_closure={refined_result['sps_acceptance_residual']:.10g}")
    print(f"refined_dps_closure={refined_result['dps_acceptance_residual']:.10g}")
    print(f"refined_dps_jackknife={refined_result['dps_jackknife_precision']:.10g}")
    print(f"dps_closure_change={metadata['difference']['dps_residual_absolute']:.10g}")


if __name__ == "__main__":
    main()
