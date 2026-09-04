#!/usr/bin/env python3

import argparse
import csv
import hashlib
import json
import math
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[5]
BUILDER = HERE / "build_efficiency_maps_23style.cpp"
SPS_CLOSURE = HERE / "count.cpp"
DPS_CLOSURE = (
    REPO / "SKIM_tightfilter/DPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc_eff/count.cpp"
)
REFERENCE_ROOT = REPO / "tests/acceptance_efficiency_recompute_20260831/results"
REFERENCE_SPS_RAW = REFERENCE_ROOT / "sps_efficiency/raw_efficiency.txt"
REFERENCE_DPS_RAW = REFERENCE_ROOT / "dps_efficiency_dedup60/raw_efficiency.txt"
REFERENCE_SPS_INPUTS = REFERENCE_ROOT / "sps_efficiency/resolved_inputs.txt"
REFERENCE_DPS_INPUTS = REFERENCE_ROOT / "dps_efficiency_dedup60/resolved_inputs.txt"
FROZEN_MIXED = (
    REPO / "Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/"
    "efficiency_sps0p8_dps0p2_dedup60_v1.txt"
)
BASELINE_ARTIFACT = (
    REPO / "GEN_nofilter/DPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc/"
    "closure_results/allcross_frozen_sps0p8_dedup60_v2_20260902"
)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def command_output(command):
    return subprocess.run(command, text=True, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, check=True).stdout.strip()


def root_string(value):
    return str(value).replace("\\", "\\\\").replace('"', '\\"')


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


def read_raw(path):
    lines = path.read_text().splitlines()
    header = [int(value) for value in lines[0].split()]
    n_pt, n_y, n_acc, n_jpsi = header[:4]
    pt_edges = [float(value) for value in lines[1].split()]
    y_edges = [float(value) for value in lines[2].split()]
    single = []
    for line in lines[3:3 + n_y]:
        values = [int(value) for value in line.split()]
        if len(values) != 4 * n_pt:
            raise RuntimeError(f"bad single-J row in {path}")
        single.append(values)
    event = []
    for line in lines[3 + n_y:3 + n_y + n_pt]:
        values = [int(value) for value in line.split()]
        if len(values) != 3 * n_pt:
            raise RuntimeError(f"bad event row in {path}")
        event.append(values)
    if len(pt_edges) != n_pt + 1 or len(y_edges) != n_y + 1:
        raise RuntimeError(f"bad edges in {path}")
    return {
        "header": header,
        "n_pt": n_pt,
        "n_y": n_y,
        "n_acc": n_acc,
        "n_jpsi": n_jpsi,
        "pt_edges": pt_edges,
        "y_edges": y_edges,
        "single": single,
        "event": event,
    }


def raw_count_content(table):
    return {
        key: table[key]
        for key in ("n_pt", "n_y", "n_acc", "n_jpsi", "pt_edges", "y_edges", "single", "event")
    }


def mix_tables(sps, dps, sps_fraction=0.8, acceptance_sps=0.23713,
               acceptance_dps=0.18933):
    if sps["n_pt"] != dps["n_pt"] or sps["n_y"] != dps["n_y"]:
        raise RuntimeError("SPS/DPS raw map dimensions differ")
    if sps["pt_edges"] != dps["pt_edges"] or sps["y_edges"] != dps["y_edges"]:
        raise RuntimeError("SPS/DPS raw map edges differ")
    weight_sps = sps_fraction * acceptance_sps / sps["n_acc"]
    weight_dps = (1 - sps_fraction) * acceptance_dps / dps["n_acc"]
    single = [[weight_sps * a + weight_dps * b for a, b in zip(row_sps, row_dps)]
              for row_sps, row_dps in zip(sps["single"], dps["single"])]
    event = [[weight_sps * a + weight_dps * b for a, b in zip(row_sps, row_dps)]
             for row_sps, row_dps in zip(sps["event"], dps["event"])]
    return {"weight_sps": weight_sps, "weight_dps": weight_dps,
            "single": single, "event": event}


def write_mixed(path, reference, mixed):
    with path.open("w") as stream:
        stream.write(f"{reference['n_pt']} {reference['n_y']}\n")
        stream.write(" ".join(format(value, "g") for value in reference["pt_edges"]) + " \n")
        stream.write(" ".join(format(value, "g") for value in reference["y_edges"]) + " \n")
        for row in mixed["single"]:
            values = []
            for index in range(reference["n_pt"]):
                values.extend((row[4 * index], row[4 * index + 3]))
            stream.write(" ".join(str(value) for value in values) + " \n")
        for row in mixed["event"]:
            values = []
            for index in range(reference["n_pt"]):
                values.extend((row[3 * index], row[3 * index + 2]))
            stream.write(" ".join(str(value) for value in values) + " \n")


def component_statistics(table):
    groups = {"single_vertex": [], "event_trigger": []}
    for iy, row in enumerate(table["single"]):
        for ip in range(table["n_pt"]):
            groups["single_vertex"].append({
                "numerator": row[4 * ip], "denominator": row[4 * ip + 3],
                "pt_low": table["pt_edges"][ip], "pt_high": table["pt_edges"][ip + 1],
                "y_low": table["y_edges"][iy], "y_high": table["y_edges"][iy + 1],
            })
    for ip2, row in enumerate(table["event"]):
        for ip1 in range(table["n_pt"]):
            groups["event_trigger"].append({
                "numerator": row[3 * ip1], "denominator": row[3 * ip1 + 2],
                "pt1_low": table["pt_edges"][ip1], "pt1_high": table["pt_edges"][ip1 + 1],
                "pt2_low": table["pt_edges"][ip2], "pt2_high": table["pt_edges"][ip2 + 1],
            })
    result = {}
    for name, values in groups.items():
        supported = [value for value in values
                     if value["numerator"] > 0 and value["denominator"] > 0]
        for value in supported:
            efficiency = value["numerator"] / value["denominator"]
            value["relative_binomial_error"] = math.sqrt(
                (1 - efficiency) / (efficiency * value["denominator"]))
        result[name] = {
            "bins": len(values),
            "zero_numerator_bins": sum(value["numerator"] == 0 for value in values),
            "zero_denominator_bins": sum(value["denominator"] == 0 for value in values),
            "minimum_positive_numerator": min(supported, key=lambda value: value["numerator"]),
            "minimum_positive_denominator": min(supported, key=lambda value: value["denominator"]),
            "maximum_relative_binomial_error": max(
                supported, key=lambda value: value["relative_binomial_error"]),
        }
    return result


def mixed_support(reference, mixed):
    single = []
    for row in mixed["single"]:
        for index in range(reference["n_pt"]):
            single.append((row[4 * index], row[4 * index + 3]))
    event = []
    for row in mixed["event"]:
        for index in range(reference["n_pt"]):
            event.append((row[3 * index], row[3 * index + 2]))
    return {
        "single_zero_numerator_bins": sum(numerator <= 0 for numerator, _ in single),
        "single_zero_denominator_bins": sum(denominator <= 0 for _, denominator in single),
        "event_zero_numerator_bins": sum(numerator <= 0 for numerator, _ in event),
        "event_zero_denominator_bins": sum(denominator <= 0 for _, denominator in event),
        "minimum_single_numerator": min(numerator for numerator, _ in single),
        "minimum_single_denominator": min(denominator for _, denominator in single),
        "minimum_event_numerator": min(numerator for numerator, _ in event),
        "minimum_event_denominator": min(denominator for _, denominator in event),
    }


def closure_result(text):
    count = parse_scalar(text, "nPassAcc")
    weight = parse_scalar(text, "totWeight")
    return {
        "count": count,
        "matched": parse_scalar(text, "nMatchTrg"),
        "weight": weight,
        "residual": (count - weight) / weight,
    }


def input_metadata(paths):
    result = []
    for path in paths:
        stat = path.stat()
        result.append({"path": str(path), "size": stat.st_size,
                       "mtime_ns": stat.st_mtime_ns})
    return result


def main():
    parser = argparse.ArgumentParser(description="Run predefined 23-style mixed-efficiency closure sensitivity.")
    parser.add_argument("--tag", required=True)
    parser.add_argument("--root", default="root")
    parser.add_argument("--output-root", type=Path, default=HERE / "closure_results")
    args = parser.parse_args()

    output = args.output_root.resolve() / args.tag
    if output.exists():
        raise RuntimeError(f"output already exists, refusing to overwrite: {output}")
    output.mkdir(parents=True)

    baseline_metadata = json.loads((BASELINE_ARTIFACT / "metadata.json").read_text())
    source_hashes = {
        "sps_efficiency": sha256(SPS_CLOSURE),
        "dps_efficiency": sha256(DPS_CLOSURE),
    }
    for key, digest in source_hashes.items():
        if digest != baseline_metadata["source_sha256"][key]:
            raise RuntimeError(f"current {key} closure source differs from baseline artifact")

    sps_inputs = [Path(line) for line in REFERENCE_SPS_INPUTS.read_text().splitlines() if line]
    dps_inputs = [Path(line) for line in REFERENCE_DPS_INPUTS.read_text().splitlines() if line]
    for path in sps_inputs + dps_inputs:
        if not path.is_file():
            raise RuntimeError(f"missing input: {path}")
    sps_list = output / "sps_inputs.list"
    dps_list = output / "dps_inputs.list"
    sps_list.write_text("".join(f"{path}\n" for path in sps_inputs))
    dps_list.write_text("".join(f"{path}\n" for path in dps_inputs))

    raw_paths = {
        "sps19": output / "raw_efficiency_sps_19.txt",
        "sps23": output / "raw_efficiency_sps_23.txt",
        "dps19": output / "raw_efficiency_dps_19.txt",
        "dps23": output / "raw_efficiency_dps_23.txt",
    }
    run_root(args.root, BUILDER, ["sps", sps_list, raw_paths["sps19"], raw_paths["sps23"]],
             output / "build_sps.log")
    run_root(args.root, BUILDER, ["dps", dps_list, raw_paths["dps19"], raw_paths["dps23"]],
             output / "build_dps.log")
    tables = {key: read_raw(path) for key, path in raw_paths.items()}
    if raw_count_content(tables["sps19"]) != raw_count_content(read_raw(REFERENCE_SPS_RAW)):
        raise RuntimeError("generated 19-bin SPS raw map does not reproduce reference counts")
    if raw_count_content(tables["dps19"]) != raw_count_content(read_raw(REFERENCE_DPS_RAW)):
        raise RuntimeError("generated 19-bin DPS raw map does not reproduce reference counts")
    expected_refined_pt = [10 + 0.5 * index for index in range(9)] + list(range(15, 25)) + [26, 28, 30, 35, 40]
    if tables["sps23"]["pt_edges"] != expected_refined_pt:
        raise RuntimeError("unexpected 23-style pT edges")

    mixed19 = mix_tables(tables["sps19"], tables["dps19"])
    mixed23 = mix_tables(tables["sps23"], tables["dps23"])
    mixed19_path = output / "efficiency_sps0p8_dps0p2_19.txt"
    mixed23_path = output / "efficiency_sps0p8_dps0p2_23style.txt"
    write_mixed(mixed19_path, tables["sps19"], mixed19)
    write_mixed(mixed23_path, tables["sps23"], mixed23)
    if sha256(mixed19_path) != sha256(FROZEN_MIXED):
        raise RuntimeError("generated 19-bin mixed map does not byte-reproduce frozen nominal")

    statistics = {
        key: component_statistics(tables[key])
        for key in ("sps19", "sps23", "dps19", "dps23")
    }
    support19 = mixed_support(tables["sps19"], mixed19)
    support23 = mixed_support(tables["sps23"], mixed23)
    if any(support23[key] for key in (
            "single_zero_numerator_bins", "single_zero_denominator_bins",
            "event_zero_numerator_bins", "event_zero_denominator_bins")):
        raise RuntimeError("23-style mixed map has unsupported efficiency bins; closure not run")

    sps_text = run_root(args.root, SPS_CLOSURE, [mixed23_path], output / "closure_sps_23style.log")
    dps_text = run_root(args.root, DPS_CLOSURE, [mixed23_path], output / "closure_dps_23style.log")
    refined_results = {
        "sps": closure_result(sps_text),
        "dps": closure_result(dps_text),
    }
    baseline_results = {
        "sps": baseline_metadata["components"]["sps_efficiency"]["total"],
        "dps": baseline_metadata["components"]["dps_efficiency"]["total"],
    }
    comparison = []
    for label, results in (("baseline19", baseline_results), ("refined23style", refined_results)):
        comparison.append({
            "map": label,
            "n_pt": 19 if label == "baseline19" else 23,
            "sps_efficiency_residual": results["sps"]["residual"],
            "dps_efficiency_residual": results["dps"]["residual"],
            "efficiency_systematic": max(abs(results["sps"]["residual"]),
                                         abs(results["dps"]["residual"])),
        })
    with (output / "comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(comparison[0]))
        writer.writeheader()
        writer.writerows(comparison)

    metadata = {
        "status": "complete",
        "tag": args.tag,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "root_version": command_output(["root-config", "--version"]),
        "repo_head": command_output(["git", "-C", str(REPO), "rev-parse", "HEAD"]),
        "method": (
            "Predefined 23-style mixed-efficiency sensitivity: only 10--14 GeV pT bins "
            "change from 1 to 0.5 GeV; y, selection, 0.8:0.2 count-level mixing and closure stay frozen."
        ),
        "mixing": {"sps_fraction": 0.8, "dps_fraction": 0.2,
                   "acceptance_sps": 0.23713, "acceptance_dps": 0.18933,
                   "coarse_weight_sps": mixed19["weight_sps"],
                   "coarse_weight_dps": mixed19["weight_dps"]},
        "regressions": {
            "sps19_raw": "exact_count_match",
            "dps19_raw": "exact_count_match",
            "mixed19": "byte_identical",
            "frozen_mixed_sha256": sha256(FROZEN_MIXED),
            "generated_mixed_sha256": sha256(mixed19_path),
        },
        "statistics": statistics,
        "mixed_support": {"baseline19": support19, "refined23style": support23},
        "baseline": baseline_results,
        "refined": refined_results,
        "difference": {
            "sps_residual_absolute": refined_results["sps"]["residual"] - baseline_results["sps"]["residual"],
            "dps_residual_absolute": refined_results["dps"]["residual"] - baseline_results["dps"]["residual"],
            "systematic_absolute": comparison[1]["efficiency_systematic"] - comparison[0]["efficiency_systematic"],
        },
        "source_sha256": {"builder": sha256(BUILDER), "runner": sha256(Path(__file__).resolve()),
                          **source_hashes},
        "input_lists": {"sps_sha256": sha256(sps_list), "dps_sha256": sha256(dps_list)},
        "inputs": {"sps": input_metadata(sps_inputs), "dps": input_metadata(dps_inputs)},
        "artifacts": {key: {"path": str(path), "sha256": sha256(path)}
                      for key, path in {**raw_paths, "mixed19": mixed19_path,
                                        "mixed23": mixed23_path}.items()},
    }
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    artifacts = sorted(path for path in output.iterdir()
                       if path.is_file() and path.name != "checksums.sha256")
    (output / "checksums.sha256").write_text(
        "".join(f"{sha256(path)}  {path.name}\n" for path in artifacts)
    )

    print(f"output={output}")
    for row in comparison:
        print(f"{row['map']}: sps={row['sps_efficiency_residual']:.10g} "
              f"dps={row['dps_efficiency_residual']:.10g} "
              f"systematic={row['efficiency_systematic']:.10g}")
    print(f"systematic_change={metadata['difference']['systematic_absolute']:.10g}")


if __name__ == "__main__":
    main()
