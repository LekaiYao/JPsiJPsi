#!/usr/bin/env python3

import argparse
import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path

from run_refined_sps_allcross_acceptance import (
    BASELINE_DIR,
    BASELINE_MAP,
    BUILDER,
    DPS_CLOSURE,
    HERE,
    SPS_CLOSURE,
    parse_scalar,
    read_table,
    run_root,
)
from run_single_jpsi_acceptance_closure import (
    DPS_PREFIX,
    REPO,
    SPS_PREFIX,
    collect_inputs,
    command_output,
    sha256,
)


PREVIOUS_REFINED_DIR = (
    HERE / "closure_results/allcross_spsmap_pt10to14_0p5_v1_20260902"
)


def expected_edges():
    coarse = list(range(10, 25)) + [26, 28, 30, 35, 40]
    refined = [10 + 0.5 * index for index in range(9)] + list(range(15, 25)) + [26, 28, 30, 35, 40]
    quarter = [10 + 0.25 * index for index in range(17)] + list(range(15, 25)) + [26, 28, 30, 35, 40]
    piecewise = (
        [10 + 0.25 * index for index in range(17)]
        + [14 + 0.5 * index for index in range(1, 33)]
        + list(range(31, 41))
    )
    return {"coarse19x10": coarse, "refined23x10": refined,
            "quarter31x10": quarter, "piecewise58x10": piecewise}


def map_statistics(table):
    values = []
    for iy, row in enumerate(table["rows"]):
        for ip in range(len(table["pt_edges"]) - 1):
            accepted = row[2 * ip]
            fiducial = row[2 * ip + 1]
            if accepted <= 0 or fiducial <= 0 or accepted > fiducial:
                raise RuntimeError("invalid SPS acceptance map count")
            acceptance = accepted / fiducial
            relative_error = math.sqrt((1 - acceptance) / (acceptance * fiducial))
            values.append({
                "accepted": accepted,
                "fiducial": fiducial,
                "failed": fiducial - accepted,
                "relative_error": relative_error,
                "pt_low": table["pt_edges"][ip],
                "pt_high": table["pt_edges"][ip + 1],
                "y_low": table["y_edges"][iy],
                "y_high": table["y_edges"][iy + 1],
            })
    minimum_accepted = min(values, key=lambda value: value["accepted"])
    minimum_fiducial = min(values, key=lambda value: value["fiducial"])
    minimum_failed = min(values, key=lambda value: value["failed"])
    maximum_relative_error = max(values, key=lambda value: value["relative_error"])
    return {
        "minimum_accepted": minimum_accepted,
        "minimum_fiducial": minimum_fiducial,
        "minimum_failed": minimum_failed,
        "maximum_relative_binomial_error": maximum_relative_error,
    }


def run_closure(root, acceptance_map, dps_list, output, label):
    sps_text = run_root(root, SPS_CLOSURE, [acceptance_map],
                        output / f"sps_acceptance_{label}.log")
    dps_text = run_root(root, DPS_CLOSURE, [acceptance_map, dps_list],
                        output / f"dps_acceptance_{label}.log")
    sps_n_count = parse_scalar(sps_text, "nEvent")
    sps_n_weight = parse_scalar(sps_text, "totWeight")
    result = {
        "sps_acceptance_residual": (sps_n_count - sps_n_weight) / sps_n_weight,
        "dps_acceptance_residual": parse_scalar(dps_text, "closure"),
        "dps_jackknife_precision": parse_scalar(dps_text, "closure_jackknife_stat"),
        "sps_n_count": sps_n_count,
        "sps_n_weight": sps_n_weight,
        "dps_n_count": parse_scalar(dps_text, "nEvent"),
        "dps_n_weight": parse_scalar(dps_text, "totWeight"),
    }
    result["acceptance_systematic"] = max(
        abs(result["sps_acceptance_residual"]),
        abs(result["dps_acceptance_residual"]),
    )
    return result


def comparison_row(label, table, result, statistics):
    return {
        "map": label,
        "n_pt": len(table["pt_edges"]) - 1,
        "n_y": len(table["y_edges"]) - 1,
        "sps_acceptance_residual": result["sps_acceptance_residual"],
        "dps_acceptance_residual": result["dps_acceptance_residual"],
        "dps_jackknife_precision": result["dps_jackknife_precision"],
        "acceptance_systematic": result["acceptance_systematic"],
        "minimum_sps_accepted": statistics["minimum_accepted"]["accepted"],
        "minimum_sps_fiducial": statistics["minimum_fiducial"]["fiducial"],
        "maximum_sps_relative_binomial_error": statistics[
            "maximum_relative_binomial_error"]["relative_error"],
    }


def main():
    parser = argparse.ArgumentParser(
        description="Compare predefined 31x10 and 58x10 SPS maps with frozen all-cross closures."
    )
    parser.add_argument("--tag", required=True)
    parser.add_argument("--root", default="root")
    parser.add_argument("--output-root", type=Path, default=HERE / "closure_results")
    args = parser.parse_args()

    output = args.output_root.resolve() / args.tag
    if output.exists():
        raise RuntimeError(f"output already exists, refusing to overwrite: {output}")
    output.mkdir(parents=True)

    baseline_metadata = json.loads((BASELINE_DIR / "metadata.json").read_text())
    previous_metadata = json.loads((PREVIOUS_REFINED_DIR / "metadata.json").read_text())
    current_sources = {
        "sps_acceptance": sha256(SPS_CLOSURE),
        "dps_acceptance": sha256(DPS_CLOSURE),
    }
    for key, digest in current_sources.items():
        if digest != baseline_metadata["source_sha256"][key]:
            raise RuntimeError(f"current {key} source differs from frozen baseline artifact")

    sps_files, sps_metadata = collect_inputs(SPS_PREFIX, "SPS_2016_JJ")
    dps_files, dps_metadata = collect_inputs(DPS_PREFIX, "DPS_2016_JJ")
    sps_list = output / "sps_gen_inputs.list"
    dps_list = output / "dps_gen_inputs.list"
    sps_list.write_text("".join(f"{path}\n" for path in sps_files))
    dps_list.write_text("".join(f"{path}\n" for path in dps_files))

    map_paths = {
        "coarse19x10": output / "acceptance_sps_pairconditioned_coarse19x10.txt",
        "refined23x10": output / "acceptance_sps_pairconditioned_refined23x10.txt",
        "quarter31x10": output / "acceptance_sps_pairconditioned_quarter31x10.txt",
        "piecewise58x10": output / "acceptance_sps_pairconditioned_piecewise58x10.txt",
    }
    builder_text = run_root(
        args.root,
        BUILDER,
        [map_paths["coarse19x10"], map_paths["refined23x10"], sps_list,
         map_paths["quarter31x10"], map_paths["piecewise58x10"]],
        output / "build_maps.log",
    )

    tables = {key: read_table(path) for key, path in map_paths.items()}
    for key, edges in expected_edges().items():
        if tables[key]["pt_edges"] != edges:
            raise RuntimeError(f"unexpected pT edges for {key}")
    if len({tuple(table["y_edges"]) for table in tables.values()}) != 1:
        raise RuntimeError("y edges differ between comparison maps")
    if tables["coarse19x10"] != read_table(BASELINE_MAP):
        raise RuntimeError("generated 19x10 map does not reproduce frozen baseline")
    previous_refined_map = PREVIOUS_REFINED_DIR / "acceptance_sps_pairconditioned_refined23x10.txt"
    if tables["refined23x10"] != read_table(previous_refined_map):
        raise RuntimeError("generated 23x10 map does not reproduce previous artifact")

    results = {
        "coarse19x10": previous_metadata["baseline"],
        "refined23x10": previous_metadata["refined"],
        "quarter31x10": run_closure(
            args.root, map_paths["quarter31x10"], dps_list, output, "quarter31x10"),
        "piecewise58x10": run_closure(
            args.root, map_paths["piecewise58x10"], dps_list, output, "piecewise58x10"),
    }
    statistics = {key: map_statistics(table) for key, table in tables.items()}
    rows = [comparison_row(key, tables[key], results[key], statistics[key])
            for key in ("coarse19x10", "refined23x10", "quarter31x10", "piecewise58x10")]
    fields = list(rows[0])
    with (output / "comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    baseline = results["coarse19x10"]
    refined = results["refined23x10"]
    differences = {}
    for key in ("quarter31x10", "piecewise58x10"):
        result = results[key]
        differences[key] = {
            "dps_residual_minus_19x10": (
                result["dps_acceptance_residual"] - baseline["dps_acceptance_residual"]),
            "dps_residual_minus_23x10": (
                result["dps_acceptance_residual"] - refined["dps_acceptance_residual"]),
            "systematic_minus_19x10": (
                result["acceptance_systematic"] - baseline["acceptance_systematic"]),
            "systematic_minus_23x10": (
                result["acceptance_systematic"] - refined["acceptance_systematic"]),
        }

    metadata = {
        "status": "complete",
        "tag": args.tag,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "root_version": command_output(["root-config", "--version"]),
        "repo_head": command_output(["git", "-C", str(REPO), "rev-parse", "HEAD"]),
        "method": (
            "Original pair-conditioned SPS acceptance map and deterministic different-source "
            "DPS all-cross-slot closure; only predefined pT edges change; all y edges remain frozen."
        ),
        "selection": {
            "sps_map": (
                "m(JJ)>7.5 GeV; both pair-selected Jpsi fiducial; each numerator applies "
                "that Jpsi's own muon acceptance"
            ),
            "closure": (
                "SPS original pairs and DPS different-source all-cross-slot pairs; both Jpsi "
                "fiducial; m(JJ)>7.5 GeV; pair product weight"
            ),
        },
        "map_regressions": {
            "coarse19x10": {
                "status": "exact_count_match",
                "reference": str(BASELINE_MAP),
                "reference_sha256": sha256(BASELINE_MAP),
                "generated_sha256": sha256(map_paths["coarse19x10"]),
            },
            "refined23x10": {
                "status": "exact_count_match",
                "reference": str(previous_refined_map),
                "reference_sha256": sha256(previous_refined_map),
                "generated_sha256": sha256(map_paths["refined23x10"]),
            },
            "builder_n_event": parse_scalar(builder_text, "nEvent"),
            "builder_n_accepted_event": parse_scalar(builder_text, "nAcceptedEvent"),
        },
        "maps": {
            key: {
                "path": str(map_paths[key]),
                "sha256": sha256(map_paths[key]),
                "pt_edges": tables[key]["pt_edges"],
                "y_edges": tables[key]["y_edges"],
                "statistics": statistics[key],
                "closure": results[key],
            }
            for key in map_paths
        },
        "differences": differences,
        "existing_artifacts": {
            "coarse19x10": str(BASELINE_DIR),
            "refined23x10": str(PREVIOUS_REFINED_DIR),
        },
        "source_sha256": {
            "builder": sha256(BUILDER),
            "runner": sha256(Path(__file__).resolve()),
            **current_sources,
        },
        "inputs": {"sps": sps_metadata, "dps": dps_metadata},
    }
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    artifacts = sorted(path for path in output.iterdir()
                       if path.is_file() and path.name != "checksums.sha256")
    (output / "checksums.sha256").write_text(
        "".join(f"{sha256(path)}  {path.name}\n" for path in artifacts)
    )

    print(f"output={output}")
    for row in rows:
        print(
            f"{row['map']}: sps={row['sps_acceptance_residual']:.10g} "
            f"dps={row['dps_acceptance_residual']:.10g} "
            f"jackknife={row['dps_jackknife_precision']:.10g} "
            f"systematic={row['acceptance_systematic']:.10g}"
        )


if __name__ == "__main__":
    main()
