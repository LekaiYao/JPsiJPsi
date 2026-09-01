#!/usr/bin/env python3
"""Create a tagged postprocess input from a completed f_DPS summary."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


def read_kv(path: pathlib.Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key] = value
    return values


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("usage: prepare_inputs.py SYSTEMATICS_TAG POSTPROCESS_TAG")
    systematics_tag, postprocess_tag = sys.argv[1:]
    repo = pathlib.Path(__file__).resolve().parents[3]
    systematics = repo / "Data_driven/results" / systematics_tag
    output = repo / "Data_driven/results" / postprocess_tag
    if output.exists():
        raise RuntimeError(f"refusing to overwrite {output}")
    summary = read_kv(systematics / "summary.txt")
    if summary.get("status") != "global_nominal_stat_and_systematics_complete_pending_user_confirmation":
        raise RuntimeError("systematics summary is incomplete or has an unexpected status")

    template_path = repo / "Data_driven/postprocess/inputs.json"
    pair_path = repo / "Data_driven/postprocess/current_pair_cross_section.json"
    inputs = json.loads(template_path.read_text(encoding="utf-8"))
    inputs["pair_cross_section"] = json.loads(pair_path.read_text(encoding="utf-8"))
    fraction = inputs["data_driven_fraction"]
    fraction.update(
        {
            "central": float(summary["nominal_f_dps"]),
            "stat": float(summary["statistical_uncertainty"]),
            "method_syst_up": float(summary["combined_method_systematic_up"]),
            "method_syst_down": float(summary["combined_method_systematic_down"]),
            "status": "candidate_pending_user_confirmation",
            "source": str((systematics / "summary.txt").relative_to(repo)),
            "bootstrap_source": str((systematics / "metadata.txt").relative_to(repo)),
            "method_groups": {
                "single_j_category": {
                    "up": float(summary["category_group_up"]),
                    "down": float(summary["category_group_down"]),
                },
                "cr_definition_and_purity": {
                    "up": float(summary["cr_group_up"]),
                    "down": float(summary["cr_group_down"]),
                },
                "pure_dps_full_workflow_closure": {
                    "up": float(summary["dps_mc_full_workflow_closure_symmetric"]),
                    "down": float(summary["dps_mc_full_workflow_closure_symmetric"]),
                },
            },
        }
    )
    inputs["producer"].update(
        {
            "git_head": subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=repo, text=True
            ).strip(),
            "runtime_for_frozen_analysis": "ROOT 6.40.02",
            "systematics_tag": systematics_tag,
            "postprocess_tag": postprocess_tag,
            "status": "candidate_pending_user_confirmation",
        }
    )
    inputs["cr_sps_leakage"].update(
        {
            "envelope_up": float(summary["cr_group_up"]),
            "envelope_down": float(summary["cr_group_down"]),
        }
    )
    output.mkdir(parents=True)
    (output / "inputs.json").write_text(
        json.dumps(inputs, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(output / "inputs.json")


if __name__ == "__main__":
    main()
