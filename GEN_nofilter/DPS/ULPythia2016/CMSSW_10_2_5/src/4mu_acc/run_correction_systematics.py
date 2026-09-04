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
REPO = next(path for path in (HERE, *HERE.parents) if (path / ".git").exists())
DEFAULT_MANIFEST = HERE / "closure_inputs.json"

VARIABLES = {
    "delta_y": [0, 0.5, 1, 1.5, 2, 2.5, 4],
    "delta_phi": [0, 0.3927, 0.7854, 1.1781, 1.5708, 1.9635, 2.3562, 2.7489, 3.1416],
    "evt_y": [0, 0.4, 0.8, 1.2, 1.6, 2],
    "evt_pt": [0, 5, 10, 15, 20, 25, 30, 35, 40, 80],
    "evt_mass": [7.5, 17.5, 27.5, 37.5, 47.5, 57.5, 67.5, 107.5],
}
TPL_ORDER = ("delta_y", "delta_phi", "evt_mass", "evt_y", "evt_pt")
COMPONENTS = ("sps_acceptance", "sps_efficiency", "dps_acceptance", "dps_efficiency")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fraction_token(value):
    return format(float(value), ".12g").replace("-", "m").replace(".", "p")


def verify_file(entry, check_size=False):
    path = Path(entry["path"])
    if not path.is_file():
        raise RuntimeError(f"missing frozen input: {path}")
    if check_size and path.stat().st_size != entry["size"]:
        raise RuntimeError(f"size mismatch for {path}")
    actual = sha256(path)
    if actual != entry["sha256"]:
        raise RuntimeError(f"SHA256 mismatch for {path}: {actual}")
    return path


def read_acceptance_table(path):
    lines = Path(path).read_text().splitlines()
    n_pt, n_y = (int(value) for value in lines[0].split()[:2])
    pt_edges = [float(value) for value in lines[1].split()]
    y_edges = [float(value) for value in lines[2].split()]
    if len(pt_edges) != n_pt + 1 or len(y_edges) != n_y + 1:
        raise RuntimeError(f"bad acceptance binning in {path}")
    rows = []
    for line in lines[3:3 + n_y]:
        values = [float(value) for value in line.split()]
        if len(values) != 2 * n_pt:
            raise RuntimeError(f"bad acceptance row in {path}")
        rows.append([(values[2 * index], values[2 * index + 1]) for index in range(n_pt)])
    if len(rows) != n_y:
        raise RuntimeError(f"missing acceptance rows in {path}")
    return {"n_pt": n_pt, "n_y": n_y, "pt_edges": pt_edges, "y_edges": y_edges, "rows": rows}


def build_mixed_acceptance(sps_path, dps_path, fractions, output_path):
    sps = read_acceptance_table(sps_path)
    dps = read_acceptance_table(dps_path)
    for key in ("n_pt", "n_y", "pt_edges", "y_edges"):
        if sps[key] != dps[key]:
            raise RuntimeError(f"SPS/DPS acceptance {key} mismatch")
    sps_den = sum(den for row in sps["rows"] for _, den in row)
    dps_den = sum(den for row in dps["rows"] for _, den in row)
    if sps_den <= 0 or dps_den <= 0:
        raise RuntimeError("non-positive SPS/DPS acceptance denominator normalization")
    sps_scale = float(fractions["sps_fraction"]) / sps_den
    dps_scale = float(fractions["dps_fraction"]) / dps_den
    mixed_rows = []
    for sps_row, dps_row in zip(sps["rows"], dps["rows"]):
        mixed_row = []
        for (sps_num, sps_bin_den), (dps_num, dps_bin_den) in zip(sps_row, dps_row):
            mixed_num = sps_scale * sps_num + dps_scale * dps_num
            mixed_den = sps_scale * sps_bin_den + dps_scale * dps_bin_den
            mixed_row.append((mixed_num, mixed_den))
        mixed_rows.append(mixed_row)
    lines = [f'{sps["n_pt"]} {sps["n_y"]}']
    lines.append(" ".join(f"{value:.17g}" for value in sps["pt_edges"]))
    lines.append(" ".join(f"{value:.17g}" for value in sps["y_edges"]))
    for row in mixed_rows:
        lines.append(" ".join(f"{value:.17g}" for pair in row for value in pair))
    output_path.write_text("\n".join(lines) + "\n")
    return {
        "sps_fraction": float(fractions["sps_fraction"]),
        "dps_fraction": float(fractions["dps_fraction"]),
        "sps_denominator_sum": sps_den,
        "dps_denominator_sum": dps_den,
        "sps_scale": sps_scale,
        "dps_scale": dps_scale,
        "definition": fractions["normalization"],
    }


def cxx_string(value):
    return str(value).replace("\\", "\\\\").replace('"', '\\"')


def run_root(root_executable, macro, arguments, log_path):
    call = str(macro) + "(" + ",".join(f'"{cxx_string(value)}"' for value in arguments) + ")"
    proc = subprocess.run(
        [root_executable, "-l", "-b", "-q", call],
        cwd=macro.parent,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    log_path.write_text(proc.stdout)
    if proc.returncode != 0:
        raise RuntimeError(f"ROOT closure failed ({macro}); see {log_path}")
    return proc.stdout


def parse_scalar(text, key):
    match = re.search(rf"(?:^|\n){re.escape(key)}=([-+0-9.eE]+)", text)
    if not match:
        raise RuntimeError(f"missing {key} in ROOT output")
    value = float(match.group(1))
    if not math.isfinite(value):
        raise RuntimeError(f"non-finite {key} in ROOT output")
    return value


def parse_residuals(text):
    result = {}
    undefined = {}
    for variable, edges in VARIABLES.items():
        match = re.search(rf"(?:^|\n){re.escape(variable)}:\s*\{{([^}}]*)\}}", text)
        if not match:
            raise RuntimeError(f"missing {variable} residuals in ROOT output")
        tokens = [token.strip() for token in match.group(1).split(",") if token.strip()]
        if len(tokens) != len(edges) - 1:
            raise RuntimeError(f"wrong bin count for {variable}: {len(tokens)}")
        values = []
        undefined[variable] = {}
        for index, token in enumerate(tokens):
            bracket = re.fullmatch(r"\[([0-9]+)\]", token)
            if bracket:
                values.append(None)
                undefined[variable][str(index)] = int(bracket.group(1))
            else:
                values.append(float(token))
        if not all(value is None or math.isfinite(value) for value in values):
            raise RuntimeError(f"non-finite residual in {variable}")
        result[variable] = values
    return result, {key: value for key, value in undefined.items() if value}


def parse_component(name, text):
    count_key = "nEvent" if name.endswith("acceptance") else "nPassAcc"
    count = parse_scalar(text, count_key)
    weight = parse_scalar(text, "totWeight")
    if weight <= 0:
        raise RuntimeError(f"non-positive total weight for {name}")
    residuals, undefined = parse_residuals(text)
    component = {
        "total": {
            "count": count,
            "weight": weight,
            "residual": (count - weight) / weight,
        },
        "bins": residuals,
        "undefined_weight_bins": undefined,
    }
    if name == "dps_acceptance":
        component["diagnostics"] = {
            key: parse_scalar(text, key)
            for key in ("pool1", "pool2", "nTried", "nPassAcc", "closure_jackknife_stat")
        }
        printed = parse_scalar(text, "closure")
        if not math.isclose(printed, component["total"]["residual"], rel_tol=1e-12, abs_tol=1e-12):
            raise RuntimeError("DPS acceptance printed closure disagrees with counts")
    return component


def combine_sys3(components):
    def combine(values):
        sps_inputs = (values["sps_acceptance"], values["sps_efficiency"])
        dps_inputs = (values["dps_acceptance"], values["dps_efficiency"])
        sps = (
            (1.0 + sps_inputs[0]) * (1.0 + sps_inputs[1]) - 1.0
            if all(value is not None for value in sps_inputs)
            else None
        )
        dps = (
            (1.0 + dps_inputs[0]) * (1.0 + dps_inputs[1]) - 1.0
            if all(value is not None for value in dps_inputs)
            else None
        )
        available = [("sps", sps), ("dps", dps)]
        available = [(name, value) for name, value in available if value is not None]
        if not available:
            raise RuntimeError("neither SPS nor DPS has a defined combined closure")
        envelope_model, envelope_value = max(
            available, key=lambda item: abs(item[1])
        )
        return {
            "sps_combined": sps,
            "dps_combined": dps,
            "sys3": abs(envelope_value),
            "envelope_model": envelope_model,
            "available_models": "+".join(name for name, _ in available),
        }

    total_values = {name: components[name]["total"]["residual"] for name in COMPONENTS}
    result = {"total": {**total_values, **combine(total_values)}, "bins": {}}
    for variable in TPL_ORDER:
        rows = []
        for index in range(len(VARIABLES[variable]) - 1):
            values = {name: components[name]["bins"][variable][index] for name in COMPONENTS}
            rows.append({**values, **combine(values)})
        result["bins"][variable] = rows
    return result


def build_acceptance_systematic(components):
    def envelope(sps, dps):
        available = [("sps", sps), ("dps", dps)]
        available = [(name, value) for name, value in available if value is not None]
        if not available:
            raise RuntimeError("neither SPS nor DPS has a defined acceptance closure")
        return {
            "sps_acceptance_residual": sps,
            "dps_acceptance_residual": dps,
            "acceptance_systematic": max(abs(value) for _, value in available),
            "available_models": "+".join(name for name, _ in available),
        }

    result = {
        "total": envelope(
            components["sps_acceptance"]["total"]["residual"],
            components["dps_acceptance"]["total"]["residual"]),
        "bins": {},
    }
    for variable in TPL_ORDER:
        result["bins"][variable] = [
            envelope(sps, dps)
            for sps, dps in zip(
                components["sps_acceptance"]["bins"][variable],
                components["dps_acceptance"]["bins"][variable])
        ]
    return result


def write_component_csv(path, components):
    fields = ["component", "scope", "variable", "bin", "low", "high", "residual"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for name in COMPONENTS:
            writer.writerow({"component": name, "scope": "total", "residual": f'{components[name]["total"]["residual"]:.17g}'})
            for variable in TPL_ORDER:
                edges = VARIABLES[variable]
                for index, residual in enumerate(components[name]["bins"][variable]):
                    writer.writerow({
                        "component": name,
                        "scope": "bin",
                        "variable": variable,
                        "bin": index,
                        "low": edges[index],
                        "high": edges[index + 1],
                        "residual": "" if residual is None else f"{residual:.17g}",
                    })


def write_sys3_csv(path, sys3):
    fields = ["scope", "variable", "bin", "low", "high", *COMPONENTS, "sps_combined", "dps_combined", "sys3", "envelope_model", "available_models"]
    def formatted(values):
        return {
            key: (
                value
                if key in {"envelope_model", "available_models"}
                else ("" if value is None else f"{value:.17g}")
            )
            for key, value in values.items()
        }
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerow({"scope": "total", **formatted(sys3["total"])})
        for variable in TPL_ORDER:
            edges = VARIABLES[variable]
            for index, values in enumerate(sys3["bins"][variable]):
                writer.writerow({
                    "scope": "bin",
                    "variable": variable,
                    "bin": index,
                    "low": edges[index],
                    "high": edges[index + 1],
                    **formatted(values),
                })


def write_acceptance_systematic_csv(path, systematic):
    fields = ["scope", "variable", "bin", "low", "high", "sps_acceptance_residual", "dps_acceptance_residual", "acceptance_systematic", "available_models"]
    def formatted(values):
        return {
            key: ("" if value is None else f"{value:.17g}") if key != "available_models" else value
            for key, value in values.items()
        }
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerow({"scope": "total", **formatted(systematic["total"])})
        for variable in TPL_ORDER:
            edges = VARIABLES[variable]
            for index, values in enumerate(systematic["bins"][variable]):
                writer.writerow({
                    "scope": "bin",
                    "variable": variable,
                    "bin": index,
                    "low": edges[index],
                    "high": edges[index + 1],
                    **formatted(values),
                })


def write_cpp(path, sys3):
    lines = ["const vector<vector<double>> sys3 = {"]
    for variable in TPL_ORDER:
        values = ", ".join(f'{row["sys3"]:.10g}' for row in sys3["bins"][variable])
        lines.append(f"    {{{values}}}, // {variable}")
    lines.append("};")
    path.write_text("\n".join(lines) + "\n")


def command_version(command):
    proc = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True)
    return proc.stdout.strip()


def main():
    parser = argparse.ArgumentParser(description="Run the four correction closure tests and build sys3 deterministically.")
    parser.add_argument("--tag", required=True, help="Unique output tag; existing directories are never overwritten.")
    parser.add_argument("--root", default="root", help="ROOT executable.")
    parser.add_argument("--output-root", type=Path, default=HERE / "closure_results")
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST,
                        help="Frozen input manifest (default: closure_inputs.json).")
    args = parser.parse_args()

    output = args.output_root.resolve() / args.tag
    if output.exists():
        raise RuntimeError(f"output already exists, refusing to overwrite: {output}")
    output.mkdir(parents=True)

    manifest_path = args.manifest.resolve()
    manifest = json.loads(manifest_path.read_text())
    sps_acceptance = verify_file(manifest["sps_acceptance_table"])
    dps_acceptance = verify_file(manifest["dps_acceptance_table"])
    efficiency = verify_file(manifest["mixed_efficiency_table"])
    if sps_acceptance.read_text().splitlines()[0].strip() != manifest["sps_acceptance_table"]["header"]:
        raise RuntimeError("frozen SPS acceptance header mismatch")
    dps_files = [verify_file(entry, check_size=True) for entry in manifest["dps_gen"]["files"]]
    input_list = output / "dps_gen_inputs.list"
    input_list.write_text("".join(f"{path}\n" for path in dps_files))
    fractions = manifest["acceptance_mixing"]
    mixed_acceptance = output / (
        f"acceptance_sps{fraction_token(fractions['sps_fraction'])}"
        f"_dps{fraction_token(fractions['dps_fraction'])}.txt")
    acceptance_mix = build_mixed_acceptance(
        sps_acceptance, dps_acceptance, fractions, mixed_acceptance)

    macros = {
        "sps_acceptance": REPO / "GEN_nofilter/SPS/CMSSW_10_2_5/src/4mu_acc/count.cpp",
        "dps_acceptance": HERE / "count.cpp",
        "sps_efficiency": REPO / "SKIM_tightfilter/SPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc_eff/count.cpp",
        "dps_efficiency": REPO / "SKIM_tightfilter/DPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc_eff/count.cpp",
    }
    calls = {
        "sps_acceptance": [mixed_acceptance],
        "dps_acceptance": [mixed_acceptance, input_list],
        "sps_efficiency": [efficiency],
        "dps_efficiency": [efficiency],
    }

    components = {}
    for name in COMPONENTS:
        text = run_root(args.root, macros[name], calls[name], output / f"{name}.log")
        components[name] = parse_component(name, text)

    acceptance_systematic = build_acceptance_systematic(components)
    sys3 = combine_sys3(components)
    write_component_csv(output / "closure_components.csv", components)
    write_acceptance_systematic_csv(output / "acceptance_systematic.csv", acceptance_systematic)
    write_sys3_csv(output / "sys3.csv", sys3)
    write_cpp(output / "sys3_cpp.txt", sys3)

    metadata = {
        "status": "complete",
        "tag": args.tag,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "root_version": command_version(["root-config", "--version"]),
        "repo_head": command_version(["git", "-C", str(REPO), "rev-parse", "HEAD"]),
        "manifest": str(manifest_path),
        "manifest_sha256": sha256(manifest_path),
        "source_sha256": {name: sha256(path) for name, path in macros.items()},
        "input_sha256": {
            "sps_acceptance_table": sha256(sps_acceptance),
            "dps_acceptance_table": sha256(dps_acceptance),
            "mixed_acceptance_table": sha256(mixed_acceptance),
            "mixed_efficiency_table": sha256(efficiency),
            "dps_gen": {str(path): sha256(path) for path in dps_files},
        },
        "combination_rule": "For SPS and DPS separately, preserve the signed closure residuals and calculate (1 + acceptance_residual) * (1 + efficiency_residual) - 1; sys3 is the larger absolute combined residual.",
        "acceptance_mixing": acceptance_mix,
        "variable_order": list(TPL_ORDER),
        "components": components,
        "acceptance_systematic": acceptance_systematic,
        "correction_systematics": sys3,
    }
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")

    artifacts = sorted(path for path in output.iterdir() if path.is_file() and path.name != "checksums.sha256")
    (output / "checksums.sha256").write_text("".join(f"{sha256(path)}  {path.name}\n" for path in artifacts))
    print(f"output={output}")
    print(f'acceptance_systematic={acceptance_systematic["total"]["acceptance_systematic"]:.10g}')
    print(f'total_sys3={sys3["total"]["sys3"]:.10g}')
    print(f'dps_acceptance_total={components["dps_acceptance"]["total"]["residual"]:.10g}')
    print(f'dps_acceptance_jackknife={components["dps_acceptance"]["diagnostics"]["closure_jackknife_stat"]:.10g}')


if __name__ == "__main__":
    main()
