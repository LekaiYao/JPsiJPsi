#!/usr/bin/env python3
import csv
import hashlib
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
from datetime import datetime


REPO = Path("/eos/home-l/leyao/26JJ/JPsiJPsi")
DD = REPO / "Data_driven"
UPSTREAM_TAG = "newdata167_accmix23_effmix19_lumi36p684_root640_v1_20260902"
UPSTREAM = REPO / "Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/fit_results" / UPSTREAM_TAG
MODEL_SOURCE = UPSTREAM / "total/Model_4D_tot.root"
FIT_RESULT_SOURCE = UPSTREAM / "total/Fit_4D_tot_native_asymptotic_unseeded.root"
WEIGHT_DATA_SOURCE = UPSTREAM / "inputs/WeightData_newdata167_accmix23_effmix19_v1_20260902.root"
DATA_MANIFEST_SOURCE = REPO / "Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/data_ntuple_manifest_chensh_newdata_v1.list"
ACCEPTANCE_SOURCE = REPO / "GEN_nofilter/DPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc/closure_results/nominal_mixed23_sps0p8_dps0p2_avgacc_v1_20260902/acceptance_sps0p8_dps0p2.txt"
EFFICIENCY_SOURCE = REPO / "GEN_nofilter/DPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc/closure_results/inputs_mixed23_sps0p8_dps0p2_avgacc_v1_20260902/efficiency_sps0p8_dps0p2_avgacc_19x10.txt"
SELECTION_SOURCE = REPO / "Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/rephrase.cpp"
DPS_BASE = REPO / "Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/DPS_ntuple"

EXPECTED = {
    MODEL_SOURCE: "9987bf734f9308f4177b290828cff82c7cb470d035cca8bbdba337534af61238",
    FIT_RESULT_SOURCE: "d0a237af9042d663cce10a217b5a869963bf6c697887a5ef6acabfb93f3cf25a",
    WEIGHT_DATA_SOURCE: "8cd1abae2a896f5ebf82a805e62009e96a7f7e4102350c1f6e18c2568c1e8957",
    DATA_MANIFEST_SOURCE: "21c2830dc65e6be490f3ed262be30800b23cc92da61ea197305fed640fb84dee",
    ACCEPTANCE_SOURCE: "71cfd7a62f4bb2c1a4bff045b96bbe78498a9798d30f0f6bd38b6bc6d262c077",
    EFFICIENCY_SOURCE: "e3ced212b79d33bbb5bc218648867c2984781efcfb4ca8ac6299baed4f92a256",
}
PHI70 = 1.2217304763960306
PHI90 = 1.5707963267948966
PHI110 = 1.9198621771937625
SEED = 20260902
WEIGHT_DATA_SNAPSHOT_NAME = "WeightData_newdata167_accmix23_effmix19.root"
ACCEPTANCE_SNAPSHOT_NAME = "acceptance_sps0p8_dps0p2.txt"
EFFICIENCY_SNAPSHOT_NAME = "efficiency_sps0p8_dps0p2_avgacc_19x10.txt"
DPS_EXPECTED_FILES = 65
DPS_EXCLUDED_INDICES = set()
DATA_LABEL_SYMMETRIZE = False
DATA_LABEL_SEED = 50
DPS_LABEL_SYMMETRIZE = False
DPS_LABEL_SEED = 50
EXPECTED_DATA_RANDOMIZATION = None
CR_DEFINITIONS = [("phi70", 1.8, PHI70), ("phi110", 1.8, PHI110),
                  ("dy2p4", 2.4, PHI90)]
EXPECTED_PP_FITS = 51
EXTRA_SNAPSHOTS = []
RUNNER_RELATIVE = "Data_driven/run_newdata_nominal_chain.py"
GENERATION_LABEL = "newdata167_accmix23_effmix19"


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stage(message):
    print(f"[{datetime.now().astimezone().isoformat(timespec='seconds')}] {message}", flush=True)


def run(command, log_path=None):
    if log_path is None:
        subprocess.run(command, cwd=REPO, check=True)
        return
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w", encoding="utf-8") as log:
        subprocess.run(command, cwd=REPO, stdout=log, stderr=subprocess.STDOUT, check=True)


def root(expression, log_path):
    run(["root", "-l", "-b", "-q", expression], log_path)


def root_eval(expression, log_path):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w", encoding="utf-8") as log:
        completed = subprocess.run(
            ["root", "-l", "-b", "-q", "-e", expression], cwd=REPO,
            stdout=log, stderr=subprocess.STDOUT, check=False)
    need(completed.returncode in (0, 255),
         f"ROOT inline audit failed with exit {completed.returncode}")


def read_kv(path):
    result = {}
    with open(path, encoding="utf-8") as source:
        for raw in source:
            line = raw.rstrip("\n")
            if "=" in line:
                key, value = line.split("=", 1)
                result[key] = value
    return result


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def validate_pp(directory):
    summary = read_kv(directory / "summary.txt")
    cells = int(summary["cells"])
    need(int(summary["accepted_fits"]) == cells, f"PP fits not all accepted: {directory}")
    need(int(summary["failed_fits"]) == 0, f"Failed PP fit: {directory}")
    with open(directory / "fit_results.csv", newline="", encoding="utf-8") as source:
        rows = list(csv.DictReader(source))
    need(len(rows) == cells, f"PP CSV cell mismatch: {directory}")
    for row in rows:
        need(int(row["accepted"]) == 1, f"Rejected PP cell: {directory}")
        need(int(row["status"]) == 0, f"Bad PP status: {directory}")
        need(int(row["covQual"]) >= 2, f"Bad PP covQual: {directory}")
        need(float(row["edm"]) < 0.01, f"Bad PP EDM: {directory}")
    need(sum(int(row["comb_comb_fixed_zero"]) for row in rows) ==
         int(summary["comb_comb_fixed_zero_fits"]), f"Comb_Comb fallback mismatch: {directory}")
    need(sum(int(row["sig_comb_fixed_zero"]) for row in rows) ==
         int(summary["sig_comb_fixed_zero_fits"]), f"Sig_Comb fallback mismatch: {directory}")
    return summary, rows


def write_inventory(paths, output):
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", encoding="utf-8", newline="") as target:
        target.write("path\tsize_bytes\tmtime_epoch\tsha256\n")
        for path in paths:
            stat = path.stat()
            target.write(f"{path}\t{stat.st_size}\t{int(stat.st_mtime)}\t{sha256(path)}\n")


def render_pngs(plot_dir):
    pdfs = sorted(plot_dir.glob("*.pdf"))
    need(len(pdfs) == 9, f"Expected 9 PDF plots in {plot_dir}, got {len(pdfs)}")
    for pdf in pdfs:
        run(["pdftoppm", "-png", "-singlefile", "-r", "120", str(pdf), str(pdf.with_suffix(""))])
    need(len(list(plot_dir.glob("*.png"))) == 9, f"Expected 9 PNG plots in {plot_dir}")


def fmt(value):
    return f"{value:.12g}"


def main():
    if len(sys.argv) != 2:
        raise SystemExit(f"Usage: {sys.argv[0]} <unique-run-tag>")
    tag = sys.argv[1]
    need(tag and all(c.isalnum() or c in "._-" for c in tag), f"Unsafe tag: {tag}")
    need(subprocess.check_output(["root-config", "--version"], text=True).strip() == "6.40.02",
         "This chain requires ROOT 6.40.02")

    output = DD / "results" / tag
    need(not output.exists(), f"Refusing to overwrite existing tag: {output}")

    stage("Validate transferred model, WeightData, manifest and corrections")
    for path, expected in EXPECTED.items():
        need(path.is_file() and path.stat().st_size > 0, f"Missing input: {path}")
        need(sha256(path) == expected, f"Checksum mismatch: {path}")

    data_paths = []
    for raw in DATA_MANIFEST_SOURCE.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line and not line.startswith("#"):
            data_paths.append(Path(line))
    need(len(data_paths) == 167 and len(set(data_paths)) == 167,
         "Data manifest is not exactly 167 unique paths")
    need(all(path.is_file() and path.stat().st_size > 0 for path in data_paths),
         "At least one raw Data file is missing")

    all_dps_paths = sorted(DPS_BASE.glob("Ntuple_2016_DPS_*.root"),
                           key=lambda p: int(p.stem.rsplit("_", 1)[1]))
    all_dps_indices = {int(path.stem.rsplit("_", 1)[1]) for path in all_dps_paths}
    need(all_dps_indices == set(range(1, 66)),
         f"DPS source inventory is not exactly indices 1..65: {sorted(all_dps_indices)}")
    need(DPS_EXCLUDED_INDICES <= all_dps_indices,
         f"Unknown excluded DPS indices: {sorted(DPS_EXCLUDED_INDICES-all_dps_indices)}")
    dps_paths = [path for path in all_dps_paths
                 if int(path.stem.rsplit("_", 1)[1]) not in DPS_EXCLUDED_INDICES]
    need(len(dps_paths) == DPS_EXPECTED_FILES,
         f"Expected {DPS_EXPECTED_FILES} DPS files, got {len(dps_paths)}")
    need(all(path.stat().st_size > 0 for path in dps_paths), "Empty DPS input file")

    inp = output / "input"
    logs = output / "logs"
    inp.mkdir(parents=True)
    logs.mkdir()
    model = inp / "Model_4D_tot.root"
    fit_result = inp / "Fit_4D_tot_native_asymptotic_unseeded.root"
    weight_data = inp / WEIGHT_DATA_SNAPSHOT_NAME
    data_manifest = inp / "data_ntuple_manifest_chensh_newdata_v1.list"
    acceptance = inp / ACCEPTANCE_SNAPSHOT_NAME
    efficiency = inp / EFFICIENCY_SNAPSHOT_NAME
    snapshots = [
        (MODEL_SOURCE, model), (FIT_RESULT_SOURCE, fit_result),
        (WEIGHT_DATA_SOURCE, weight_data), (DATA_MANIFEST_SOURCE, data_manifest),
        (ACCEPTANCE_SOURCE, acceptance), (EFFICIENCY_SOURCE, efficiency),
        (SELECTION_SOURCE, inp / "rephrase.cpp"),
        (UPSTREAM / "total/Fit_4D_tot.cpp", inp / "Fit_4D_tot.cpp"),
    ] + EXTRA_SNAPSHOTS
    for source, target in snapshots:
        shutil.copy2(source, target)

    stage("Build raw Data and DPS input inventories")
    data_inventory = inp / "data_ntuple_inventory.tsv"
    dps_inventory = inp / "dps_ntuple_inventory.tsv"
    dps_manifest = inp / "dps_ntuple_manifest.list"
    write_inventory(data_paths, data_inventory)
    write_inventory(dps_paths, dps_inventory)
    dps_manifest.write_text("".join(f"{path}\n" for path in dps_paths), encoding="utf-8")

    stage("Audit snapshotted WeightData")
    audit_expr = (
        f'TFile f("{weight_data}"); auto *t=(TTree*)f.Get("data"); '
        'if(!t) gSystem->Exit(1); '
        'std::cout<<"weightdata_entries="<<t->GetEntries()<<std::endl; '
        'std::cout<<"fit_range_entries="<<t->GetEntries('
        '"Jpsi_mass1>=2.95&&Jpsi_mass1<=3.25&&Jpsi_mass2>=2.95&&Jpsi_mass2<=3.25&&'
        'Jpsi_ctau1>=-0.03&&Jpsi_ctau1<=0.16&&Jpsi_ctau2>=-0.03&&Jpsi_ctau2<=0.16")'
        '<<std::endl;'
    )
    root_eval(audit_expr, logs / "audit_weightdata.log")
    audit_text = (logs / "audit_weightdata.log").read_text(encoding="utf-8")
    need("weightdata_entries=2954" in audit_text, "WeightData entry mismatch")
    need("fit_range_entries=2900" in audit_text, "WeightData fit-range entry mismatch")

    data_candidates = inp / "jpsi12_candidates_newdata167.root"
    stage("Build single-J/psi candidates from the 167-file manifest")
    data_sym = "true" if DATA_LABEL_SYMMETRIZE else "false"
    root(f'Data_driven/build_jpsi12_candidates.cpp+("{data_candidates}","data","","{data_manifest}",{data_sym},{DATA_LABEL_SEED})',
         logs / "build_data_candidates.log")
    data_candidate_info = read_kv(logs / "build_data_candidates.log")
    need(int(data_candidate_info["input_files"]) == 167, "Candidate Data file count mismatch")
    need(int(data_candidate_info["selected_events"]) == 2954, "Candidate selected count mismatch")
    need(data_candidate_info["label_randomization"] ==
         ("enabled" if DATA_LABEL_SYMMETRIZE else "disabled"),
         "Candidate Data label-randomization mode mismatch")
    if EXPECTED_DATA_RANDOMIZATION is not None:
        for key, expected in EXPECTED_DATA_RANDOMIZATION.items():
            need(int(data_candidate_info[key]) == expected,
                 f"Candidate Data randomization mismatch for {key}")

    pp_nominal = output / "pp_data_2d"
    stage("Run nominal 12-cell PP fits")
    root(f'Data_driven/fit_pp_2d_adaptive.cpp+("{pp_nominal}","{weight_data}",true,false,true,1.8,{PHI90:.17g},"{model}")',
         logs / "fit_pp_nominal.log")
    nominal_pp_summary, nominal_pp_rows = validate_pp(pp_nominal)
    need(int(nominal_pp_summary["cells"]) == 12, "Nominal PP geometry is not 12 cells")

    cr_base = output / "cr"
    cr_pp = {}
    for label, dy_min, phi_max in CR_DEFINITIONS:
        stage(f"Run PP fits for CR variation {label}")
        directory = cr_base / f"{label}_pp"
        root(f'Data_driven/fit_pp_2d_adaptive.cpp+("{directory}","{weight_data}",true,false,true,{dy_min:.17g},{phi_max:.17g},"{model}")',
             logs / f"fit_pp_{label}.log")
        cr_pp[label] = (directory, *validate_pp(directory))

    stage("Fit single-J/psi sWeights and build Data event mixing")
    splot = output / "event_mixing/splot"
    splot.mkdir(parents=True)
    for slot in (1, 2):
        target = splot / f"jpsi{slot}_sweights.root"
        root(f'Data_driven/fit_single_jpsi_global_calibrated_splot.cpp+({slot},"{data_candidates}","{model}","{target}")',
             logs / f"fit_jpsi{slot}.log")
        summary = read_kv(splot / f"global_calibrated_splot_summary_jpsi{slot}.txt")
        need(summary["fit_mode"] == "single_global_calibrated_two_stage", "Wrong single-J fit mode")
        need(int(summary["calibration_status"]) == 0 and int(summary["calibration_covQual"]) >= 2
             and float(summary["calibration_edm"]) < 0.01, "Single-J calibration fit failed QA")
        need(int(summary["yield_fit_status"]) == 0 and int(summary["yield_fit_covQual"]) >= 2
             and float(summary["yield_fit_edm"]) < 0.01, "Single-J yield fit failed QA")

    event_mix = output / "event_mixing/mixed_dps_allpairs.root"
    root(f'Data_driven/build_crossslot_splot_mixed.cpp+(0,0,"{splot}/jpsi1_sweights.root","{splot}/jpsi2_sweights.root","{event_mix}","{acceptance}","{efficiency}")',
         logs / "build_data_event_mixing.log")
    event_mix_info = read_kv(logs / "build_data_event_mixing.log")
    need(event_mix_info["mixing_mode"] == "all_cross_event_pairs", "Wrong Data mixing mode")
    need(int(event_mix_info["correction_rejected"]) == 0, "Data mixing correction rejection")
    event_template = output / "event_mixing/templates_2d"
    root(f'Data_driven/build_2d_templates_adaptive.cpp+("{event_template}","{pp_nominal}/fit_results.csv","{event_mix}",1.8,{PHI90:.17g},0.0,true,"data_single_jpsi_event_mixing_crosscheck")',
         logs / "build_event_mixing_template.log")

    stage("Rebuild WeightDPS and direct DPS template with new corrections")
    dps_reference = output / "dps_reference"
    dps_reference.mkdir()
    weight_dps = dps_reference / "WeightDPS_current.root"
    full_dps_input = dps_reference / "nominal_dps_template_input.root"
    dps_sym = "true" if DPS_LABEL_SYMMETRIZE else "false"
    root(f'Data_driven/build_dps_reference.cpp+("{weight_dps}","{dps_manifest}","{acceptance}","{efficiency}",{dps_sym},{DPS_LABEL_SEED},{DPS_EXPECTED_FILES})',
         logs / "build_weightdps.log")
    need(weight_dps.is_file() and weight_dps.stat().st_size > 0, "WeightDPS not built")
    root(f'Data_driven/prepare_dps_template_input.cpp+("{weight_dps}","{full_dps_input}")',
         logs / "prepare_dps_template.log")
    prepared = read_kv(logs / "prepare_dps_template.log")
    full_entries = int(prepared["accepted_entries"])
    need(full_entries > 1, "DPS selected entries are not positive")
    need(int(prepared["rejected_entries"]) == 0, "DPS template adapter rejected entries")
    half_entries = full_entries // 2
    half_policy = ("exact_half_without_replacement" if full_entries % 2 == 0
                   else "floor_half_without_replacement")

    full_template = output / "full_dpsmc_reference/templates_2d"
    stage("Build full-DPS-MC reference")
    root(f'Data_driven/build_2d_templates_adaptive.cpp+("{full_template}","{pp_nominal}/fit_results.csv","{full_dps_input}",1.8,{PHI90:.17g},0.0,false,"dps_reconstruction_mc")',
         logs / "build_full_dpsmc_template.log")

    half_base = output / f"nominal_fixed_half_seed{SEED}"
    half_input = half_base / "dps_mc_half.root"
    half_sampling = half_base / "sampling_summary.txt"
    half_selection = half_base / "selected_entries.csv"
    stage("Draw exact half of DPS MC and build nominal template")
    root(f'Data_driven/subsample_dps_template_fixed.cpp+("{full_dps_input}","{half_input}","{half_sampling}","{half_selection}",{SEED})',
         logs / "build_fixed_half.log")
    sampling = read_kv(half_sampling)
    need(int(sampling["input_entries"]) == full_entries, "Half-sample input count mismatch")
    need(int(sampling["selected_entries"]) == half_entries, "Half-sample count mismatch")
    need(int(sampling["seed"]) == SEED and int(sampling["priority_collisions"]) == 0,
         "Half-sample seed/collision mismatch")
    half_template = half_base / "templates_2d"
    root(f'Data_driven/build_2d_templates_adaptive.cpp+("{half_template}","{pp_nominal}/fit_results.csv","{half_input}",1.8,{PHI90:.17g},0.0,false,"dps_reconstruction_mc_fixed_seed_half")',
         logs / "build_half_nominal_template.log")
    half_summary = read_kv(half_template / "summary.txt")
    half_fdps = float(half_summary["f_dps"])
    half_fsps = float(half_summary["f_sps"])

    variations = []
    for label, dy_min, phi_max in CR_DEFINITIONS:
        stage(f"Build same-half template for CR variation {label}")
        pp_dir = cr_pp[label][0]
        target = half_base / "cr" / label
        root(f'Data_driven/build_2d_templates_adaptive.cpp+("{target}","{pp_dir}/fit_results.csv","{half_input}",{dy_min:.17g},{phi_max:.17g},0.0,false,"dps_reconstruction_mc_fixed_seed_half")',
             logs / f"build_half_{label}.log")
        info = read_kv(target / "summary.txt")
        value = float(info["f_dps"])
        variations.append((label, dy_min, phi_max, int(info["cells"]), value,
                           value - half_fdps, int(info["negative_sps_bins"]),
                           int(cr_pp[label][1]["accepted_fits"])))
    with open(output / "cr_variations.csv", "w", newline="", encoding="utf-8") as target:
        writer = csv.writer(target)
        writer.writerow(["variation", "dy_min", "phi_max", "cells", "f_dps",
                         "signed_shift", "negative_sps_bins", "pp_accepted"])
        writer.writerow(["nominal", 1.8, PHI90, 12, fmt(half_fdps), 0,
                         half_summary["negative_sps_bins"], 12])
        for row in variations:
            writer.writerow([row[0], fmt(row[1]), fmt(row[2]), row[3], fmt(row[4]),
                             fmt(row[5]), row[6], row[7]])

    stage("Build DPS-MC self-mixing usage systematic and 1D QA")
    usage = output / "dps_mc_usage_systematic"
    usage.mkdir()
    usage_candidates = usage / "dps_mc_candidates.root"
    root(f'Data_driven/build_jpsi12_candidates.cpp+("{usage_candidates}","dps_mc","","{dps_manifest}",false,{DPS_LABEL_SEED})',
         logs / "build_dps_candidates.log")
    dps_candidate_info = read_kv(logs / "build_dps_candidates.log")
    need(int(dps_candidate_info["input_files"]) == DPS_EXPECTED_FILES,
         "DPS candidate file count mismatch")
    need(dps_candidate_info["label_randomization"] == "disabled",
         "DPS candidates must not receive extra label randomization")
    need(int(dps_candidate_info["selected_events"]) == full_entries,
         "DPS candidates and direct selected entries differ")
    usage_mix = usage / "dps_mc_half_crossslot_mixed.root"
    usage_mix_summary = usage / "mixing_summary.txt"
    root(f'Data_driven/mix_dps_mc_half_crossslot.cpp+("{usage_candidates}","{half_selection}","{half_input}","{usage_mix}","{usage_mix_summary}","{acceptance}","{efficiency}")',
         logs / "build_dps_self_mixing.log")
    mixing = read_kv(usage_mix_summary)
    need(mixing["status"] == "complete", "DPS self-mixing did not complete")
    need(int(mixing["selected_pair_events"]) == half_entries, "DPS self-mixing source count mismatch")
    need(int(mixing["direct_half_validation_mismatches"]) == 0, "Direct/self-mix source mismatch")
    need(int(mixing["correction_rejected"]) == 0, "DPS self-mixing correction rejection")
    usage_template = usage / "templates_2d"
    root(f'Data_driven/build_2d_templates_adaptive.cpp+("{usage_template}","{pp_nominal}/fit_results.csv","{usage_mix}",1.8,{PHI90:.17g},0.0,false,"dps_mc_fixed_half_crossslot_event_mixing")',
         logs / "build_usage_template.log")
    usage_shape = usage / "shape_comparison_1d"
    root(f'Data_driven/compare_mixed_dps_to_mc_1d.cpp+("{usage_shape}","{usage_mix}","{half_input}","dps_mc","DPS MC event mixing","DPS MC direct","DPS MC mixing closure","shared_dps_mc_source_events_cross_covariance_not_propagated")',
         logs / "compare_usage_1d.log")
    need(len(list(csv.DictReader(open(usage_shape / "shape_metrics.csv", encoding="utf-8")))) == 9,
         "Usage shape metric count mismatch")
    render_pngs(usage_shape / "plots")

    stage("Build Data event-mixing versus full DPS-MC 1D cross-check")
    event_shape = output / "event_mixing/shape_comparison_1d"
    root(f'Data_driven/compare_mixed_dps_to_mc_1d.cpp+("{event_shape}","{event_mix}","{full_dps_input}","original_data","single-J/psi event mixing","DPS MC","DPS shape validation","independent_samples")',
         logs / "compare_event_mixing_1d.log")
    need(len(list(csv.DictReader(open(event_shape / "shape_metrics.csv", encoding="utf-8")))) == 9,
         "Event-mixing shape metric count mismatch")
    render_pngs(event_shape / "plots")

    full_summary = read_kv(full_template / "summary.txt")
    full_fdps = float(full_summary["f_dps"])
    full_fsps = float(full_summary["f_sps"])
    stat_signed = half_fdps - full_fdps
    stat_uncertainty = abs(stat_signed)
    cr_up = max(0.0, *(row[5] for row in variations))
    cr_down = max(0.0, *(-row[5] for row in variations))
    usage_summary = read_kv(usage_template / "summary.txt")
    usage_fdps = float(usage_summary["f_dps"])
    usage_signed = usage_fdps - half_fdps
    usage_uncertainty = abs(usage_signed)
    stat_percent = 100.0 * stat_uncertainty / half_fdps
    cr_up_percent = 100.0 * cr_up / half_fdps
    cr_down_percent = 100.0 * cr_down / half_fdps
    usage_percent = 100.0 * usage_uncertainty / half_fdps
    total_systematic_up = math.hypot(cr_up, usage_uncertainty)
    total_systematic_down = math.hypot(cr_down, usage_uncertainty)
    total_systematic_up_percent = 100.0 * total_systematic_up / half_fdps
    total_systematic_down_percent = 100.0 * total_systematic_down / half_fdps
    event_summary = read_kv(event_template / "summary.txt")
    event_fdps = float(event_summary["f_dps"])
    event_signed = event_fdps - half_fdps
    event_relative = 100.0 * event_signed / half_fdps

    pp_groups = [(nominal_pp_summary, nominal_pp_rows)]
    pp_groups.extend((item[1], item[2]) for item in cr_pp.values())
    total_pp = sum(int(summary["cells"]) for summary, _ in pp_groups)
    comb_fixed_zero = sum(int(summary["comb_comb_fixed_zero_fits"]) for summary, _ in pp_groups)
    sig_fixed_zero = sum(int(summary["sig_comb_fixed_zero_fits"]) for summary, _ in pp_groups)
    min_covqual = min(int(row["covQual"]) for _, rows in pp_groups for row in rows)
    max_edm = max(float(row["edm"]) for _, rows in pp_groups for row in rows)
    need(total_pp == EXPECTED_PP_FITS,
         f"Expected {EXPECTED_PP_FITS} PP fits, got {total_pp}")

    summary_lines = [
        "status=candidate_complete_pending_user_confirmation",
        "artifact=newdata_fixed_half_dpsmc_nominal_fdps",
        f"upstream_data_tag={UPSTREAM_TAG}",
        "upstream_partial_generation_gate=user_message_accepted_as_data_driven_start_gate",
        "method_nominal_status=user_confirmed_method_new_values_pending_confirmation",
        f"method=fixed_seed_{half_policy}_DPS_reconstruction_MC_template_CR_normalized",
        f"f_dps={fmt(half_fdps)}", f"f_sps={fmt(half_fsps)}",
        "stat_definition=absolute_shift_fixed_half_vs_full_DPS_MC",
        f"stat_uncertainty={fmt(stat_uncertainty)}",
        f"stat_uncertainty_percent={fmt(stat_percent)}",
        f"stat_signed_half_minus_full={fmt(stat_signed)}",
        f"stat_full_reference_f_dps={fmt(full_fdps)}",
        f"stat_full_reference_f_sps={fmt(full_fsps)}",
        f"stat_seed={SEED}", "stat_engine=std_mt19937_64",
        f"stat_sampling={half_policy}",
        f"stat_input_entries={full_entries}", f"stat_selected_entries={half_entries}",
        f"stat_selected_fraction={fmt(half_entries / full_entries)}",
        "cr_systematic_definition=zero_leakage_" +
        "_".join(row[0] for row in CR_DEFINITIONS) +
        "_envelope_same_fixed_half",
        f"cr_systematic_up={fmt(cr_up)}", f"cr_systematic_down={fmt(cr_down)}",
        f"cr_systematic_symmetric_max={fmt(max(cr_up, cr_down))}",
        f"cr_systematic_up_percent={fmt(cr_up_percent)}",
        f"cr_systematic_down_percent={fmt(cr_down_percent)}",
    ]
    for row in variations:
        summary_lines.extend([
            f"cr_{row[0]}_f_dps={fmt(row[4])}",
            f"cr_{row[0]}_signed_shift={fmt(row[5])}",
            f"cr_{row[0]}_negative_sps_bins={row[6]}",
        ])
    summary_lines.extend([
        "dps_mc_usage_systematic_definition=absolute_shift_self_mixed_half_vs_direct_half",
        f"dps_mc_self_mixed_f_dps={fmt(usage_fdps)}",
        f"dps_mc_usage_signed_shift={fmt(usage_signed)}",
        f"dps_mc_usage_systematic={fmt(usage_uncertainty)}",
        f"dps_mc_usage_systematic_percent={fmt(usage_percent)}",
        "dps_mc_event_mixing_systematic_alias=dps_mc_usage_systematic",
        "total_systematic_definition=quadrature_CR_and_DPS_MC_event_mixing_usage",
        f"total_systematic_up={fmt(total_systematic_up)}",
        f"total_systematic_down={fmt(total_systematic_down)}",
        f"total_systematic_up_percent={fmt(total_systematic_up_percent)}",
        f"total_systematic_down_percent={fmt(total_systematic_down_percent)}",
        "event_mixing_role=cross_check_not_nominal",
        f"event_mixing_f_dps={fmt(event_fdps)}",
        f"event_mixing_signed_shift_from_nominal={fmt(event_signed)}",
        f"event_mixing_relative_shift_percent={fmt(event_relative)}",
        "sps_leakage_fraction=0_for_central_and_all_CR_variations",
        "weightdata_selected_entries=2954", "weightdata_fit_range_entries=2900",
        f"data_label_randomization={'enabled' if DATA_LABEL_SYMMETRIZE else 'disabled'}",
        f"data_label_randomization_seed={DATA_LABEL_SEED}",
        f"dps_input_label_randomization={'enabled' if DPS_LABEL_SYMMETRIZE else 'disabled'}",
        "dps_raw_already_pair_symmetrized=true",
        f"dps_ntuple_files={DPS_EXPECTED_FILES}",
        "dps_excluded_indices=" +
        (",".join(str(index) for index in sorted(DPS_EXCLUDED_INDICES))
         if DPS_EXCLUDED_INDICES else "none"),
        f"dps_reco_selected_entries={full_entries}", f"pp_fits_total={total_pp}",
        "pp_fit_failed=0", f"pp_fit_min_covQual={min_covqual}",
        f"pp_fit_max_edm={fmt(max_edm)}", f"pp_fit_comb_comb_fixed_zero={comb_fixed_zero}",
        f"pp_fit_sig_comb_fixed_zero={sig_fixed_zero}",
        "pp_fit_error_convention=ROOT_native_AsymptoticError_true",
        "pp_comb_comb_boundary_strategy=fix_zero_and_refit_on_root_range_exception",
        "pp_sig_comb_boundary_strategy=fix_shared_zero_for_Sig_Comb_and_Comb_Sig_and_refit_on_root_range_exception",
        "usage_shape_plots=9", "event_mixing_shape_plots=9",
        "shape_metrics=max_CDF_Jensen_Shannon_and_covariance_chi2",
        "sigma_dps_status=pending_new_total_cross_section_not_calculated",
        "sigma_eff_status=pending_new_total_cross_section_not_calculated",
        "bootstrap_status=not_run_by_user_gate",
        "category_systematic_status=not_run_by_user_gate",
        "pure_dps_closure_status=not_run_by_user_gate",
    ])
    (output / "summary.txt").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    root_version = subprocess.check_output(["root-config", "--version"], text=True).strip()
    git_head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip()
    metadata_lines = [
        "status=candidate_complete_pending_user_confirmation", f"output_tag={tag}",
        f"generated_at={datetime.now().astimezone().isoformat(timespec='seconds')}",
        f"root_version={root_version}", f"git_head={git_head}",
        f"upstream_data_tag={UPSTREAM_TAG}",
        f"model_snapshot={model}", f"model_sha256={sha256(model)}",
        f"fit_result_snapshot={fit_result}", f"fit_result_sha256={sha256(fit_result)}",
        f"weightdata_snapshot={weight_data}", f"weightdata_sha256={sha256(weight_data)}",
        f"data_manifest_snapshot={data_manifest}", f"data_manifest_sha256={sha256(data_manifest)}",
        "data_ntuple_files=167", f"data_ntuple_inventory={data_inventory}",
        f"data_ntuple_inventory_sha256={sha256(data_inventory)}",
        f"dps_ntuple_files={DPS_EXPECTED_FILES}",
        "dps_ntuple_excluded_indices=" +
        (",".join(str(index) for index in sorted(DPS_EXCLUDED_INDICES))
         if DPS_EXCLUDED_INDICES else "none"),
        f"dps_ntuple_manifest={dps_manifest}",
        f"dps_ntuple_inventory={dps_inventory}",
        f"dps_ntuple_inventory_sha256={sha256(dps_inventory)}",
        f"acceptance_snapshot={acceptance}", f"acceptance_sha256={sha256(acceptance)}",
        f"efficiency_snapshot={efficiency}", f"efficiency_sha256={sha256(efficiency)}",
        f"selection_source_snapshot={inp / 'rephrase.cpp'}",
        f"selection_source_sha256={sha256(inp / 'rephrase.cpp')}",
        f"nominal_definition=fixed_seed_{half_policy}_DPS_reconstruction_MC_template_CR_normalized",
        f"seed_policy=std_mt19937_64_seed_{SEED}_{half_policy}",
        "cr_variations=" + "_".join(row[0] for row in CR_DEFINITIONS) +
        "_zero_SPS_leakage",
        f"data_label_randomization={'enabled' if DATA_LABEL_SYMMETRIZE else 'disabled'}",
        f"data_label_randomization_seed={DATA_LABEL_SEED}",
        f"dps_input_label_randomization={'enabled' if DPS_LABEL_SYMMETRIZE else 'disabled'}",
        f"dps_input_label_randomization_seed={DPS_LABEL_SEED}",
        "dps_raw_already_pair_symmetrized=true",
        "data_event_mixing=deterministic_all_cross_Jpsi1_A_x_Jpsi2_B_A_not_equal_B",
        "dps_self_mixing=deterministic_all_cross_Jpsi1_A_x_Jpsi2_B_A_not_equal_B",
        "all_pair_weights=same_snapshotted_acceptance_and_efficiency",
        "fit_policy=ROOT_6.40.02_native_AsymptoticError_true_nCombComb_and_shared_nSigComb_boundary_fix0",
        "systematic_presentation=CR_and_DPS_MC_event_mixing_usage_percentages_plus_quadrature_total",
        "final_numeric_status=pending_user_confirmation",
        "sigma_postprocess=not_run_waiting_new_total_cross_section",
    ]
    code_paths = [
        "Data_driven/build_jpsi12_candidates.cpp", "Data_driven/build_dps_reference.cpp",
        "Data_driven/prepare_dps_template_input.cpp", "Data_driven/fit_pp_2d_adaptive.cpp",
        "Data_driven/fit_single_jpsi_category_splot.cpp",
        "Data_driven/fit_single_jpsi_global_calibrated_splot.cpp",
        "Data_driven/build_mixed_dps.cpp", "Data_driven/build_crossslot_splot_mixed.cpp",
        "Data_driven/build_2d_templates_adaptive.cpp",
        "Data_driven/subsample_dps_template_fixed.cpp",
        "Data_driven/mix_dps_mc_half_crossslot.cpp",
        "Data_driven/compare_mixed_dps_to_mc_1d.cpp",
        "Data_driven/run_newdata_nominal_chain.py",
        RUNNER_RELATIVE,
    ]
    for relative in code_paths:
        metadata_lines.append(f"{sha256(REPO / relative)}  {relative}")
    (output / "run.metadata.txt").write_text("\n".join(metadata_lines) + "\n", encoding="utf-8")

    stage("Finalize and verify artifact checksums")
    checksum_path = output / "artifact_checksums.txt"
    artifacts = sorted(path for path in output.rglob("*")
                       if path.is_file() and path.name not in
                       {"artifact_checksums.txt", "complete.marker"})
    with open(checksum_path, "w", encoding="utf-8") as target:
        for path in artifacts:
            target.write(f"{sha256(path)}  {path}\n")
    for line in checksum_path.read_text(encoding="utf-8").splitlines():
        expected, name = line.split("  ", 1)
        need(sha256(Path(name)) == expected, f"Artifact checksum failed: {name}")
    (output / "complete.marker").write_text(
        "candidate_complete_pending_user_confirmation\n", encoding="utf-8")

    stage(
        "NEW_DATA_FDPS_CANDIDATE_COMPLETE "
        f"tag={tag} f_DPS={fmt(half_fdps)} stat={fmt(stat_uncertainty)} "
        f"cr_up={fmt(cr_up)} cr_down={fmt(cr_down)} "
        f"usage={fmt(usage_uncertainty)}"
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"NEW_DATA_FDPS_CHAIN_FAILED {error}", file=sys.stderr, flush=True)
        raise
