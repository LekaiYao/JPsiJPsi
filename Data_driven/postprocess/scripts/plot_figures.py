#!/usr/bin/env python3
"""Create H016 result-summary figures from results.json."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


PACKAGE = Path(__file__).resolve().parents[1]
if len(sys.argv) > 3:
    raise SystemExit("usage: plot_figures.py [INPUT_JSON] [OUTPUT_DIR]")
INPUT_PATH = Path(sys.argv[1]).resolve() if len(sys.argv) >= 2 else PACKAGE / "inputs.json"
OUTPUT = Path(sys.argv[2]).resolve() if len(sys.argv) >= 3 else PACKAGE / "output"
FIGURES = OUTPUT / "figures"
RESULTS = json.loads((OUTPUT / "results.json").read_text(encoding="utf-8"))
INPUTS = json.loads(INPUT_PATH.read_text(encoding="utf-8"))
FROZEN_H016_REFERENCE = (
    RESULTS.get("status") == "producer_ready_pending_consumer_review"
)

BLUE = "#3b6ea8"
ORANGE = "#d97732"
GREEN = "#2a8f6a"
PURPLE = "#7655a6"
DARK = "#263238"
LIGHT = "#e9eef5"


def setup_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 16,
            "axes.titlesize": 23,
            "axes.labelsize": 19,
            "xtick.labelsize": 15,
            "ytick.labelsize": 15,
            "legend.fontsize": 14,
            "axes.linewidth": 1.3,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save(fig: plt.Figure, stem: str) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES / f"{stem}.pdf", facecolor="white")
    fig.savefig(FIGURES / f"{stem}.png", dpi=180, facecolor="white")
    plt.close(fig)


def result_line(ax: plt.Axes, text: str, color: str = DARK, y: float = 0.72) -> None:
    ax.text(0.5, y, text, transform=ax.transAxes, ha="center", va="center",
            fontsize=22, fontweight="bold", color=color)


def fdps_line(r: dict) -> str:
    if FROZEN_H016_REFERENCE:
        return (
            rf"$f_{{DPS}}={r['central']:.3f}\pm {r['stat']:.3f}\,(\mathrm{{stat.}})"
            rf"\,{{}}^{{+{r['method_syst_up']:.3f}}}_{{-{r['method_syst_down']:.3f}}}"
            rf"\,(\mathrm{{method\ syst.}})$"
        )
    return (
        f"f_DPS = {r['central']:.3f} +/- {r['stat']:.3f} (stat.) "
        f"+{r['method_syst_up']:.3f}/-{r['method_syst_down']:.3f} (method syst.)"
    )


def cross_section_line(r: dict, symbol: str, unit: str) -> str:
    if FROZEN_H016_REFERENCE:
        return (
            rf"${symbol}={r['central']:.2f}\pm {r['stat_rho_scan']['0.0']:.2f}"
            rf"\,(\mathrm{{stat.}})\,{{}}^{{+{r['syst_up']:.2f}}}_{{-{r['syst_down']:.2f}}}"
            rf"\,(\mathrm{{syst.}})\ \mathrm{{{unit}}}$"
        )
    plain_symbol = "sigma_DPS" if "DPS" in symbol else "sigma_eff"
    return (
        f"{plain_symbol} = {r['central']:.2f} +/- "
        f"{r['stat_rho_scan']['0.0']:.2f} (stat.) "
        f"+{r['syst_up']:.2f}/-{r['syst_down']:.2f} (syst.) {unit}"
    )


def plot_fdps() -> None:
    r = RESULTS["f_DPS"]
    groups = INPUTS["data_driven_fraction"]["method_groups"]
    fig, (left, right) = plt.subplots(
        1, 2, figsize=(13.333, 7.5), gridspec_kw={"width_ratios": [1.03, 0.97]},
        constrained_layout=True,
    )
    fig.suptitle(
        r"Current data-driven $f_{DPS}$ result"
        if FROZEN_H016_REFERENCE
        else r"Candidate data-driven $f_{DPS}$ result",
        fontsize=26,
        fontweight="bold",
    )

    central = r["central"]
    left.axvline(central, color=DARK, linestyle="--", linewidth=1.5, alpha=0.7)
    left.errorbar(
        central, 0.62, xerr=r["stat"], fmt="o", markersize=11, color=BLUE,
        linewidth=3.0, capsize=8, label="Statistical (fixed-500 bootstrap)",
    )
    left.errorbar(
        central, 0.38,
        xerr=np.array([[r["method_syst_down"]], [r["method_syst_up"]]]),
        fmt="s", markersize=9, color=ORANGE, linewidth=3.0, capsize=8,
        label="Method systematic",
    )
    left.set_ylim(0.12, 0.88)
    span = 1.35 * max(r["stat"], r["method_syst_up"], r["method_syst_down"])
    left.set_xlim(max(0.0, central - span), min(1.0, central + span))
    left.set_yticks([0.38, 0.62])
    left.set_yticklabels(["Method syst.", "Stat."])
    left.set_xlabel(r"$f_{DPS}$")
    left.grid(axis="x", alpha=0.25)
    left.legend(loc="lower center", bbox_to_anchor=(0.5, -0.24), frameon=False)
    left.text(
        0.5, 0.96,
        fdps_line(r),
        transform=left.transAxes, ha="center", va="top", fontsize=18,
    )

    labels = ["Single-J category", "CR definition / purity", "Pure-DPS closure"]
    keys = ["single_j_category", "cr_definition_and_purity", "pure_dps_full_workflow_closure"]
    y = np.arange(3)[::-1]
    up = np.array([groups[key]["up"] for key in keys])
    down = np.array([groups[key]["down"] for key in keys])
    right.barh(y, up, color=GREEN, height=0.34, label="Up")
    right.barh(y, -down, color=PURPLE, height=0.34, label="Down")
    for yi, u, d in zip(y, up, down):
        right.text(u + 0.0002, yi, f"+{u:.4f}", va="center", fontsize=13)
        right.text(-d - 0.0002, yi, f"-{d:.4f}", ha="right", va="center", fontsize=13)
    right.axvline(0.0, color=DARK, linewidth=1.2)
    right.set_yticks(y)
    right.set_yticklabels(labels)
    group_span = 1.25 * max(max(up), max(down), 1e-6)
    right.set_xlim(-group_span, group_span)
    right.set_xlabel(r"Absolute effect on $f_{DPS}$")
    right.set_title(
        f"Three method-systematic groups\n"
        f"Separate quadrature: +{r['method_syst_up']:.4f}/-{r['method_syst_down']:.4f}",
        fontsize=15, pad=8,
    )
    right.grid(axis="x", alpha=0.22)
    right.legend(loc="lower center", bbox_to_anchor=(0.5, -0.22), ncol=2, frameon=False)
    save(fig, "fdps_current_result")


def plot_sigma_dps() -> None:
    r = RESULTS["sigma_DPS_mJJ_gt_7p5_pb"]
    no_mass = RESULTS["sigma_DPS_no_mJJ_cut_pb"]
    c = r["components"]
    fig, (left, right) = plt.subplots(
        1, 2, figsize=(13.333, 7.5), gridspec_kw={"width_ratios": [1.05, 0.95]},
        constrained_layout=True,
    )
    fig.suptitle(r"Data-driven DPS cross section", fontsize=26, fontweight="bold")

    left.axis("off")
    left.text(0.5, 0.90, r"$\sigma_{DPS}=f_{DPS}\,\sigma_{pair}$", ha="center", fontsize=25)
    pair_pb = float(INPUTS["pair_cross_section"]["central_pb"])
    left.text(0.5, 0.78, f"= {RESULTS['f_DPS']['central']:.6f} x {pair_pb:.1f} pb" if FROZEN_H016_REFERENCE else f"= {RESULTS['f_DPS']['central']:.6f} x {pair_pb:.3f} pb",
              ha="center", fontsize=20)
    result_line(
        left,
        cross_section_line(r, r"\sigma_{DPS}", "pb"),
        BLUE,
        0.61,
    )
    left.text(0.5, 0.50, r"Fiducial pair region: $m(J/\psi J/\psi)>7.5$ GeV",
              ha="center", fontsize=17)
    left.annotate("", xy=(0.50, 0.31), xytext=(0.50, 0.44), xycoords="axes fraction",
                  arrowprops={"arrowstyle": "-|>", "lw": 2.0, "color": DARK})
    left.text(0.5, 0.25, f"For sigma_eff: divide by epsilon_m={RESULTS['mass_extrapolation']['epsilon_m']:.6f}",
              ha="center", fontsize=17)
    left.text(0.5, 0.14, rf"$\sigma_{{DPS}}^{{\mathrm{{no\ mass\ cut}}}}={no_mass['central']:.3f}$ pb",
              ha="center", fontsize=20, color=GREEN, fontweight="bold")
    left.text(0.5, 0.03, "Preliminary result" if FROZEN_H016_REFERENCE else "Candidate result; pending user confirmation",
              ha="center", fontsize=14, color="#5b6470")

    f_group_effects = c["f_DPS_method_groups"]
    pair_effects = c["AN_pair_sources"]
    f_method_up = np.sqrt(sum(effect["up"] ** 2 for effect in f_group_effects.values()))
    f_method_down = np.sqrt(sum(effect["down"] ** 2 for effect in f_group_effects.values()))
    labels = [
        r"", "Pair total stat.", r" syst.",
        "Pair BR", "CMS luminosity", "Acceptance / efficiency", "Lifetime variable", "Fitter stability",
    ]
    effects = [
        c["stat_from_f_DPS"], c["stat_from_pair_total"], max(f_method_up, f_method_down),
        pair_effects["branching_fraction"], pair_effects["cms_luminosity"],
        pair_effects["acceptance_efficiency"], pair_effects["lifetime_variable"],
        pair_effects["fitter_stability"],
    ]
    colors = [BLUE, BLUE] + [ORANGE] * 6
    y = np.arange(len(labels))[::-1]
    right.barh(y, effects, color=colors, height=0.52)
    for yi, value in zip(y, effects):
        right.text(value + 0.04, yi, f"{value:.3f} pb", va="center", fontsize=14)
    right.set_yticks(y)
    right.set_yticklabels(labels)
    right.set_ylim(-0.5, len(labels) + 0.2)
    right.tick_params(axis="y", labelsize=12)
    right.set_xlim(0, 2.85)
    right.set_xlabel("Absolute propagated component")
    right.set_title(
        f"Source-level components\nTotal stat.: {r['stat_rho_scan']['0.0']:.3f} pb\n"
        f"Total syst.: +{r['syst_up']:.3f}/-{r['syst_down']:.3f} pb",
        fontsize=14, pad=6,
    )
    right.grid(axis="x", alpha=0.24)
    save(fig, "sigma_dps_preliminary_breakdown")


def plot_sigma_eff() -> None:
    r = RESULTS["sigma_eff_mb"]
    c = r["components"]
    mass = RESULTS["mass_extrapolation"]
    fig, (left, right) = plt.subplots(
        1, 2, figsize=(13.333, 7.5), gridspec_kw={"width_ratios": [1.08, 0.92]},
        constrained_layout=True,
    )
    fig.suptitle(r"Effective DPS cross section", fontsize=26, fontweight="bold")

    left.axis("off")
    left.text(0.5, 0.91, r"$\sigma_{eff}=\frac{1}{2}\frac{\sigma^2(J/\psi)}{\sigma_{DPS}^{\mathrm{no\ mass\ cut}}}$",
              ha="center", fontsize=25)
    left.text(0.5, 0.76, r"$X=B\,\sigma(J/\psi)=10.63$ nb", ha="center", fontsize=18)
    left.text(0.5, 0.66, rf"$\epsilon_m={mass['epsilon_m']:.6f}$  from DPS GEN-only event mixing",
              ha="center", fontsize=17)
    left.text(0.5, 0.56, r"$B$ cancels between $X^2$ and the pair-visible denominator",
              ha="center", fontsize=16, color=GREEN)
    result_line(
        left,
        cross_section_line(r, r"\sigma_{eff}", "mb"),
        BLUE,
        0.38,
    )
    left.text(0.5, 0.22, "Prompt J/psi: 10 < pT < 40 GeV, |y| < 2.0",
              ha="center", fontsize=16)
    left.text(0.5, 0.10, "Preliminary result" if FROZEN_H016_REFERENCE else "Candidate result; pending user confirmation",
              ha="center", fontsize=14, color="#5b6470")

    f_endpoints = c["f_DPS_method_endpoint_sources"]
    pair_endpoints = c["CMS_pair_endpoint_sources"]
    labels = [
        r"", "Pair total stat.", "ATLAS single stat.",
        "Single-J category", "CR definition / purity", "Pure-DPS closure",
        "CMS luminosity", "Acceptance / efficiency", "Lifetime variable", "Fitter stability",
        "ATLAS bin syst.", "ATLAS luminosity",
    ]
    effects = [
        c["stat_from_f_DPS"], c["stat_from_pair_total"], c["stat_from_ATLAS_single"],
        max(f_endpoints["single_j_category"].values()),
        max(f_endpoints["cr_definition_and_purity"].values()),
        max(f_endpoints["pure_dps_full_workflow_closure"].values()),
        max(pair_endpoints["cms_luminosity"].values()),
        max(pair_endpoints["acceptance_efficiency"].values()),
        max(pair_endpoints["lifetime_variable"].values()),
        max(pair_endpoints["fitter_stability"].values()),
        max(c["ATLAS_single_bin_endpoint"].values()),
        max(c["ATLAS_luminosity_endpoint"].values()),
    ]
    colors = [BLUE, BLUE, BLUE] + [ORANGE] * 9
    y = np.arange(len(labels))[::-1]
    right.barh(y, effects, color=colors, height=0.50)
    for yi, value in zip(y, effects):
        right.text(value + 0.012, yi, f"{value:.3f}", va="center", fontsize=10.5)
    right.set_yticks(y)
    right.set_yticklabels(labels)
    right.set_xlim(0, 1.08)
    right.set_xlabel("Absolute effect [mb]")
    right.set_ylim(-0.5, len(labels) + 0.8)
    right.tick_params(axis="y", labelsize=11)
    right.set_title(
        f"Exact source-by-source endpoints\nBR cancels; totals: "
        f"{r['stat_rho_scan']['0.0']:.3f} stat., +{r['syst_up']:.3f}/-{r['syst_down']:.3f} syst.",
        fontsize=13, pad=8,
    )
    right.grid(axis="x", alpha=0.24)
    save(fig, "sigma_eff_preliminary_breakdown")


def plot_summary() -> None:
    f = RESULTS["f_DPS"]
    d = RESULTS["sigma_DPS_mJJ_gt_7p5_pb"]
    e = RESULTS["sigma_eff_mb"]
    fig, ax = plt.subplots(figsize=(13.333, 7.5), constrained_layout=True)
    ax.axis("off")
    ax.set_title(
        "Current data-driven results" if FROZEN_H016_REFERENCE else "Candidate data-driven results",
        fontsize=28,
        fontweight="bold",
        pad=18,
    )

    y = [0.74, 0.49, 0.24]
    colors = [GREEN, BLUE, PURPLE]
    names = [r"DPS fraction  $f_{DPS}$", r"DPS cross section  $\sigma_{DPS}$", r"Effective cross section  $\sigma_{eff}$"]
    values = [fdps_line(f), cross_section_line(d, r"\sigma_{DPS}", "pb"), cross_section_line(e, r"\sigma_{eff}", "mb")]
    notes = (
        [
            "Accepted fixed-500 statistical and three-group method uncertainty",
            r"Fiducial $m(J/\psi J/\psi)>7.5$ GeV; preliminary propagation",
            rf"No pair-mass cut after $\epsilon_m={RESULTS['mass_extrapolation']['epsilon_m']:.6f}$; preliminary propagation",
        ]
        if FROZEN_H016_REFERENCE
        else [
            "Fixed-500 statistical and three-group method uncertainty; candidate",
            r"Fiducial $m(J/\psi J/\psi)>7.5$ GeV; candidate propagation",
            f"No pair-mass cut after epsilon_m={RESULTS['mass_extrapolation']['epsilon_m']:.6f}; candidate propagation",
        ]
    )
    for yi, color, name, value, note in zip(y, colors, names, values, notes):
        ax.plot([0.055, 0.055], [yi - 0.085, yi + 0.085], color=color, linewidth=8,
                transform=ax.transAxes, solid_capstyle="round")
        ax.text(0.09, yi + 0.055, name, transform=ax.transAxes, fontsize=19,
                fontweight="bold", color=color, va="center")
        ax.text(0.50, yi + 0.01, value, transform=ax.transAxes, fontsize=22,
                fontweight="bold", ha="center", va="center", color=DARK)
        ax.text(0.09, yi - 0.065, note, transform=ax.transAxes, fontsize=14.5,
                color="#5b6470", va="center")
    ax.text(
        0.5, 0.055,
        "Cross-section results are preliminary" if FROZEN_H016_REFERENCE else "All results are candidates pending user confirmation",
        transform=ax.transAxes, ha="center", va="center", fontsize=14,
        bbox={"boxstyle": "round,pad=0.38", "facecolor": LIGHT, "edgecolor": "none"},
    )
    save(fig, "data_driven_result_summary")


def main() -> None:
    setup_style()
    plot_fdps()
    plot_sigma_dps()
    plot_sigma_eff()
    plot_summary()
    print("Created four H016 figures in PDF and PNG")


if __name__ == "__main__":
    main()
