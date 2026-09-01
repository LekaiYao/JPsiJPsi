#!/usr/bin/env python3
"""Create H015 presentation figures from frozen CSV/JSON inputs."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


PACKAGE = Path(__file__).resolve().parents[1]
DATA = PACKAGE / "data"
FIGURES = PACKAGE / "figures"


def read_csv(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def setup_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 16,
            "axes.titlesize": 24,
            "axes.labelsize": 20,
            "xtick.labelsize": 16,
            "ytick.labelsize": 16,
            "legend.fontsize": 14,
            "figure.titlesize": 26,
            "axes.linewidth": 1.4,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save(fig: plt.Figure, stem: str) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES / f"{stem}.pdf", facecolor="white")
    fig.savefig(FIGURES / f"{stem}.png", dpi=180, facecolor="white")
    plt.close(fig)


def float_cells() -> list[dict[str, float]]:
    numeric = {
        "cell",
        "dy_low",
        "dy_high",
        "dphi_low",
        "dphi_high",
        "in_control",
        "data",
        "data_error",
        "dps",
        "dps_error",
        "sps",
        "sps_error",
        "sps_significance",
    }
    result: list[dict[str, float]] = []
    for row in read_csv("cells_current12.csv"):
        result.append({key: float(row[key]) for key in numeric})
    return result


def configure_phase_space(ax: plt.Axes) -> None:
    ax.set_xlim(0.0, 4.0)
    ax.set_ylim(0.0, np.pi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"$|\Delta y|$")
    ax.set_ylabel(r"$|\Delta \phi|$")
    ax.set_xticks([0.0, 0.6, 1.2, 1.8, 4.0])
    ax.set_xticklabels(["0", "0.6", "1.2", "1.8", "4.0"])
    ax.set_yticks([0.0, np.pi / 2.0, 3.0 * np.pi / 4.0, np.pi])
    ax.set_yticklabels(["0", r"$\pi/2$", r"$3\pi/4$", r"$\pi$"])


def draw_cell_map(
    ax: plt.Axes,
    cells: list[dict[str, float]],
    value_key: str,
    error_key: str | None,
    cmap_name: str,
    norm: Normalize,
    annotate_format: str,
    text_size: float = 12,
) -> mpl.cm.ScalarMappable:
    cmap = plt.get_cmap(cmap_name)
    for row in cells:
        x0, x1 = row["dy_low"], row["dy_high"]
        y0, y1 = row["dphi_low"], row["dphi_high"]
        value = row[value_key]
        face = cmap(norm(value))
        rect = Rectangle(
            (x0, y0), x1 - x0, y1 - y0, facecolor=face, edgecolor="black", linewidth=1.6
        )
        ax.add_patch(rect)
        rgb = mpl.colors.to_rgb(face)
        luminance = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]
        color = "black" if luminance > 0.55 else "white"
        if error_key:
            label = f"{format(value, annotate_format)}\n+/-{format(row[error_key], annotate_format)}"
        else:
            label = format(value, annotate_format)
        ax.text(
            0.5 * (x0 + x1),
            0.5 * (y0 + y1),
            label,
            ha="center",
            va="center",
            color=color,
            fontsize=text_size,
            fontweight="semibold",
        )
    configure_phase_space(ax)
    return mpl.cm.ScalarMappable(norm=norm, cmap=cmap)


def plot_f01(cells: list[dict[str, float]]) -> None:
    fig, ax = plt.subplots(figsize=(13.333, 7.5), constrained_layout=True)
    values = np.array([row["data"] for row in cells])
    sm = draw_cell_map(ax, cells, "data", "data_error", "YlGnBu", Normalize(values.min(), values.max()), ".1f")
    ax.set_title("Prompt-prompt yield in Data", pad=14)
    cbar = fig.colorbar(sm, ax=ax, pad=0.03, fraction=0.045)
    cbar.set_label(r"Fitted $N_{PP}^{i}$")
    ax.text(
        0.99,
        0.015,
        "12/12 cell fits accepted",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=14,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.9},
    )
    save(fig, "pp_yield_map_current12")


def plot_f02(cells: list[dict[str, float]]) -> None:
    fig, ax = plt.subplots(figsize=(13.333, 7.5), constrained_layout=True)
    for row in cells:
        x0, x1 = row["dy_low"], row["dy_high"]
        y0, y1 = row["dphi_low"], row["dphi_high"]
        in_control = bool(row["in_control"])
        rect = Rectangle(
            (x0, y0),
            x1 - x0,
            y1 - y0,
            facecolor="#f4a261" if in_control else "#e9ecef",
            edgecolor="#9c3d10" if in_control else "#495057",
            linewidth=3.0 if in_control else 1.8,
        )
        ax.add_patch(rect)
        if in_control:
            ax.text(
                0.5 * (x0 + x1),
                0.5 * (y0 + y1),
                "Control region\nDPS normalization",
                ha="center",
                va="center",
                fontsize=18,
                fontweight="bold",
                color="#6a2c00",
            )
    configure_phase_space(ax)
    ax.set_title("Data-driven control region in the 12-cell phase space", pad=14)
    ax.text(
        0.02,
        0.98,
        r"$1.8 \leq |\Delta y| < 4.0$ and $|\Delta \phi| \leq \pi/2$",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=17,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "alpha": 0.92},
    )
    save(fig, "control_region_geometry_current12")


def plot_f03(cells: list[dict[str, float]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13.333, 7.5), constrained_layout=True)
    specs = [
        ("dps", "dps_error", "CR-normalized data-derived DPS yield", "YlGnBu", False),
        ("sps", "sps_error", "Data-derived SPS residual", "RdBu_r", True),
    ]
    for ax, (value, error, title, cmap, diverging) in zip(axes, specs):
        values = np.array([row[value] for row in cells])
        if diverging:
            extent = max(abs(values.min()), abs(values.max()))
            norm: Normalize = TwoSlopeNorm(vmin=-extent, vcenter=0.0, vmax=extent)
        else:
            norm = Normalize(values.min(), values.max())
        sm = draw_cell_map(ax, cells, value, error, cmap, norm, ".1f", 10.5)
        display_title = "CR-normalized data-derived\nDPS yield" if value == "dps" else title
        ax.set_title(display_title, pad=12, fontsize=18)
        cbar = fig.colorbar(sm, ax=ax, pad=0.025, fraction=0.05)
        cbar.set_label("Yield")
    save(fig, "dps_sps_decomposition_current12")


def plot_f05() -> None:
    rows = read_csv("single_jpsi_category_dependence.csv")
    labels = [row["label"] for row in rows]
    values = np.array([float(row["f_dps"]) for row in rows])
    shifts = np.array([float(row["signed_shift"]) for row in rows])
    central = values[0]
    fig, ax = plt.subplots(figsize=(13.333, 7.5), constrained_layout=True)
    y = np.arange(len(rows))[::-1]
    ax.axvspan(central - 0.006246760932, central + 0.000778391595, color="#457b9d", alpha=0.16)
    ax.axvline(central, color="#1d3557", linewidth=2.0, linestyle="--", label="Central value")
    ax.scatter(values, y, s=130, color=["#1d3557", "#e76f51", "#2a9d8f"], zorder=3)
    for yi, value, shift in zip(y, values, shifts):
        ax.text(value + 0.00025, yi, f"{value:.6f}  ({shift:+.6f})", va="center", fontsize=16)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel(r"$f_{DPS}$  (signed shift in parentheses)")
    ax.set_title("Single-J/psi category dependence", pad=34)
    ax.grid(axis="x", alpha=0.25)
    ax.legend(loc="lower right")
    ax.set_xlim(values.min() - 0.0025, values.max() + 0.0055)
    ax.text(
        0.50,
        1.01,
        r"Category envelope: $+0.000778/-0.006247$",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=17,
    )
    save(fig, "single_jpsi_category_dependence")


def plot_f06() -> None:
    rows = read_csv("bootstrap_fdps_fixed500.csv")
    values = np.array([float(row["f_dps"]) for row in rows])
    meta = json.loads((DATA / "bootstrap_plot_metadata.json").read_text(encoding="utf-8"))
    edges = np.linspace(meta["histogram_range"][0], meta["histogram_range"][1], meta["histogram_bins"] + 1)
    fig, ax = plt.subplots(figsize=(13.333, 7.5), constrained_layout=True)
    ax.hist(values, bins=edges, color="#5e81ac", edgecolor="white", linewidth=1.2)
    ax.axvspan(meta["q16"], meta["q84"], color="#88c0d0", alpha=0.24, label="16%-84% interval")
    ax.axvline(meta["central"], color="#bf616a", linewidth=2.4, label=f"Central = {meta['central']:.6f}")
    ax.axvline(meta["mean"], color="#2e3440", linewidth=2.2, linestyle="--", label=f"Mean = {meta['mean']:.6f}")
    ax.set_xlabel(r"Bootstrap $f_{DPS}$")
    ax.set_ylabel("Replicas")
    ax.set_title("Full-workflow bootstrap distribution", pad=14)
    ax.grid(axis="y", alpha=0.25)
    ax.legend(loc="upper right")
    ax.text(
        0.02,
        0.95,
        f"Sample SD = {meta['sample_sd']:.6f}\n500/500 successful replicas; no replacement attempts",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=17,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "alpha": 0.92},
    )
    save(fig, "bootstrap_fdps_fixed500")


def plot_f07() -> None:
    rows = read_csv("control_region_variations.csv")
    central = float(rows[0]["f_dps"])
    fig, ax = plt.subplots(figsize=(13.333, 7.5), constrained_layout=True)
    y = np.arange(len(rows))[::-1]
    styles = {
        "zero leakage": ("o", "#1d3557"),
        "MC leakage": ("s", "#e76f51"),
        "zero or MC leakage": ("D", "#2a9d8f"),
    }
    ax.axvspan(central - 0.006494652437, central + 0.007830845807, color="#457b9d", alpha=0.14)
    ax.axvline(central, color="black", linewidth=1.8, linestyle="--")
    for yi, row in zip(y, rows):
        marker, color = styles[row["leakage_model"]]
        value = float(row["f_dps"])
        shift = float(row["signed_shift"])
        ax.scatter(value, yi, marker=marker, color=color, s=115, zorder=3)
        ax.text(value + 0.00030, yi, f"{value:.6f}  ({shift:+.6f})", va="center", fontsize=14)
    ax.set_yticks(y)
    display_labels = {
        "|Delta phi| <= 70 deg": r"$|\Delta \phi| \leq 70^{\circ}$",
        "|Delta phi| <= 110 deg": r"$|\Delta \phi| \leq 110^{\circ}$",
    }
    ax.set_yticklabels([display_labels.get(row["label"], row["label"]) for row in rows])
    ax.set_xlabel(r"$f_{DPS}$  (signed shift in parentheses)")
    ax.set_title("Control-region definition and purity variations", pad=34)
    ax.grid(axis="x", alpha=0.25)
    handles = [
        Line2D([0], [0], marker=marker, color="none", markerfacecolor=color, markeredgecolor=color, markersize=10, label=label)
        for label, (marker, color) in styles.items()
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=3)
    ax.set_xlim(min(float(row["f_dps"]) for row in rows) - 0.004, max(float(row["f_dps"]) for row in rows) + 0.007)
    ax.text(
        0.50,
        1.01,
        r"Group envelope: $+0.007831/-0.006495$",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=17,
    )
    save(fig, "control_region_variations")


def plot_f08() -> None:
    values = json.loads((DATA / "pure_dps_full_workflow_closure.json").read_text(encoding="utf-8"))
    fig, ax = plt.subplots(figsize=(13.333, 7.5), constrained_layout=True)
    x = np.array([0.0, 1.0])
    y = np.array([values["expected_f_dps"], values["recovered_f_dps"]])
    ax.plot(x, y, color="#6c757d", linewidth=2.0, zorder=1)
    ax.scatter(x, y, s=260, color=["#457b9d", "#e76f51"], zorder=3)
    ax.axhline(1.0, color="#457b9d", linestyle="--", linewidth=1.8)
    ax.set_xticks(x)
    ax.set_xticklabels(["Expected pure DPS", "Recovered by full workflow"])
    ax.set_ylabel(r"$f_{DPS}$")
    ax.set_ylim(0.96, 1.065)
    ax.set_title("Pure-DPS full-workflow closure", pad=14)
    ax.grid(axis="y", alpha=0.25)
    ax.text(0.0, y[0] + 0.006, f"{y[0]:.6f}", ha="center", fontsize=18, fontweight="bold")
    ax.text(1.0, y[1] + 0.006, f"{y[1]:.6f}", ha="center", fontsize=18, fontweight="bold")
    mapping = (
        f"Relative closure bias = {values['relative_closure_bias']:.6f}\n"
        f"Mapped Data uncertainty = {values['data_f_dps']:.6f} x {values['relative_closure_bias']:.6f}"
        f" = {values['mapped_absolute_uncertainty']:.6f}"
    )
    ax.text(
        0.50,
        0.08,
        mapping,
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=17,
        bbox={"boxstyle": "round,pad=0.45", "facecolor": "white", "edgecolor": "#6c757d"},
    )
    ax.text(
        0.01,
        0.01,
        "Conditional on the current DPS simulation model.",
        transform=ax.transAxes,
        fontsize=13,
        color="#495057",
    )
    save(fig, "pure_dps_full_workflow_closure")


def plot_f09() -> None:
    rows = read_csv("method_systematic_breakdown.csv")
    labels = [row["group"] for row in rows]
    up = np.array([float(row["up"]) for row in rows])
    down = np.array([float(row["down"]) for row in rows])
    y = np.arange(len(rows))[::-1]
    colors = ["#457b9d", "#2a9d8f", "#e76f51", "#1d3557"]
    fig, ax = plt.subplots(figsize=(13.333, 7.5), constrained_layout=True)
    ax.axvline(0.0, color="black", linewidth=1.6)
    for index, (yi, left, right, color) in enumerate(zip(y, down, up, colors)):
        linewidth = 5.0 if index == len(rows) - 1 else 3.2
        ax.plot([-left, right], [yi, yi], color=color, linewidth=linewidth, solid_capstyle="round")
        ax.scatter([-left, right], [yi, yi], color=color, s=75, zorder=3)
        ax.text(right + 0.00035, yi, f"+{right:.6f}", va="center", fontsize=15)
        ax.text(-left - 0.00035, yi, f"-{left:.6f}", va="center", ha="right", fontsize=15)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel(r"Signed uncertainty on $f_{DPS}$")
    ax.set_title("Data-driven method uncertainty", pad=14)
    ax.grid(axis="x", alpha=0.25)
    ax.set_xlim(-0.0145, 0.0135)
    ax.text(
        0.98,
        0.43,
        "Combined separately in upward/downward quadrature",
        transform=ax.transAxes,
        ha="right",
        fontsize=16,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.9},
    )
    save(fig, "method_systematic_breakdown")


def plot_f10(cells: list[dict[str, float]]) -> None:
    fig, ax = plt.subplots(figsize=(13.333, 7.5), constrained_layout=True)
    values = np.array([row["sps_significance"] for row in cells])
    extent = max(abs(values.min()), abs(values.max()))
    sm = draw_cell_map(
        ax,
        cells,
        "sps_significance",
        None,
        "RdBu_r",
        TwoSlopeNorm(vmin=-extent, vcenter=0.0, vmax=extent),
        ".2f",
    )
    ax.set_title("SPS residual significance by cell", pad=14)
    cbar = fig.colorbar(sm, ax=ax, pad=0.03, fraction=0.045)
    cbar.set_label(r"$N_{SPS}/\sigma(N_{SPS})$")
    ax.text(
        0.99,
        0.015,
        "Per-cell diagnostic only",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=14,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.9},
    )
    save(fig, "sps_significance_map_current12")


def main() -> None:
    setup_style()
    cells = float_cells()
    plot_f01(cells)
    plot_f02(cells)
    plot_f03(cells)
    plot_f05()
    plot_f06()
    plot_f07()
    plot_f08()
    plot_f09()
    plot_f10(cells)
    print("Created F01-F03 and F05-F10 in PDF and PNG")


if __name__ == "__main__":
    main()
