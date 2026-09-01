#!/usr/bin/env python3
"""Build H015 figure manifests, provenance metadata, and checksum index."""

from __future__ import annotations

import csv
import hashlib
import json
import struct
import subprocess
from datetime import datetime, timezone
from pathlib import Path


PACKAGE = Path(__file__).resolve().parents[1]
REPO = Path(__file__).resolve().parents[3]
FIGURES = PACKAGE / "figures"
DATA = PACKAGE / "data"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def png_dimensions(path: Path) -> tuple[int, int]:
    with path.open("rb") as handle:
        header = handle.read(24)
    if header[:8] != b"\x89PNG\r\n\x1a\n":
        raise RuntimeError(f"not a PNG: {path}")
    return struct.unpack(">II", header[16:24])


def pdf_info(path: Path) -> tuple[int, float, float, int]:
    output = subprocess.check_output(["pdfinfo", str(path)], text=True)
    pages = 0
    width = height = 0.0
    rotation = 0
    for line in output.splitlines():
        if line.startswith("Pages:"):
            pages = int(line.split(":", 1)[1])
        elif line.startswith("Page size:"):
            words = line.split()
            width, height = float(words[2]), float(words[4])
        elif line.startswith("Page rot:"):
            rotation = int(line.split(":", 1)[1]) % 360
    if rotation in (90, 270):
        width, height = height, width
    return pages, width, height, rotation


FIGURE_SPECS = [
    {
        "id": "F01",
        "stem": "pp_yield_map_current12",
        "title": "Prompt-prompt yield in Data",
        "sources": ["data/cells_current12.csv"],
        "definition": "N_PP and fitted uncertainty in each accepted cell",
        "caveat": "12/12 accepted cell fits; errors are fitted prompt-prompt uncertainties",
    },
    {
        "id": "F02",
        "stem": "control_region_geometry_current12",
        "title": "Data-driven control region in the 12-cell phase space",
        "sources": ["data/cells_current12.csv"],
        "definition": "Accepted 12-cell geometry with the single low-|Delta phi| control cell",
        "caveat": "Geometry only; the highlighted control cell supplies DPS normalization",
    },
    {
        "id": "F03",
        "stem": "dps_sps_decomposition_current12",
        "title": "CR-normalized data-derived DPS yield; Data-derived SPS residual",
        "sources": ["data/cells_current12.csv"],
        "definition": "Stored DPS/SPS yields and propagated per-cell uncertainties",
        "caveat": "Signed SPS residuals are retained without truncation",
    },
    {
        "id": "F04",
        "stem": "single_jpsi_fit_summary",
        "title": "Prompt single-J/psi extraction in Data",
        "sources": [
            "data/single_jpsi_fit_qa.csv",
            "Data_driven/reference/preupdates_20260829/h015/jpsi1_sweights.root",
            "Data_driven/reference/preupdates_20260829/h015/jpsi2_sweights.root",
            "Data_driven/reference/preupdates_20260829/h015/Model_4D_tot.root",
        ],
        "definition": "Unweighted mass/decay-length projections with saved two-stage yield-fit parameters",
        "caveat": "Projection only; no fitTo call and no fit result was regenerated",
    },
    {
        "id": "F05",
        "stem": "single_jpsi_category_dependence",
        "title": "Single-J/psi category dependence",
        "sources": ["data/single_jpsi_category_dependence.csv"],
        "definition": "Accepted central value and pT/|y| category variations",
        "caveat": "Alternatives are variations; neither is assigned a method-status label",
    },
    {
        "id": "F06",
        "stem": "bootstrap_fdps_fixed500",
        "title": "Full-workflow bootstrap distribution",
        "sources": ["data/bootstrap_fdps_fixed500.csv", "data/bootstrap_plot_metadata.json"],
        "definition": "Accepted fixed 500/500 full-workflow replica f_DPS distribution",
        "caveat": "No failed or replacement attempts; 20 fixed bins over [0.02,0.40]",
    },
    {
        "id": "F07",
        "stem": "control_region_variations",
        "title": "Control-region definition and purity variations",
        "sources": ["data/control_region_variations.csv"],
        "definition": "Accepted CR-boundary and zero/MC-leakage f_DPS variations",
        "caveat": "Envelope is the accepted group envelope, not an additional fit",
    },
    {
        "id": "F08",
        "stem": "pure_dps_full_workflow_closure",
        "title": "Pure-DPS full-workflow closure",
        "sources": ["data/pure_dps_full_workflow_closure.json"],
        "definition": "Expected/recovered pure-DPS fraction and accepted mapping to Data uncertainty",
        "caveat": "Conditional on the current DPS simulation model being treated as correct",
    },
    {
        "id": "F09",
        "stem": "method_systematic_breakdown",
        "title": "Data-driven method uncertainty",
        "sources": ["data/method_systematic_breakdown.csv"],
        "definition": "Three included groups and separate upward/downward quadrature",
        "caveat": "Deferred and diagnostic sources are intentionally absent, not set to zero",
    },
    {
        "id": "F10",
        "stem": "sps_significance_map_current12",
        "title": "SPS residual significance by cell",
        "sources": ["data/cells_current12.csv"],
        "definition": "Stored signed per-cell SPS residual significance",
        "caveat": "Per-cell diagnostic only; not a global significance or extra systematic",
    },
]


def resolve_source(source: str) -> Path:
    path = Path(source)
    return PACKAGE / path if source.startswith("data/") else REPO / path


def build_manifest() -> None:
    config = json.loads((DATA / "figure_config.json").read_text(encoding="utf-8"))
    now = datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")
    rows: list[dict[str, str]] = []
    for spec in FIGURE_SPECS:
        pdf = FIGURES / f"{spec['stem']}.pdf"
        png = FIGURES / f"{spec['stem']}.png"
        if not pdf.is_file() or not png.is_file():
            raise RuntimeError(f"missing figure pair for {spec['id']}")
        pages, width, height, rotation = pdf_info(pdf)
        png_width, png_height = png_dimensions(png)
        if pages != 1 or abs(width / height - 16.0 / 9.0) > 0.015:
            raise RuntimeError(f"{spec['id']} PDF is not a single 16:9 page: {width}x{height}")
        if abs(png_width / png_height - 16.0 / 9.0) > 0.015:
            raise RuntimeError(f"{spec['id']} PNG is not 16:9: {png_width}x{png_height}")
        source_paths = [resolve_source(source) for source in spec["sources"]]
        if not all(path.is_file() for path in source_paths):
            raise RuntimeError(f"missing machine source for {spec['id']}")
        rows.append(
            {
                "figure_id": spec["id"],
                "title": spec["title"],
                "status": "complete",
                "pdf": str(pdf.relative_to(PACKAGE)),
                "pdf_sha256": sha256(pdf),
                "pdf_page_points": f"{width:g}x{height:g}",
                "pdf_page_rotation": str(rotation),
                "png": str(png.relative_to(PACKAGE)),
                "png_sha256": sha256(png),
                "png_pixels": f"{png_width}x{png_height}",
                "machine_sources": ";".join(spec["sources"]),
                "machine_source_sha256": ";".join(sha256(path) for path in source_paths),
                "selection": config["physics_selection"],
                "normalization_weight_error_definition": spec["definition"],
                "generation_script": (
                    "scripts/plot_single_jpsi_fit_summary.cpp"
                    if spec["id"] == "F04"
                    else "scripts/plot_figures.py"
                ),
                "reproduction_command": "bash scripts/run_all.sh",
                "producer_version": f"{config['producer_branch']}@{config['producer_head']}",
                "generated_at": now,
                "caveat": spec["caveat"],
            }
        )
    fieldnames = list(rows[0])
    with (PACKAGE / "figure_manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    lines = ["# H015 figure manifest", "", "Status: `ready` for consumer review.", ""]
    for row in rows:
        lines.extend(
            [
                f"## {row['figure_id']} - {row['title']}",
                "",
                f"- PDF: `{row['pdf']}` (`{row['pdf_sha256']}`)",
                f"- PNG: `{row['png']}` (`{row['png_sha256']}`)",
                f"- Machine source: `{row['machine_sources']}`",
                f"- Definition: {row['normalization_weight_error_definition']}",
                f"- Caveat: {row['caveat']}",
                "",
            ]
        )
    (PACKAGE / "figure_manifest.md").write_text("\n".join(lines), encoding="utf-8")


def build_metadata_and_checksums() -> None:
    branch = subprocess.check_output(["git", "rev-parse", "--abbrev-ref", "HEAD"], cwd=REPO, text=True).strip()
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip()
    dirty = subprocess.check_output(
        ["git", "status", "--porcelain", "--untracked-files=all"], cwd=REPO, text=True
    ).splitlines()
    metadata = {
        "artifact": "H015_data_driven_slide_figure_package",
        "status": "ready",
        "consumer_status": "pending",
        "producer_branch": branch,
        "producer_head": head,
        "accepted_physics_source": "H014 accepted package",
        "legacy_manifest_sha256": sha256(
            REPO / "Data_driven/reference/preupdates_20260829/h015/figure_manifest.csv"
        ),
        "analysis_fit_bootstrap_systematic_rerun": False,
        "F04_runtime": "CMSSW_10_6_20 / ROOT 6.14/09 / cmssw-el7",
        "F04_operation": "plot projections from saved global_yield_fit_result; no fitTo call",
        "figures_complete": "10/10 PDF and 10/10 PNG",
        "tracked_dirty_files": dirty,
        "dirty_file_ownership": "generated presentation outputs are ignored by repository rules",
        "generated_at": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
        "real_newlines": True,
    }
    (PACKAGE / "package_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    checksum_path = PACKAGE / "package_checksums.sha256"
    files = sorted(
        path for path in PACKAGE.rglob("*") if path.is_file() and path != checksum_path
    )
    checksum_path.write_text(
        "".join(f"{sha256(path)}  {path.relative_to(PACKAGE)}\n" for path in files),
        encoding="utf-8",
    )


def main() -> None:
    build_manifest()
    build_metadata_and_checksums()
    print("Built H015 ready-package manifests and checksums")


if __name__ == "__main__":
    main()
