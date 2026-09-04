#!/usr/bin/env python3
"""Build a temporary pair cross-section input from a total-only 4D fit."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from pathlib import Path


EXPECTED = {
    "total_summary":
        "37f415988f93ffcafde46553d791776d1b831456d18ee16d0da68a4315cbc2ee",
    "total_metadata":
        "ef1520c2ecb9fc7851f277f14e907bf331f8f213a4fdfb552dcf82aea97811dc",
    "sys3":
        "eabea32c595a895a599632967a1b401aafce7cab083c81515ed47889b498e63d",
}
LUMINOSITY_FB = 36.684
BRANCHING_FRACTION = 0.05961


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    if len(sys.argv) != 5:
        raise SystemExit(
            "usage: build_total_only_pair_cross_section.py "
            "TOTAL_SUMMARY TOTAL_METADATA SYS3_CSV OUTPUT_JSON"
        )
    summary_path, metadata_path, sys3_path, output_path = map(
        Path, sys.argv[1:]
    )
    if output_path.exists():
        raise RuntimeError(f"refusing to overwrite {output_path}")
    for label, path in (
        ("total_summary", summary_path),
        ("total_metadata", metadata_path),
        ("sys3", sys3_path),
    ):
        if sha256(path) != EXPECTED[label]:
            raise RuntimeError(f"{label} checksum mismatch: {path}")

    with summary_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1 or rows[0]["fit"] != "nominal":
        raise RuntimeError("unexpected total fit summary")
    row = rows[0]
    if int(row["status"]) != 0 or int(row["covQual"]) != 3:
        raise RuntimeError("total fit failed QA")
    if float(row["edm"]) >= 0.01:
        raise RuntimeError("total fit EDM failed QA")
    n_pp = float(row["n_P_P"])
    n_pp_error = float(row["n_P_P_error"])

    with sys3_path.open(newline="", encoding="utf-8") as handle:
        total_rows = [
            item for item in csv.DictReader(handle)
            if item["scope"] == "total"
        ]
    if len(total_rows) != 1:
        raise RuntimeError("sys3 CSV does not contain exactly one total row")
    correction_relative = float(total_rows[0]["sys3"])

    denominator = (
        LUMINOSITY_FB * 1000.0 * BRANCHING_FRACTION * BRANCHING_FRACTION
    )
    central_pb = n_pp / denominator
    stat_pb = n_pp_error / denominator
    sources = {
        "branching_fraction": {
            "relative": 0.011,
            "label": "Branching fraction",
        },
        "cms_luminosity": {
            "relative": 0.012,
            "label": "CMS integrated luminosity",
        },
        "acceptance_efficiency": {
            "relative": correction_relative,
            "label": "0.85:0.15 acceptance and efficiency correction",
        },
        "lifetime_variable": {
            "relative": 0.003,
            "label": "Lifetime variable",
        },
    }
    relative_systematic = math.sqrt(
        sum(source["relative"] ** 2 for source in sources.values())
    )
    result = {
        "central_pb": central_pb,
        "stat_pb": stat_pb,
        "included_syst_pb": central_pb * relative_systematic,
        "reported_syst_pb": None,
        "systematic_sources_relative": sources,
        "pending_systematic_sources": {
            "fitter_stability": {
                "relative": None,
                "status": "pending_not_included_not_zero",
                "reason": (
                    "the transferred total-only generation has no accepted "
                    "nominal/reference fitter variation or Data cross-section result"
                ),
            }
        },
        "systematic_combination": (
            "Known numeric pair sources combined in quadrature; fitter stability "
            "remains pending and is not treated as zero."
        ),
        "reconstructed_relative_systematic": relative_systematic,
        "reconstructed_absolute_systematic_pb":
            central_pb * relative_systematic,
        "phase_space": (
            "prompt Jpsi pair; each Jpsi 10<=pT<=40 GeV and abs(y)<2.0; "
            "m(JJ)>=7.5 GeV"
        ),
        "visible_convention": (
            "physical pair cross section after division by B(Jpsi->mumu)^2"
        ),
        "derivation": (
            "n_P_P/(luminosity_fb*1000*B(Jpsi->mumu)^2) from the transferred "
            "weighted total-only 4D fit"
        ),
        "n_P_P": n_pp,
        "n_P_P_error": n_pp_error,
        "luminosity_fb": LUMINOSITY_FB,
        "branching_fraction": BRANCHING_FRACTION,
        "source": (
            "0.85:0.15 pair-symmetrized total-only fit; temporary Data_driven "
            "candidate pending the Data cross-section campaign"
        ),
        "source_files": [
            str(summary_path.resolve()),
            str(metadata_path.resolve()),
            str(sys3_path.resolve()),
        ],
        "source_sha256": {
            "total_summary": sha256(summary_path),
            "total_metadata": sha256(metadata_path),
            "sys3": sha256(sys3_path),
        },
        "status": (
            "temporary_total_only_pending_Data_cross_section_and_user_confirmation"
        ),
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        f"PAIR_TOTAL_ONLY_CANDIDATE central_pb={central_pb:.12g} "
        f"stat_pb={stat_pb:.12g} included_syst_pb="
        f"{central_pb * relative_systematic:.12g}"
    )


if __name__ == "__main__":
    main()
