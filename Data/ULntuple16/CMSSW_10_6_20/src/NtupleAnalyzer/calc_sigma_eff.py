#!/usr/bin/env python3
"""Calculate sigma_eff from the current AN inputs."""

import argparse
import math


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--single-j-times-br", type=float, default=10.63, help="nb")
    parser.add_argument("--single-j-stat", type=float, default=0.01, help="nb")
    parser.add_argument("--single-j-syst", type=float, default=0.09, help="nb")
    parser.add_argument("--br-mumu", type=float, default=0.05961)
    parser.add_argument("--br-mumu-error", type=float, default=0.00033)
    parser.add_argument("--sigma-dps", type=float, default=12.6, help="pb")
    parser.add_argument("--sigma-dps-stat", type=float, default=3.7, help="pb")
    parser.add_argument("--sigma-dps-syst-up", type=float, default=1.1, help="pb")
    parser.add_argument("--sigma-dps-syst-down", type=float, default=6.4, help="pb")
    parser.add_argument("--output", default="sigma_eff.txt")
    return parser.parse_args()


def main():
    args = parse_args()
    if args.br_mumu <= 0.0 or args.sigma_dps <= 0.0:
        raise ValueError("BR(J/psi -> mumu) and sigma_DPS must be positive.")

    # 0.001 converts nb^2/pb to mb.  The ATLAS input already contains
    # one BR(J/psi -> mumu), so the numerator requires division by BR^2.
    sigma_eff = (
        0.001
        * 0.5
        * args.single_j_times_br**2
        / (args.br_mumu**2 * args.sigma_dps)
    )

    single_j_stat = sigma_eff * 2.0 * args.single_j_stat / args.single_j_times_br
    single_j_syst = sigma_eff * 2.0 * args.single_j_syst / args.single_j_times_br
    br_syst = sigma_eff * 2.0 * args.br_mumu_error / args.br_mumu
    dps_stat = sigma_eff * args.sigma_dps_stat / args.sigma_dps

    # sigma_eff is inversely proportional to sigma_DPS, so directions swap.
    dps_syst_up = sigma_eff * args.sigma_dps_syst_down / args.sigma_dps
    dps_syst_down = sigma_eff * args.sigma_dps_syst_up / args.sigma_dps

    total_stat = math.hypot(single_j_stat, dps_stat)
    common_syst = math.hypot(single_j_syst, br_syst)
    total_syst_up = math.hypot(common_syst, dps_syst_up)
    total_syst_down = math.hypot(common_syst, dps_syst_down)

    lines = [
        f"single_j_times_br = {args.single_j_times_br:.6g} nb",
        f"br_mumu = {args.br_mumu:.6g} +/- {args.br_mumu_error:.6g}",
        f"sigma_DPS = {args.sigma_dps:.6g} +/- {args.sigma_dps_stat:.6g} (stat.) +{args.sigma_dps_syst_up:.6g}/-{args.sigma_dps_syst_down:.6g} (syst.) pb",
        f"sigma_eff = {sigma_eff:.6g} +/- {total_stat:.6g} (stat.) +{total_syst_up:.6g}/-{total_syst_down:.6g} (syst.) mb",
        "",
        "Uncertainty components (mb):",
        f"single_j_stat = {single_j_stat:.6g}",
        f"sigma_DPS_stat = {dps_stat:.6g}",
        f"single_j_syst = {single_j_syst:.6g}",
        f"br_syst = {br_syst:.6g}",
        f"sigma_DPS_syst_effect = +{dps_syst_up:.6g}/-{dps_syst_down:.6g}",
    ]
    output = "\n".join(lines) + "\n"
    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write(output)
    print(output, end="")
    print("Saved to", args.output)


if __name__ == "__main__":
    main()
