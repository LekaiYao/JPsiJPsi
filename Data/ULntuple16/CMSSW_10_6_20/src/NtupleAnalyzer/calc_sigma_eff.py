#!/usr/bin/env python3
import math
import re


unit = 1000.0  # convert nb^2/pb to mb

# 1) Hard-coded total cross section and uncertainties, unit pb
sigma_tot = 39.2
sigma_tot_stat = 0.6
sigma_tot_syst = 6.1

# 2) Hard-coded input file from Tpl3_Fit.cpp output
fractions_file = "fig/temp3/fractions3_summary.txt"

# 4) Hard-coded single-J/psi cross section(unit nb) and fiducial efficiency correction
# Assume no uncertainty on these two values for now.
sigma_J = 268.0
eff_fid = 0.9862


def read_f_dps_values(path):
    with open(path, "r") as f:
        content = f.read()

    def extract(name):
        m = re.search(rf"^{re.escape(name)}\s*=\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*$", content, re.MULTILINE)
        if not m:
            raise ValueError(f"Cannot find '{name}' in {path}")
        return float(m.group(1))

    f_dps = extract("2D f_DPS")
    f_dps_stat = extract("2D f_DPS_stat_err")
    f_dps_syst = extract("2D f_DPS_sys_err")
    return f_dps, f_dps_stat, f_dps_syst


def main():
    # 2) read f_DPS and its stat/syst errors
    f_dps, f_dps_stat, f_dps_syst = read_f_dps_values(fractions_file)

    # 3) sigma_DPS and uncertainties
    sigma_dps = sigma_tot * f_dps
    sigma_dps_stat = math.sqrt((f_dps * sigma_tot_stat) ** 2 + (sigma_tot * f_dps_stat) ** 2)
    sigma_dps_syst = math.sqrt((f_dps * sigma_tot_syst) ** 2 + (sigma_tot * f_dps_syst) ** 2)

    if sigma_dps == 0:
        raise ZeroDivisionError("sigma_DPS is zero, cannot compute sigma_eff.")

    # 5) sigma_eff and uncertainties
    sigma_eff = (sigma_J * sigma_J) / (2.0 * sigma_dps * eff_fid * unit)
    sigma_eff_stat = abs(sigma_eff) * abs(sigma_dps_stat / sigma_dps)
    sigma_eff_syst = abs(sigma_eff) * abs(sigma_dps_syst / sigma_dps)

    # 6) output to sigma_eff.txt
    out_file = "sigma_eff.txt"
    with open(out_file, "w") as f:
        f.write(f"sigma_tot = {sigma_tot} pb\n")
        f.write(f"sigma_tot_stat = {sigma_tot_stat} pb\n")
        f.write(f"sigma_tot_syst = {sigma_tot_syst} pb\n")
        f.write(f"f_DPS = {f_dps}\n")
        f.write(f"f_DPS_stat = {f_dps_stat}\n")
        f.write(f"f_DPS_syst = {f_dps_syst}\n")
        f.write(f"sigma_DPS = {sigma_dps} pb\n")
        f.write(f"sigma_DPS_stat = {sigma_dps_stat} pb\n")
        f.write(f"sigma_DPS_syst = {sigma_dps_syst} pb\n")
        f.write(f"sigma_J = {sigma_J} nb\n")
        f.write(f"eff_fid = {eff_fid}\n")
        f.write(f"sigma_eff = {sigma_eff} mb\n")
        f.write(f"sigma_eff_stat = {sigma_eff_stat} mb\n")
        f.write(f"sigma_eff_syst = {sigma_eff_syst} mb\n")

    print("Saved to", out_file)


if __name__ == "__main__":
    main()
