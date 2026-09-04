#!/usr/bin/env python3
"""Generate one direct H018 total-fit source from the frozen nominal source."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


BASE_SHA256 = "9cccff5a8235a576a9fadc20dc936bdbb1ca1cae8f0dd1f60f7ae6916ab812c1"
VARIATIONS = {
    "f1_double_gaussian",
    "f1_double_crystal_ball",
    "f1_triple_gaussian",
    "f2_linear",
    "f3_single_gaussian",
    "f3_triple_gaussian",
    "f4_exp_double_gaussian",
    "f4_double_exp_gaussian",
}


def replace_once(source: str, old: str, new: str) -> str:
    if source.count(old) != 1:
        raise RuntimeError(f"expected one exact replacement target, found {source.count(old)}")
    return source.replace(old, new, 1)


def replace_block(source: str, start: str, end: str, replacement: str) -> str:
    start_index = source.find(start)
    if start_index < 0:
        raise RuntimeError(f"missing block start: {start.strip()}")
    end_index = source.find(end, start_index)
    if end_index < 0:
        raise RuntimeError(f"missing block end: {end.strip()}")
    end_index += len(end)
    return source[:start_index] + replacement + source[end_index:]


def add_freeze(source: str, marker: str, additions: str) -> str:
    return replace_once(source, marker, marker + additions)


def generate(source: str, variation: str) -> str:
    if variation == "f1_double_gaussian":
        # The frozen source already implements this exact legacy variation as isRef=true.
        return source

    if variation == "f1_double_crystal_ball":
        source = replace_block(
            source,
            '    RooRealVar Jpsi_devia2("Jpsi_devia2"',
            '    RooAddPdf JpsiMassSig2("JpsiMassSig2", "JpsiMassSig2", RooArgList(Jpsi_core_2, Jpsi_gaussian_2), Jpsi_ratio);\n',
            '''    RooRealVar Jpsi_devia2("Jpsi_devia2", "Jpsi_devia2", 0.05, 0.02, 0.08);
    RooRealVar Jpsi_alpha2("Jpsi_alpha2", "Jpsi_alpha2", 1.5, 0.1, 3.5);
    RooRealVar Jpsi_nx2("Jpsi_nx2", "Jpsi_nx2", 1, 0, 100);
    RooRealVar Jpsi_ratio("Jpsi_ratio", "Jpsi_ratio", 0.6, 0, 1);
    RooCBShape Jpsi_crysBall_1("Jpsi_crysBall_1", "Jpsi_crysBall_1", Jpsi_mass1, Jpsi_mean, Jpsi_devia1, Jpsi_alpha1, Jpsi_nx1);
    RooCBShape Jpsi_crysBallB_1("Jpsi_crysBallB_1", "Jpsi_crysBallB_1", Jpsi_mass1, Jpsi_mean, Jpsi_devia2, Jpsi_alpha2, Jpsi_nx2);
    RooCBShape Jpsi_crysBall_2("Jpsi_crysBall_2", "Jpsi_crysBall_2", Jpsi_mass2, Jpsi_mean, Jpsi_devia1, Jpsi_alpha1, Jpsi_nx1);
    RooCBShape Jpsi_crysBallB_2("Jpsi_crysBallB_2", "Jpsi_crysBallB_2", Jpsi_mass2, Jpsi_mean, Jpsi_devia2, Jpsi_alpha2, Jpsi_nx2);
    RooAddPdf JpsiMassSig1("JpsiMassSig1", "JpsiMassSig1", RooArgList(Jpsi_crysBall_1, Jpsi_crysBallB_1), Jpsi_ratio);
    RooAddPdf JpsiMassSig2("JpsiMassSig2", "JpsiMassSig2", RooArgList(Jpsi_crysBall_2, Jpsi_crysBallB_2), Jpsi_ratio);
''',
        )
        source = add_freeze(
            source,
            "    Jpsi_nx1.setConstant(kTRUE);\n",
            "    Jpsi_alpha2.setConstant(kTRUE);\n    Jpsi_nx2.setConstant(kTRUE);\n",
        )
    elif variation == "f1_triple_gaussian":
        source = replace_block(
            source,
            '    RooRealVar Jpsi_devia2("Jpsi_devia2"',
            '    RooAddPdf JpsiMassSig2("JpsiMassSig2", "JpsiMassSig2", RooArgList(Jpsi_core_2, Jpsi_gaussian_2), Jpsi_ratio);\n',
            '''    RooRealVar Jpsi_devia2("Jpsi_devia2", "Jpsi_devia2", 0.05, 0.02, 0.08);
    RooRealVar Jpsi_devia3("Jpsi_devia3", "Jpsi_devia3", 0.06, 0.001, 0.1);
    RooRealVar Jpsi_ratio("Jpsi_ratio", "Jpsi_ratio", 0.6, 0, 1);
    RooRealVar Jpsi_ratio2("Jpsi_ratio2", "Jpsi_ratio2", 0.3, 0, 1);
    RooGaussian Jpsi_gaussianA_1("Jpsi_gaussianA_1", "Jpsi_gaussianA_1", Jpsi_mass1, Jpsi_mean, Jpsi_devia1);
    RooGaussian Jpsi_gaussianB_1("Jpsi_gaussianB_1", "Jpsi_gaussianB_1", Jpsi_mass1, Jpsi_mean, Jpsi_devia2);
    RooGaussian Jpsi_gaussianC_1("Jpsi_gaussianC_1", "Jpsi_gaussianC_1", Jpsi_mass1, Jpsi_mean, Jpsi_devia3);
    RooGaussian Jpsi_gaussianA_2("Jpsi_gaussianA_2", "Jpsi_gaussianA_2", Jpsi_mass2, Jpsi_mean, Jpsi_devia1);
    RooGaussian Jpsi_gaussianB_2("Jpsi_gaussianB_2", "Jpsi_gaussianB_2", Jpsi_mass2, Jpsi_mean, Jpsi_devia2);
    RooGaussian Jpsi_gaussianC_2("Jpsi_gaussianC_2", "Jpsi_gaussianC_2", Jpsi_mass2, Jpsi_mean, Jpsi_devia3);
    RooAddPdf JpsiMassSig1("JpsiMassSig1", "JpsiMassSig1", RooArgList(Jpsi_gaussianA_1, Jpsi_gaussianB_1, Jpsi_gaussianC_1), RooArgList(Jpsi_ratio, Jpsi_ratio2));
    RooAddPdf JpsiMassSig2("JpsiMassSig2", "JpsiMassSig2", RooArgList(Jpsi_gaussianA_2, Jpsi_gaussianB_2, Jpsi_gaussianC_2), RooArgList(Jpsi_ratio, Jpsi_ratio2));
''',
        )
        source = add_freeze(
            source,
            "    Jpsi_ratio.setConstant(kTRUE);\n",
            "    Jpsi_devia3.setConstant(kTRUE);\n    Jpsi_ratio2.setConstant(kTRUE);\n",
        )
    elif variation == "f2_linear":
        source = replace_once(
            source,
            '''    RooChebychev JpsiMassComb1("JpsiMassComb1", "JpsiMassComb1", Jpsi_mass1, RooArgList());
    RooChebychev JpsiMassComb2("JpsiMassComb2", "JpsiMassComb2", Jpsi_mass2, RooArgList());
''',
            '''    RooRealVar Jpsi_a("Jpsi_a", "Jpsi_a", 0, -1, 1);
    RooChebychev JpsiMassComb1("JpsiMassComb1", "JpsiMassComb1", Jpsi_mass1, RooArgList(Jpsi_a));
    RooChebychev JpsiMassComb2("JpsiMassComb2", "JpsiMassComb2", Jpsi_mass2, RooArgList(Jpsi_a));
''',
        )
        source = add_freeze(
            source,
            "    Jpsi_ratio.setConstant(kTRUE);\n",
            "    Jpsi_a.setConstant(kTRUE);\n",
        )
    elif variation in {"f3_single_gaussian", "f3_triple_gaussian"}:
        if variation == "f3_single_gaussian":
            prompt_block = '''    RooRealVar Jpsi_mu1("Jpsi_mu1", "Jpsi_mu1", 0, -0.005, 0.005);
    RooRealVar Jpsi_sigma1("Jpsi_sigma1", "Jpsi_sigma1", 0.001, 0, 0.003);
    RooRealVar Jpsi_sigma2("Jpsi_sigma2", "Jpsi_sigma2", 0.004, 0.002, 0.008);
    RooRealVar Jpsi_prop1("Jpsi_prop1", "Jpsi_prop1", 0.5, 0, 1);
    RooGaussian JpsiCtauSig1("JpsiCtauSig1", "JpsiCtauSig1", Jpsi_ctau1, Jpsi_mu1, Jpsi_sigma1);
    RooGaussian JpsiCtauSig2("JpsiCtauSig2", "JpsiCtauSig2", Jpsi_ctau2, Jpsi_mu1, Jpsi_sigma1);
'''
        else:
            prompt_block = '''    RooRealVar Jpsi_mu1("Jpsi_mu1", "Jpsi_mu1", 0, -0.005, 0.005);
    RooRealVar Jpsi_sigma1("Jpsi_sigma1", "Jpsi_sigma1", 0.001, 0, 0.003);
    RooRealVar Jpsi_sigma2("Jpsi_sigma2", "Jpsi_sigma2", 0.004, 0.002, 0.008);
    RooRealVar Jpsi_sigma4("Jpsi_sigma4", "Jpsi_sigma4", 0.006, 0.0001, 0.01);
    RooRealVar Jpsi_prop1("Jpsi_prop1", "Jpsi_prop1", 0.5, 0, 1);
    RooRealVar Jpsi_prop3("Jpsi_prop3", "Jpsi_prop3", 0.5, 0, 1);
    RooGaussian Jpsi_gauss1_1("Jpsi_gauss1_1", "Jpsi_gauss1_1", Jpsi_ctau1, Jpsi_mu1, Jpsi_sigma1);
    RooGaussian Jpsi_gauss2_1("Jpsi_gauss2_1", "Jpsi_gauss2_1", Jpsi_ctau1, Jpsi_mu1, Jpsi_sigma2);
    RooGaussian Jpsi_gauss3_1("Jpsi_gauss3_1", "Jpsi_gauss3_1", Jpsi_ctau1, Jpsi_mu1, Jpsi_sigma4);
    RooGaussian Jpsi_gauss1_2("Jpsi_gauss1_2", "Jpsi_gauss1_2", Jpsi_ctau2, Jpsi_mu1, Jpsi_sigma1);
    RooGaussian Jpsi_gauss2_2("Jpsi_gauss2_2", "Jpsi_gauss2_2", Jpsi_ctau2, Jpsi_mu1, Jpsi_sigma2);
    RooGaussian Jpsi_gauss3_2("Jpsi_gauss3_2", "Jpsi_gauss3_2", Jpsi_ctau2, Jpsi_mu1, Jpsi_sigma4);
    RooAddPdf JpsiCtauSig1("JpsiCtauSig1", "JpsiCtauSig1", RooArgList(Jpsi_gauss1_1, Jpsi_gauss2_1, Jpsi_gauss3_1), RooArgList(Jpsi_prop1, Jpsi_prop3));
    RooAddPdf JpsiCtauSig2("JpsiCtauSig2", "JpsiCtauSig2", RooArgList(Jpsi_gauss1_2, Jpsi_gauss2_2, Jpsi_gauss3_2), RooArgList(Jpsi_prop1, Jpsi_prop3));
'''
        source = replace_block(
            source,
            '    RooRealVar Jpsi_mu1("Jpsi_mu1"',
            '    RooAddPdf JpsiCtauSig2("JpsiCtauSig2", "JpsiCtauSig2", RooArgList(Jpsi_gauss1_2, Jpsi_gauss2_2), Jpsi_prop1);\n',
            prompt_block,
        )
        if variation == "f3_triple_gaussian":
            source = add_freeze(
                source,
                "    Jpsi_prop1.setConstant(kTRUE);\n",
                "    Jpsi_sigma4.setConstant(kTRUE);\n    Jpsi_prop3.setConstant(kTRUE);\n",
            )
    elif variation in {"f4_exp_double_gaussian", "f4_double_exp_gaussian"}:
        if variation == "f4_exp_double_gaussian":
            nonprompt_block = '''    RooRealVar Jpsi_sigma3("Jpsi_sigma3", "Jpsi_sigma3", 0.001, 0, 0.01);
    RooRealVar Jpsi_sigma5("Jpsi_sigma5", "Jpsi_sigma5", 0.005, 0.001, 0.02);
    RooRealVar Jpsi_coef1("Jpsi_coef1", "Jpsi_coef1", 0.06, 0.01, 0.1);
    RooRealVar Jpsi_prop4("Jpsi_prop4", "Jpsi_prop4", 0.5, 0, 1);
    RooGExpModel Jpsi_expgs1_1("Jpsi_expgs1_1", "Jpsi_expgs1_1", Jpsi_ctau1, Jpsi_sigma3, Jpsi_coef1, false, RooGExpModel::Type::Flipped);
    RooGExpModel Jpsi_expgs2_1("Jpsi_expgs2_1", "Jpsi_expgs2_1", Jpsi_ctau1, Jpsi_sigma5, Jpsi_coef1, false, RooGExpModel::Type::Flipped);
    RooGExpModel Jpsi_expgs1_2("Jpsi_expgs1_2", "Jpsi_expgs1_2", Jpsi_ctau2, Jpsi_sigma3, Jpsi_coef1, false, RooGExpModel::Type::Flipped);
    RooGExpModel Jpsi_expgs2_2("Jpsi_expgs2_2", "Jpsi_expgs2_2", Jpsi_ctau2, Jpsi_sigma5, Jpsi_coef1, false, RooGExpModel::Type::Flipped);
    RooAddPdf JpsiCtauBkg1("JpsiCtauBkg1", "JpsiCtauBkg1", RooArgList(Jpsi_expgs1_1, Jpsi_expgs2_1), Jpsi_prop4);
    RooAddPdf JpsiCtauBkg2("JpsiCtauBkg2", "JpsiCtauBkg2", RooArgList(Jpsi_expgs1_2, Jpsi_expgs2_2), Jpsi_prop4);
'''
            extra_freeze = "    Jpsi_sigma5.setConstant(kTRUE);\n    Jpsi_prop4.setConstant(kTRUE);\n"
        else:
            nonprompt_block = '''    RooRealVar Jpsi_sigma3("Jpsi_sigma3", "Jpsi_sigma3", 0.001, 0, 0.01);
    RooRealVar Jpsi_coef1("Jpsi_coef1", "Jpsi_coef1", 0.06, 0.01, 0.1);
    RooRealVar Jpsi_coef2("Jpsi_coef2", "Jpsi_coef2", 0.03, 0.02, 0.5);
    RooRealVar Jpsi_prop4("Jpsi_prop4", "Jpsi_prop4", 0.5, 0, 1);
    RooGExpModel Jpsi_expgs1_1("Jpsi_expgs1_1", "Jpsi_expgs1_1", Jpsi_ctau1, Jpsi_sigma3, Jpsi_coef1, false, RooGExpModel::Type::Flipped);
    RooGExpModel Jpsi_expgs2_1("Jpsi_expgs2_1", "Jpsi_expgs2_1", Jpsi_ctau1, Jpsi_sigma3, Jpsi_coef2, false, RooGExpModel::Type::Flipped);
    RooGExpModel Jpsi_expgs1_2("Jpsi_expgs1_2", "Jpsi_expgs1_2", Jpsi_ctau2, Jpsi_sigma3, Jpsi_coef1, false, RooGExpModel::Type::Flipped);
    RooGExpModel Jpsi_expgs2_2("Jpsi_expgs2_2", "Jpsi_expgs2_2", Jpsi_ctau2, Jpsi_sigma3, Jpsi_coef2, false, RooGExpModel::Type::Flipped);
    RooAddPdf JpsiCtauBkg1("JpsiCtauBkg1", "JpsiCtauBkg1", RooArgList(Jpsi_expgs1_1, Jpsi_expgs2_1), Jpsi_prop4);
    RooAddPdf JpsiCtauBkg2("JpsiCtauBkg2", "JpsiCtauBkg2", RooArgList(Jpsi_expgs1_2, Jpsi_expgs2_2), Jpsi_prop4);
'''
            extra_freeze = "    Jpsi_coef2.setConstant(kTRUE);\n    Jpsi_prop4.setConstant(kTRUE);\n"
        source = replace_block(
            source,
            '    RooRealVar Jpsi_sigma3("Jpsi_sigma3"',
            '    RooGExpModel JpsiCtauBkg2("JpsiCtauBkg2", "JpsiCtauBkg2", Jpsi_ctau2, Jpsi_sigma3, Jpsi_coef1, false, RooGExpModel::Type::Flipped);\n',
            nonprompt_block,
        )
        source = add_freeze(source, "    Jpsi_coef1.setConstant(kTRUE);\n", extra_freeze)
    else:
        raise RuntimeError(f"unsupported variation: {variation}")

    return source


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("variation", choices=sorted(VARIATIONS))
    parser.add_argument("base", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    source_bytes = args.base.read_bytes()
    digest = hashlib.sha256(source_bytes).hexdigest()
    if digest != BASE_SHA256:
        raise RuntimeError(f"frozen Fit_4D_tot.cpp hash mismatch: {digest}")
    generated = generate(source_bytes.decode(), args.variation)
    if args.output.exists():
        raise RuntimeError(f"refusing to overwrite {args.output}")
    args.output.write_text(generated)


if __name__ == "__main__":
    main()
