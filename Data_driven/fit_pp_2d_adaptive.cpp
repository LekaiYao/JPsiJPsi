#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH2D.h"
#include "TH2Poly.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TTreeFormula.h"

#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

using namespace RooFit;

namespace {

constexpr double kPi = 3.14159265358979323846;
const std::vector<double> kDyEdges = {0.0, 0.6, 1.2, 1.8,
                                      2.4, 4.0};
const std::vector<double> kDphiEdges = {0.0, kPi / 2.0,
                                        3.0 * kPi / 4.0, kPi};

struct FitCell {
    int iy = -1;
    int iphi = -1;
    double dyLow = 0;
    double dyHigh = 0;
    double dphiLow = 0;
    double dphiHigh = 0;
    Long64_t treeEntries = 0;
    int fitEntries = 0;
    double sumWeights = 0;
    double sumWeights2 = 0;
    double effectiveEntries = 0;
    double ppYield = std::numeric_limits<double>::quiet_NaN();
    double ppError = std::numeric_limits<double>::quiet_NaN();
    double pNpYield = std::numeric_limits<double>::quiet_NaN();
    double npNpYield = std::numeric_limits<double>::quiet_NaN();
    double sigCombYield = std::numeric_limits<double>::quiet_NaN();
    double combCombYield = std::numeric_limits<double>::quiet_NaN();
    int status = 999;
    int covQual = -1;
    double edm = std::numeric_limits<double>::quiet_NaN();
    double minNll = std::numeric_limits<double>::quiet_NaN();
    int attempts = 0;
    bool accepted = false;
};

std::string makeCut(double dyLow, double dyHigh, double dphiLow,
                    double dphiHigh, bool lastDy, bool lastDphi) {
    std::ostringstream out;
    out << std::setprecision(17)
        << "delta_y >= " << dyLow << " && delta_y "
        << (lastDy ? "<= " : "< ") << dyHigh << " && delta_phi >= "
        << dphiLow << " && delta_phi " << (lastDphi ? "<= " : "< ")
        << dphiHigh;
    return out.str();
}

void configureYields(RooWorkspace &workspace, double scale, int attempt) {
    struct YieldConfig {
        const char *name;
        double minimum;
    };
    const YieldConfig yields[] = {{"n_P_P", 1.0},      {"n_P_NP", 0.0},
                                  {"n_NP_NP", 0.0},    {"n_Sig_Comb", 0.0},
                                  {"n_Comb_Comb", 0.0}};
    const double deterministicShift[] = {1.0, 0.8, 1.2, 0.6, 1.5, 0.4};
    const double shift = deterministicShift[attempt % 6];
    for (const YieldConfig &cfg : yields) {
        RooRealVar *var = workspace.var(cfg.name);
        if (!var)
            continue;
        const double global = var->getVal();
        var->setMin(cfg.minimum);
        var->setMax(std::max(1000.0, global * std::max(2.0, 10.0 * scale)));
        var->setVal(std::max(cfg.minimum + 1e-3, global * scale * shift));
        var->setConstant(false);
    }
}

FitCell fitCell(TTree &tree, const FitCell &definition,
                const char *modelPath) {
    FitCell out = definition;
    const bool lastDy = std::abs(out.dyHigh - kDyEdges.back()) < 1e-12;
    const bool lastDphi = std::abs(out.dphiHigh - kPi) < 1e-12;
    const std::string cut =
        makeCut(out.dyLow, out.dyHigh, out.dphiLow, out.dphiHigh, lastDy,
                lastDphi);
    out.treeEntries = tree.GetEntries(cut.c_str());

    RooRealVar mass1("Jpsi_mass1", "Jpsi_mass1", 2.95, 3.25);
    RooRealVar mass2("Jpsi_mass2", "Jpsi_mass2", 2.95, 3.25);
    RooRealVar ctau1("Jpsi_ctau1", "Jpsi_ctau1", -0.03, 0.16);
    RooRealVar ctau2("Jpsi_ctau2", "Jpsi_ctau2", -0.03, 0.16);
    RooRealVar weight("evt_weight", "evt_weight", 0, 1000);
    RooRealVar dy("delta_y", "delta_y", 0, 4.0);
    RooRealVar dphi("delta_phi", "delta_phi", 0, kPi);
    RooArgSet inputVars(mass1, mass2, ctau1, ctau2, weight, dy, dphi);
    RooDataSet data("data", "data", inputVars, Import(tree),
                    Cut(cut.c_str()), WeightVar("evt_weight"));
    std::unique_ptr<RooAbsData> fitData(
        data.reduce(RooArgSet(mass1, mass2, ctau1, ctau2)));
    out.fitEntries = fitData ? fitData->numEntries() : 0;
    out.sumWeights = fitData ? fitData->sumEntries() : 0;

    double evtWeight = 0;
    tree.SetBranchAddress("evt_weight", &evtWeight);
    std::unique_ptr<TTreeFormula> formula(
        new TTreeFormula("cell_cut", cut.c_str(), &tree));
    for (Long64_t i = 0; i < tree.GetEntries(); ++i) {
        tree.GetEntry(i);
        if (formula->EvalInstance())
            out.sumWeights2 += evtWeight * evtWeight;
    }
    tree.ResetBranchAddresses();
    if (out.sumWeights2 > 0)
        out.effectiveEntries =
            out.sumWeights * out.sumWeights / out.sumWeights2;
    if (!fitData || out.fitEntries < 10 || out.sumWeights <= 0)
        return out;

    TFile modelFile(modelPath, "READ");
    RooWorkspace *workspace =
        dynamic_cast<RooWorkspace *>(modelFile.Get("wsp"));
    if (!workspace)
        return out;
    RooAbsPdf *pdf = workspace->pdf("pdf_all");
    if (!pdf)
        return out;

    const double globalExpected =
        workspace->var("n_P_P")->getVal() +
        2.0 * workspace->var("n_P_NP")->getVal() +
        workspace->var("n_NP_NP")->getVal() +
        2.0 * workspace->var("n_Sig_Comb")->getVal() +
        workspace->var("n_Comb_Comb")->getVal();
    const double scale =
        globalExpected > 0 ? out.sumWeights / globalExpected : 0.05;

    double bestScore = std::numeric_limits<double>::infinity();
    for (int attempt = 0; attempt < 6; ++attempt) {
        configureYields(*workspace, scale, attempt);
        std::unique_ptr<RooFitResult> result(pdf->fitTo(
            *fitData, Extended(true), Save(true), AsymptoticError(true), Offset(true),
            PrintLevel(-1), Strategy(1), Optimize(false),
            Minimizer("Minuit2", "migrad")));
        ++out.attempts;
        if (!result)
            continue;
        const double score =
            1000.0 * std::abs(result->status()) +
            100.0 * std::max(0, 3 - result->covQual()) +
            std::min(100.0, std::abs(result->edm()));
        if (score < bestScore) {
            bestScore = score;
            out.status = result->status();
            out.covQual = result->covQual();
            out.edm = result->edm();
            out.minNll = result->minNll();
            RooRealVar *pp = workspace->var("n_P_P");
            out.ppYield = pp ? pp->getVal()
                             : std::numeric_limits<double>::quiet_NaN();
            out.ppError = pp ? pp->getError()
                             : std::numeric_limits<double>::quiet_NaN();
            RooRealVar *pNp = workspace->var("n_P_NP");
            RooRealVar *npNp = workspace->var("n_NP_NP");
            RooRealVar *sigComb = workspace->var("n_Sig_Comb");
            RooRealVar *combComb = workspace->var("n_Comb_Comb");
            out.pNpYield = pNp ? pNp->getVal() : NAN;
            out.npNpYield = npNp ? npNp->getVal() : NAN;
            out.sigCombYield = sigComb ? sigComb->getVal() : NAN;
            out.combCombYield = combComb ? combComb->getVal() : NAN;
        }
        if (result->status() == 0 && result->covQual() >= 2 &&
            result->edm() < 0.01) {
            out.accepted = true;
            break;
        }
    }
    return out;
}

std::vector<FitCell> makeCells(bool singleControlCell,
                               bool mergeHighDphiRows,
                               bool mergeHighDyColumns,
                               double controlDyMin,
                               double controlPhiMax) {
    std::vector<FitCell> cells;
    if (!singleControlCell) {
        for (int iy = 0; iy + 1 < static_cast<int>(kDyEdges.size()); ++iy) {
            for (int iphi = 0;
                 iphi + 1 < static_cast<int>(kDphiEdges.size()); ++iphi) {
                if (iy == 3 && iphi == 2)
                    continue;
                FitCell cell;
                cell.iy = iy;
                cell.iphi = iphi;
                cell.dyLow = kDyEdges[iy];
                cell.dyHigh = kDyEdges[iy + 1];
                cell.dphiLow = kDphiEdges[iphi];
                cell.dphiHigh =
                    (iy == 3 && iphi == 1) ? kPi : kDphiEdges[iphi + 1];
                cells.push_back(cell);
            }
        }
        return cells;
    }

    std::vector<double> dphiEdges = {0.0, controlPhiMax,
                                     3.0 * kPi / 4.0, kPi};
    int controlDyIndex = -1;
    for (int iy = 0; iy + 1 < static_cast<int>(kDyEdges.size()); ++iy) {
        if (std::abs(kDyEdges[iy] - controlDyMin) < 1e-12) {
            controlDyIndex = iy;
            break;
        }
    }
    if (controlDyIndex < 0 || !(controlPhiMax > 0.0) ||
        !(controlPhiMax < 3.0 * kPi / 4.0)) {
        std::cerr << "Unsupported single-control-cell boundaries: dy="
                  << controlDyMin << " phi=" << controlPhiMax << std::endl;
        return cells;
    }

    // Keep a 3-column grid below the control-region dy boundary. Moving the
    // CR boundary therefore repartitions the full plane without gaps.
    for (int iy = 0; iy < controlDyIndex; ++iy) {
        for (int iphi = 0; iphi < 3; ++iphi) {
            FitCell cell;
            cell.iy = iy;
            cell.iphi = iphi;
            cell.dyLow = kDyEdges[iy];
            cell.dyHigh = kDyEdges[iy + 1];
            cell.dphiLow = dphiEdges[iphi];
            cell.dphiHigh = dphiEdges[iphi + 1];
            cells.push_back(cell);
        }
    }

    // The nominal assumption is that this entire region is DPS, so extract one
    // prompt-prompt yield for the full control region.
    FitCell control;
    control.iy = 3;
    control.iphi = 0;
    control.dyLow = controlDyMin;
    control.dyHigh = 4.0;
    control.dphiLow = 0.0;
    control.dphiHigh = controlPhiMax;
    cells.push_back(control);

    // Alternative 12-cell validation: retain the two high-dphi slices but
    // merge their dy extent across the full 1.8--4.0 interval.
    if (mergeHighDyColumns) {
        for (int iphi = 1; iphi < 3; ++iphi) {
            FitCell cell;
            cell.iy = 3;
            cell.iphi = iphi;
            cell.dyLow = controlDyMin;
            cell.dyHigh = 4.0;
            cell.dphiLow = dphiEdges[iphi];
            cell.dphiHigh = dphiEdges[iphi + 1];
            cells.push_back(cell);
        }
        return cells;
    }

    // Retain dy=2.4 in the high-dphi region.  The 14-cell validation keeps
    // both dphi slices; the lower-statistics 12-cell validation merges them.
    for (int iy = controlDyIndex;
         iy + 1 < static_cast<int>(kDyEdges.size()); ++iy) {
        if (mergeHighDphiRows) {
            FitCell cell;
            cell.iy = iy;
            cell.iphi = 1;
            cell.dyLow = kDyEdges[iy];
            cell.dyHigh = kDyEdges[iy + 1];
            cell.dphiLow = controlPhiMax;
            cell.dphiHigh = kPi;
            cells.push_back(cell);
            continue;
        }
        for (int iphi = 1; iphi < 3; ++iphi) {
            FitCell cell;
            cell.iy = iy;
            cell.iphi = iphi;
            cell.dyLow = kDyEdges[iy];
            cell.dyHigh = kDyEdges[iy + 1];
            cell.dphiLow = dphiEdges[iphi];
            cell.dphiHigh = dphiEdges[iphi + 1];
            cells.push_back(cell);
        }
    }
    return cells;
}

void drawMap(TH2D &hist, const std::string &path, const char *label) {
    TCanvas canvas("map_canvas", "", 900, 700);
    canvas.SetRightMargin(0.14);
    hist.SetStats(false);
    hist.SetTitle("");
    hist.GetXaxis()->SetTitle("|#Delta y|");
    hist.GetYaxis()->SetTitle("|#Delta#phi|");
    hist.Draw("COLZ TEXT E");
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.12, 0.94, label);
    canvas.SaveAs(path.c_str());
}

void drawMap(TH2Poly &hist, const std::string &path, const char *label) {
    TCanvas canvas("poly_map_canvas", "", 900, 700);
    canvas.SetRightMargin(0.14);
    hist.SetStats(false);
    hist.SetTitle("");
    hist.GetXaxis()->SetTitle("|#Delta y|");
    hist.GetYaxis()->SetTitle("|#Delta#phi|");
    hist.Draw("COLZ TEXT E");
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.12, 0.94, label);
    canvas.SaveAs(path.c_str());
}

} // namespace

void fit_pp_2d_adaptive(
    const char *outputDir =
        "Data_driven/results/pp_data_2d_adaptive14",
    const char *dataPath =
        "Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/WeightData.root",
    bool singleControlCell = false, bool mergeHighDphiRows = false,
    bool mergeHighDyColumns = false, double controlDyMin = 1.8,
    double controlPhiMax = kPi / 2.0,
    const char *modelPath = "Data_driven/inputs/Model_4D_tot.root") {
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    gSystem->mkdir(outputDir, true);
    TFile dataFile(dataPath, "READ");
    TTree *tree = dynamic_cast<TTree *>(dataFile.Get("data"));
    if (!tree) {
        std::cerr << "Cannot open input data tree." << std::endl;
        return;
    }

    TH2D yieldMap("h_data_pp_dy_dphi",
                  "corrected prompt-prompt yield;|#Delta y|;|#Delta#phi|",
                  kDyEdges.size() - 1, kDyEdges.data(),
                  kDphiEdges.size() - 1, kDphiEdges.data());
    TH2D rawMap("h_data_raw_entries_dy_dphi",
                "raw selected entries;|#Delta y|;|#Delta#phi|",
                kDyEdges.size() - 1, kDyEdges.data(),
                kDphiEdges.size() - 1, kDphiEdges.data());
    TH2D statusMap("h_fit_status_dy_dphi",
                   "fit status;|#Delta y|;|#Delta#phi|",
                   kDyEdges.size() - 1, kDyEdges.data(),
                   kDphiEdges.size() - 1, kDphiEdges.data());
    TH2D covMap("h_fit_covqual_dy_dphi",
                "fit covariance quality;|#Delta y|;|#Delta#phi|",
                kDyEdges.size() - 1, kDyEdges.data(),
                kDphiEdges.size() - 1, kDphiEdges.data());

    const std::vector<FitCell> definitions =
        makeCells(singleControlCell, mergeHighDphiRows,
                  mergeHighDyColumns, controlDyMin, controlPhiMax);
    std::vector<FitCell> cells;
    for (const FitCell &definition : definitions) {
            std::cout << "Fit cell dy=[" << definition.dyLow << ","
                      << definition.dyHigh << "] dphi=["
                      << definition.dphiLow << ","
                      << definition.dphiHigh << "]" << std::endl;
            FitCell cell = fitCell(*tree, definition, modelPath);
            cells.push_back(cell);
            const int bx = cell.iy + 1, by = cell.iphi + 1;
            if (std::isfinite(cell.ppYield)) {
                yieldMap.SetBinContent(bx, by, cell.ppYield);
                yieldMap.SetBinError(bx, by, cell.ppError);
            }
            rawMap.SetBinContent(bx, by, cell.treeEntries);
            statusMap.SetBinContent(bx, by, cell.status);
            covMap.SetBinContent(bx, by, cell.covQual);
            std::cout << "  raw=" << cell.treeEntries
                      << " fit_entries=" << cell.fitEntries
                      << " sumw=" << cell.sumWeights
                      << " neff=" << cell.effectiveEntries
                      << " PP=" << cell.ppYield << " +/- " << cell.ppError
                      << " status=" << cell.status
                      << " covQual=" << cell.covQual
                      << " edm=" << cell.edm
                      << " attempts=" << cell.attempts << std::endl;
    }

    const std::string base(outputDir);
    TFile output((base + "/pp_data_2d.root").c_str(), "RECREATE");
    if (!singleControlCell) {
        yieldMap.Write();
        rawMap.Write();
        statusMap.Write();
        covMap.Write();
    } else {
        TH2Poly yieldPoly, rawPoly, statusPoly, covPoly;
        yieldPoly.SetName("h_data_pp_dy_dphi");
        rawPoly.SetName("h_data_raw_entries_dy_dphi");
        statusPoly.SetName("h_fit_status_dy_dphi");
        covPoly.SetName("h_fit_covqual_dy_dphi");
        for (const FitCell &cell : cells) {
            yieldPoly.AddBin(cell.dyLow, cell.dphiLow, cell.dyHigh,
                             cell.dphiHigh);
            rawPoly.AddBin(cell.dyLow, cell.dphiLow, cell.dyHigh,
                           cell.dphiHigh);
            statusPoly.AddBin(cell.dyLow, cell.dphiLow, cell.dyHigh,
                              cell.dphiHigh);
            covPoly.AddBin(cell.dyLow, cell.dphiLow, cell.dyHigh,
                           cell.dphiHigh);
        }
        for (size_t i = 0; i < cells.size(); ++i) {
            const int bin = static_cast<int>(i) + 1;
            const FitCell &cell = cells[i];
            if (std::isfinite(cell.ppYield)) {
                yieldPoly.SetBinContent(bin, cell.ppYield);
                yieldPoly.SetBinError(bin, cell.ppError);
            }
            rawPoly.SetBinContent(bin, cell.treeEntries);
            statusPoly.SetBinContent(bin, cell.status);
            covPoly.SetBinContent(bin, cell.covQual);
        }
        yieldPoly.Write();
        rawPoly.Write();
        statusPoly.Write();
        covPoly.Write();
        drawMap(yieldPoly, base + "/pp_yield_map.pdf",
                "Prompt-prompt corrected yield");
        drawMap(rawPoly, base + "/raw_entries_map.pdf",
                "Raw selected entries");
        drawMap(statusPoly, base + "/fit_status_map.pdf", "Fit status");
        drawMap(covPoly, base + "/fit_covqual_map.pdf",
                "Fit covariance quality");
    }
    output.Close();

    std::ofstream csv(base + "/fit_results.csv");
    csv << "iy,iphi,dy_low,dy_high,dphi_low,dphi_high,tree_entries,"
           "fit_entries,sum_weights,sum_weights2,effective_entries,pp_yield,"
           "pp_error,status,covQual,edm,minNll,attempts,accepted,p_np_yield,"
           "np_np_yield,sig_comb_yield,comb_comb_yield\n";
    csv << std::setprecision(12);
    for (const FitCell &c : cells)
        csv << c.iy << "," << c.iphi << "," << c.dyLow << ","
            << c.dyHigh << "," << c.dphiLow << "," << c.dphiHigh << ","
            << c.treeEntries << "," << c.fitEntries << "," << c.sumWeights
            << "," << c.sumWeights2 << "," << c.effectiveEntries << ","
            << c.ppYield << "," << c.ppError << "," << c.status << ","
            << c.covQual << "," << c.edm << "," << c.minNll << ","
            << c.attempts << "," << (c.accepted ? 1 : 0) << ","
            << c.pNpYield << "," << c.npNpYield << "," << c.sigCombYield
            << "," << c.combCombYield << "\n";

    int accepted = 0;
    for (const FitCell &c : cells)
        accepted += c.accepted;
    std::ofstream summary(base + "/summary.txt");
    summary << "input=" << dataPath << "\n"
            << "model=" << modelPath << "\n"
            << "data_tree_entries=" << tree->GetEntries() << "\n"
            << "cells=" << cells.size() << "\n"
            << "accepted_fits=" << accepted << "\n"
            << "failed_fits=" << cells.size() - accepted << "\n"
            << "single_control_cell=" << (singleControlCell ? 1 : 0) << "\n"
            << "merge_high_dphi_rows=" << (mergeHighDphiRows ? 1 : 0) << "\n"
            << "merge_high_dy_columns=" << (mergeHighDyColumns ? 1 : 0) << "\n"
            << "control_dy_min=" << controlDyMin << "\n"
            << "control_phi_max=" << controlPhiMax << "\n"
            << "acceptance_criterion=status==0 && covQual>=2 && edm<0.01\n"
            << "error_convention=ROOT_native_AsymptoticError_true\n"
            << "likelihood_offset=true\n";
    if (!singleControlCell) {
        drawMap(yieldMap, base + "/pp_yield_map.pdf",
                "Prompt-prompt corrected yield");
        drawMap(rawMap, base + "/raw_entries_map.pdf", "Raw selected entries");
        drawMap(statusMap, base + "/fit_status_map.pdf", "Fit status");
        drawMap(covMap, base + "/fit_covqual_map.pdf",
                "Fit covariance quality");
    }

    std::cout << "accepted_fits=" << accepted << "/" << cells.size()
              << "\noutput_dir=" << outputDir << std::endl;
}
