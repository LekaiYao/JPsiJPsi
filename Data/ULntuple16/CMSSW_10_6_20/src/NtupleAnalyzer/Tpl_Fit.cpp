#include "TCanvas.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "RooExtendPdf.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooFitResult.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"

using namespace std;
using namespace RooFit;

constexpr long double operator"" _PI(long double f) {
    return 3.14159265359 * f;
}
double combError(double sys1, double sys2, double sys3, double sys4) {
    return sqrt(sys1 * sys1 + sys2 * sys2 + sys3 * sys3 + sys4 * sys4);
}
void plot_temp(string varName, RooRealVar &var, RooDataHist *dh, RooDataHist *dhSPS, RooDataHist *dhDPS, double frac, double n_sps_dps) {
    TCanvas *canvas = new TCanvas("canvas", "canvas", 1500, 1300);
    RooPlot* frame = var.frame();
    RooHistFunc func_SPS("func_SPS", "func_SPS", var, *dhSPS), func_DPS("func_DPS", "func_DPS", var, *dhDPS);
    RooRealVar coef_SPS("coef_SPS", "coef_SPS", frac * n_sps_dps);
    RooRealVar coef_DPS("coef_DPS", "coef_DPS", (1 - frac) * n_sps_dps);
    RooRealSumFunc func("func", "func", {func_SPS, func_DPS}, {coef_SPS, coef_DPS});
    dh->plotOn(frame, Name("Data"));
    func.plotOn(frame, LineColor(kBlack), Name("All")); 
    func_SPS.plotOn(frame, Normalization(coef_SPS.getVal()), LineColor(kBlue), LineStyle(kDashed), Name("SPS"));
    func_DPS.plotOn(frame, Normalization(coef_DPS.getVal()), LineColor(kRed), LineStyle(kDotted), Name("DPS"));
    TLegend *legend = new TLegend(.65, .60, .85, .85);
    legend->AddEntry(frame->findObject("Data"), "RunII 2018", "L");
    legend->AddEntry(frame->findObject("All"), "Total p.d.f.", "L");
    legend->AddEntry(frame->findObject("SPS"), "SPS p.d.f.", "L");
    legend->AddEntry(frame->findObject("DPS"), "DPS p.d.f.", "L");
    frame->SetXTitle(varName.c_str());
    frame->SetYTitle(("d#sigma/d(" + varName + ")").c_str());
    frame->Draw();
    legend->DrawClone();
    canvas->SaveAs(("fig/temp/Template_" + varName + ".png").c_str());
}
void Tpl_Fit() {
    // Names of kinematic variables
    const int varNum = 5;
    const string varName[] = {"delta_y", "delta_phi", "evt_mass", "evt_y", "evt_pt"};
    // Binning configurations
    const int binNum[] = {6, 8, 7, 5, 9};
    const vector<vector<double>> bins = {
        {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 4.0},
        {0.0_PI, 0.125_PI, 0.25_PI, 0.375_PI, 0.5_PI, 0.625_PI, 0.75_PI, 0.875_PI, 1.0_PI},
        {7.5, 17.5, 27.5, 37.5, 47.5, 57.5, 67.5, 107.5},
        {0.0, 0.4, 0.8, 1.2, 1.6, 2.0},
        {0, 5, 10, 15, 20, 25, 30, 35, 40, 80}
    };
    // Differential cross sections(event yield)
    vector<vector<double>> xSec = {
        {2693.07, 897.703, 540.821, 428.172, 285.705, 213.566},
        {1374.97, 794.097, 232.649, 357.403, 303.331, 347.662, 651.977, 1041.69},
        {916.264, 907.639, 1143.22, 836.383, 436.256, 175.983, 179.425},
        {1300.21, 1365.8, 1318.3, 794.153, 263.802},
        {678.57, 696.301, 574.91, 334.685, 1169.18, 835.068, 286.57, 237.961, 284.646}
    };
    // Statistical uncertainties(error of event yield)
    vector<vector<double>> sta = {
        {55.546, 33.2241, 26.6067, 23.8414, 19.3493, 16.9672},
        {38.8299, 30.0436, 17.0393, 20.5816, 19.1955, 21.0125, 28.8214, 36.9809},
        {33.2635, 32.3605, 37.5913, 31.8681, 22.7294, 15.3576, 14.8315},
        {40.2643, 40.9217, 39.7716, 30.4744, 17.9819},
        {28.7697, 30.4292, 26.7814, 20.8368, 36.4629, 31.4697, 18.0875, 16.4044, 18.0145}
    };
    for(int i = 0; i < varNum; i++) {
        for(int j = 0; j < binNum[i]; j++) {
            double x = 1e-3 / (36.3 * 0.05961 * 0.05961) / (bins[i][j + 1] - bins[i][j]);
            xSec[i][j] *= x;
            sta[i][j] *= x;
        }
    }
    // Systematic uncertainties(in percentage)
    const double sys1 = 0.015, sys2 = 0.026;
    const vector<vector<double>> sys3 = {
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0}
    }, sys4 = {
        {0.0125, 0.0359, 0.0492, 0.0239, 0.0973, 0.1003},
        {0.0162, 0.04378, 0.0778, 0.1242, 0.077, 0.0483, 0.0493, 0.0212},
        {0.0594, 0.039, 0.0293, 0.0681, 0.0402, 0.1205, 0.101},
        {0.0499, 0.0138, 0.0349, 0.0473, 0.0965},
        {0.0983, 0.0523, 0.0653, 0.0399, 0.0128, 0.0947, 0.091, 0.077, 0.1225}
    };
    // Loop on all kinematic variables
    for(int i = 0; i < varNum; i++) {
        // Draw distribution only, with statistical and systematic errors separated
        TH1D *h1 = new TH1D(("d_" + varName[i]).c_str(), ("d_" + varName[i]).c_str(), binNum[i], bins[i].data());
        TH1D *h2 = new TH1D(("e_" + varName[i]).c_str(), ("e_" + varName[i]).c_str(), binNum[i], bins[i].data());
        for(int j = 0; j < binNum[i]; j++) {
            h1->Fill(bins[i][j] + 0.01, xSec[i][j]);
            h1->SetBinError(j + 1, sta[i][j]);
            h2->Fill(bins[i][j] + 0.01, xSec[i][j]);
            h2->SetBinError(j + 1, sta[i][j] + xSec[i][j] * combError(sys1, sys2, sys3[i][j], sys4[i][j]));
        }
        TCanvas *canvas2 = new TCanvas("canvas2", "canvas2", 1000, 1000);
        h1->SetMarkerSize(0.);
        h1->SetFillStyle(3003);
        h1->SetMarkerColor(kRed);
        h1->SetFillColor(kRed);
        h1->GetYaxis()->SetRangeUser(0, 1.2 * h1->GetMaximum());
        h1->GetXaxis()->SetTitle(varName[i].c_str());
        h1->GetYaxis()->SetTitle(("d#sigma/d(" + varName[i] + ")").c_str());
        h1->SetStats(0);
        h1->Draw("e2");
        h2->SetStats(0);
        h2->Draw("same");
        canvas2->SaveAs(("fig/temp/Xsec_" + varName[i] + ".png").c_str());
    }
    return;
}