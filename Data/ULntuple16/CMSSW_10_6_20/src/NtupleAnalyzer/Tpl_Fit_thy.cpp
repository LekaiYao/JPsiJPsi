#include "TCanvas.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TLatex.h"
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
double combError(double sys1, double sys2, double sys3, double sys4, double sys5) {
    return sqrt(sys1 * sys1 + sys2 * sys2 + sys3 * sys3 + sys4 * sys4 + sys5 * sys5);
}
void add_latex(TLatex &latex) {
    latex.SetTextSize(0.05);
    latex.SetTextFont(62);
    latex.DrawLatex(0.19, 0.85, "CMS");
    latex.SetTextSize(0.04);
    latex.SetTextFont(52);
    latex.DrawLatex(0.30, 0.85, "Preliminary");
    latex.SetTextSize(0.03);
    latex.SetTextFont(42);
    latex.DrawLatex(0.68, 0.91, "36.31 fb^{-1} (13TeV)");
    return;
}
// Plot the result of template fit
void plot_temp(string varName, RooRealVar &var, RooDataHist *dh, RooDataHist *dhSPS, RooDataHist *dhDPS, RooDataHist *dhThy, double frac, double n_sps_dps, string xTitle, string yTitle, double yMax, double yMin) {
    TCanvas *canvas = new TCanvas("canvas", "canvas", 1500, 1500);
    gPad->SetLeftMargin(0.15);
    RooPlot *frame = var.frame();
    RooHistFunc func_SPS("func_SPS", "func_SPS", var, *dhSPS), func_DPS("func_DPS", "func_DPS", var, *dhDPS);
    RooRealVar coef_SPS("coef_SPS", "coef_SPS", frac * n_sps_dps);
    RooRealVar coef_DPS("coef_DPS", "coef_DPS", (1 - frac) * n_sps_dps);
    RooRealSumFunc func("func", "func", {func_SPS, func_DPS}, {coef_SPS, coef_DPS});
    dhThy->plotOn(frame, Name("Theory"), DrawOption("E2"), FillColor(45), LineWidth(4));
    // func_SPS.plotOn(frame, Normalization(coef_SPS.getVal()), LineColor(kRed), LineStyle(kDashed), LineWidth(2), Name("SPS"));
    dh->plotOn(frame, Name("Data"), LineColor(kBlack), LineWidth(2), MarkerSize(1), MarkerStyle(20));//, DrawOption("E1")
    // func.plotOn(frame, LineColor(kBlue), Name("All"));
    func_SPS.plotOn(frame, Normalization(coef_SPS.getVal()), LineColor(kRed), LineStyle(kDashed), LineWidth(4), Name("SPS"));
    func_SPS.plotOn(frame, Normalization(coef_SPS.getVal() * 1.41), LineColor(kBlue), LineStyle(kDashed), LineWidth(4), Name("SPS+"));
    func_SPS.plotOn(frame, Normalization(coef_SPS.getVal() * 0.82), LineColor(kMagenta), LineStyle(kDashed), LineWidth(4), Name("SPS-"));
    RooHist* pull = frame->pullHist("Data", "All");
    TLegend *legend = new TLegend(.55, .60, .85, .85);
    legend->AddEntry(frame->findObject("Data"), "RunII 2016", "L");
    // legend->AddEntry(frame->findObject("All"), "Total p.d.f.", "L");
    legend->AddEntry(frame->findObject("SPS"), "SPS contribution", "L");
    legend->AddEntry(frame->findObject("SPS+"), "SPS upper limit", "L");
    legend->AddEntry(frame->findObject("SPS-"), "SPS lower limit", "L");
    legend->AddEntry(frame->findObject("Theory"), "SPS prediction", "F");
    frame->SetTitle("");
    frame->SetYTitle(yTitle.c_str());
    frame->SetXTitle(xTitle.c_str());
    frame->SetAxisRange(yMin, yMax, "Y");
    gPad->SetLogy();
    frame->Draw();
    legend->DrawClone();
    TLatex latex;
    latex.SetNDC();
    add_latex(latex);
    canvas->SaveAs(("fig/temp/Theory_" + varName + ".png").c_str());
}

void Tpl_Fit_thy() {
    // Names of kinematic variables
    const int varNum = 2;
    const string varName[] = {"delta_y", "delta_phi"};
    const string varLatx[] = {"|#Delta y(J/#psi_{1},J/#psi_{2})|", "|#Delta#phi(J/#psi_{1},J/#psi_{2})|"};
    const string varUnit[] = {"", ""};
    // Binning configurations
    const int binNum[] = {6, 6};
    const vector<vector<double>> bins = {
        {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 4.0},
        {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0},
    };
    // Differential cross sections(event yield)
    vector<vector<double>> xSec = {
        {2585.85, 848.6, 514.317, 390.036, 264.352, 207.22},
        {1603.09, 646.887, 302.983, 362.62, 563.976, 882.638},
    };
    // Statistical uncertainties(error of event yield)
    vector<vector<double>> sta = {
        {54.3052, 32.3728, 25.6471, 22.7905, 18.318, 16.6171},
        {42.1743, 27.3204, 19.1803, 21.0966, 25.9222, 34.0112},
    };
    for(int i = 0; i < varNum; i++) {
        for(int j = 0; j < binNum[i]; j++) {
            double x = 1e-3 / (36.31 * 0.05961 * 0.05961) / (bins[i][j + 1] - bins[i][j]);
            xSec[i][j] *= x;
            sta[i][j] *= x;
        }
    }
    // Systematic uncertainties
    // sys1 - branching ratio; sys2 - luminosity; sys3 - correction; sys4 - fitting; sys5 - lifetime variable
    const double sys1 = 0.0055, sys2 = 0.012, sys5 = 0.003;
    const vector<vector<double>> sys3 = {
        {0.1052, 0.1094, 0.1328, 0.1506, 0.1709, 0.1126},
        {0.0754, 0.0534, 0.0668, 0.1883, 0.1128, 0.1063},
    }, sys4 = {
        {0.0156, 0.0382, 0.0234, 0.0361, 0.0461, 0.0377},
        {0.0160, 0.0215, 0.0338, 0.0339, 0.0080, 0.0422},
    };
    // Theoretical predictions
    //// Central value (not used!)
    vector<vector<double>> predMid = {
        {139.735, 24.8229, 3.59431, 0.831184, 0.246979, 0.055994},//0.105625, 0.048838, 0.0135205},
        {36.1687, 11.3241, 2.71414, 1.3498, 1.11811, 1.48746},
    };
    //// Lower limit
    vector<vector<double>> predLow = {
        {88.8285, 14.7747, 1.89624, 0.382532, 0.105109, 0.022624},//0.0427407, 0.0198011, 0.00533146},
        {23.1626, 7.03733, 1.61862, 0.766264, 0.577607, 0.643568},
    };
    //// Upper limit
    vector<vector<double>> predUp = {
        {225.85, 42.2799, 7.03847, 1.85164, 0.600244, 0.13948},//0.2615, 0.122751, 0.0342023},
        {56.3082, 18.5907, 4.51577, 2.36708, 2.12218, 3.52344},
    };
    // Store event yields
    double frac = 0.660;
    vector<double> yVar = {93.0, 74.4};
    // Set plotting limit mannually
    vector<double> yMax = {1000, 1000}, yMin = {0.01, 0.1};
    // Read SPS and DPS data for later use
    TFile spsFile("WeightNLO.root", "READ"), dpsFile("UnweightDPS.root", "READ");
    TTree *spsTree = (TTree *)spsFile.Get("data"), *dpsTree = (TTree *)dpsFile.Get("data");
    int nSPSEtr = spsTree->GetEntries(), nDPSEtr = dpsTree->GetEntries();
    // Store histograms and p.d.f.s in arrays
    RooRealVar **var = new RooRealVar*[varNum];
    TH1D **h = new TH1D*[varNum], **hSPS = new TH1D*[varNum], **hDPS = new TH1D*[varNum], **hThy = new TH1D*[varNum];
    RooDataHist **dh = new RooDataHist*[varNum], **dhSPS = new RooDataHist*[varNum], **dhDPS = new RooDataHist*[varNum], **dhThy = new RooDataHist*[varNum];
    RooHistPdf **hpdfSPS = new RooHistPdf*[varNum], **hpdfDPS = new RooHistPdf*[varNum];
    Double_t *spsVar = new Double_t[varNum], *dpsVar = new Double_t[varNum];
    // Loop on all kinematic variables
    for(int i = 0; i < varNum; i++) {
        // Fill differential cross section results
        h[i] = new TH1D(("h_" + varName[i]).c_str(), ("h_" + varName[i]).c_str(), binNum[i], bins[i].data());
        for(int j = 0; j < binNum[i]; j++) {
            h[i]->Fill(bins[i][j] + 0.01, xSec[i][j]);
            h[i]->SetBinError(j + 1, sta[i][j] + xSec[i][j] * combError(sys1, sys2, sys3[i][j], sys4[i][j], sys5));
        }
        var[i] = new RooRealVar(varName[i].c_str(), varName[i].c_str(), bins[i][0], bins[i][binNum[i]]);
        dh[i] = new RooDataHist(("dh_" + varName[i]).c_str(), ("dh_" + varName[i]).c_str(), *(var[i]), h[i]);
        // Fill SPS and DPS events
        double nSPSEvt = 0, nSPSWgt = 0, nDPSEvt = 0, nDPSWgt = 0;
        spsTree->SetBranchAddress(varName[i].c_str(), &(spsVar[i]));
        dpsTree->SetBranchAddress(varName[i].c_str(), &(dpsVar[i]));
        hSPS[i] = new TH1D(("hSPS_" + varName[i]).c_str(), ("hSPS_" + varName[i]).c_str(), binNum[i], bins[i].data());
        hDPS[i] = new TH1D(("hDPS_ " + varName[i]).c_str(), ("hDPS_" + varName[i]).c_str(), binNum[i], bins[i].data());
        //// SPS
        for(int j = 0; j < nSPSEtr; j++) {
            spsTree->GetEntry(j);
            if(spsVar[i] < bins[i][0]) continue;
            nSPSEvt += 1;
            if(spsVar[i] > bins[i][binNum[i]]) continue;
            int pos = upper_bound(bins[i].begin(), bins[i].end(), spsVar[i]) - bins[i].begin();
            hSPS[i]->Fill(spsVar[i], 1.0 / (bins[i][pos] - bins[i][pos - 1]));
            nSPSWgt += 1.0 / (bins[i][pos] - bins[i][pos - 1]);
        }
        hSPS[i]->Scale(1.0 / nSPSWgt);
        dhSPS[i] = new RooDataHist(("dhSPS_" + varName[i]).c_str(), ("dhSPS_" + varName[i]).c_str(), *(var[i]), hSPS[i]);
        hpdfSPS[i] = new RooHistPdf(("hpdfSPS_" + varName[i]).c_str(), ("hpdfSPS_" + varName[i]).c_str(), *(var[i]), *(dhSPS[i]));
        //// DPS
        for(int j = 0; j < nDPSEtr; j++) {
            dpsTree->GetEntry(j);
            if(dpsVar[i] < bins[i][0]) continue;
            nDPSEvt += 1;
            if(dpsVar[i] > bins[i][binNum[i]]) continue;
            int pos = upper_bound(bins[i].begin(), bins[i].end(), dpsVar[i]) - bins[i].begin();
            hDPS[i]->Fill(dpsVar[i], 1.0 / (bins[i][pos] - bins[i][pos - 1]));
            nDPSWgt += 1.0 / (bins[i][pos] - bins[i][pos - 1]);
        }
        hDPS[i]->Scale(1.0 / nDPSWgt);
        dhDPS[i] = new RooDataHist(("dhDPS_" + varName[i]).c_str(), ("dhDPS_" + varName[i]).c_str(), *(var[i]), hDPS[i]);
        hpdfDPS[i] = new RooHistPdf(("hpdfDPS_" + varName[i]).c_str(), ("hpdfDPS_" + varName[i]).c_str(), *(var[i]), *(dhDPS[i]));
        // Fill theoretical predictions
        hThy[i] = new TH1D(("hThy_" + varName[i]).c_str(), ("hThy_" + varName[i]).c_str(), binNum[i], bins[i].data());
        for(int j = 0; j < binNum[i]; j++) {
            hThy[i]->Fill(bins[i][j] + 0.01, (predLow[i][j] + predUp[i][j]) / 2);
            hThy[i]->SetBinError(j + 1, (predUp[i][j] - predLow[i][j]) / 2);
        }
        dhThy[i] = new RooDataHist(("dhThy_" + varName[i]).c_str(), ("dhThy_" + varName[i]).c_str(), *(var[i]), hThy[i]);
        // Fit and plot
        string xTitle = varLatx[i]+(varUnit[i].length() ? ("("+varUnit[i]+")") : ""), yTitle = "d#sigma/d"+varLatx[i]+" ("+(varUnit[i].length() ? ("pb/"+varUnit[i]) : "pb")+")";
        plot_temp(varName[i], *(var[i]), dh[i], dhSPS[i], dhDPS[i], dhThy[i], frac, yVar[i], xTitle, yTitle, yMax[i], yMin[i]);
    }
    return;
}
