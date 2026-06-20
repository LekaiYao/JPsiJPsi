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
    const int binNum[] = {6, 5, 5, 5, 6};
    const vector<vector<double>> bins = {
        {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 4.0}, {0.0_PI, 0.2_PI, 0.4_PI, 0.6_PI, 0.8_PI, 1.0_PI}, {7.5, 17.5, 27.5, 37.5, 47.5, 107.5},
        {0.0, 0.4, 0.8, 1.2, 1.6, 2.0}
        , {0, 8, 16, 24, 32, 40, 80}
    };
    // Differential cross sections
    vector<vector<double>> xSec = {
        {632.201, 129.542, 72.1842, 83.3564, 12.3596, 19.4357}, {653.633, 60.9626, 46.7251, 60.4103, 138.908}, {796.426, 251.365, 236.285, 61.7321, 76.9698},
        {472.951, 296.483, 354.145, 222.521, 79.273}
        , {250.092, 204.513, 386.399, 357.97, 138.75, 100.456}
    };
    // Statistical uncertainties
    vector<vector<double>> sta = {
        {28.0054, 15.3417, 11.7034, 10.0083, 5.56255, 7.32121}, {30.6741, 9.6205, 9.08513, 10.8043, 16.3403}, {35.0803, 22.6198, 20.7751, 15.9625, 15.5374},
        {28.6379, 25.4596, 25.508, 19.7277, 13.0565}
        , {20.6896, 22.3792, 26.5792, 26.183, 15.1611, 12.4401}
    };
    for(int i = 0; i < varNum; i++) {
        for(int j = 0; j < binNum[i]; j++) {
            xSec[i][j] *= 1e-3 / (61.31 * 0.05961 * 0.008) / (bins[i][j + 1] - bins[i][j]);
            sta[i][j] *= 1e-3 / (61.31 * 0.05961 * 0.008) / (bins[i][j + 1] - bins[i][j]);
        }
    }
    // Systematic uncertainties
    const double sys1 = 0.0283, sys2 = 0.025;
    const vector<vector<double>> sys3 = {
        {0.00210645, 0.0180457, 0.0126121, 0.0126625, -0.00438329, 0.0115217}, {0.00959917, 0.0238948, 0.0133999, -0.00404647, 0.00467304},
        {0.0232711, 0.00839292, -0.000209678, 0.0291752, -0.00383268}, 
        {0.00551673, 0.00308563, 0.0243919, 0.00997864, 0.0170281}
        , {-0.000628592, 0.00329366, 0.013584, 0.0224671, -0.0165127, 0.0444984}
    }, sys4 = {
        {0.074069, 0.079177, 0.123122, 0.108043, 0.0917305, 0.433594}, {0.0469346, 0.0603078, 0.0501868, 0.0936282, 0.088843},
        {0.0501369, 0.106403, 0.104383, 0.0947623, 0.0210585}, 
        {0.136552, 0.0843728, 0.0157296, 0.0357353, 0.0430499}
        , {0.0432611, 0.0779885, 0.0144281, 0.0646486, 0.111898, 0.138225}
    };
    // Specify which variables participate in simultaneous fit
    const bool simul[] = {true, true, false, false, false};
    // Store results of f_SPS
    double *fSPS = new double[varNum];
    // Read SPS and DPS data for later use
    TFile spsFile("UnweightSPS.root", "READ"), dpsFile("UnweightDPS.root", "READ");
    TTree *spsTree = (TTree *)spsFile.Get("data"), *dpsTree = (TTree *)dpsFile.Get("data");
    int nSPSEtr = spsTree->GetEntries(), nDPSEtr = dpsTree->GetEntries();
    // Store histograms and p.d.f.s in arrays
    RooRealVar **var = new RooRealVar*[varNum];
    TH1D **h = new TH1D*[varNum], **hSPS = new TH1D*[varNum], **hDPS = new TH1D*[varNum];
    RooDataHist **dh = new RooDataHist*[varNum], **dhSPS = new RooDataHist*[varNum], **dhDPS = new RooDataHist*[varNum];
    RooHistPdf **hpdfSPS = new RooHistPdf*[varNum], **hpdfDPS = new RooHistPdf*[varNum];
    Double_t *spsVar = new Double_t[varNum], *dpsVar = new Double_t[varNum];
    // Loop on all kinematic variables
    for(int i = 0; i < varNum; i++) {
        // Fill differential cross section results
        h[i] = new TH1D(("h_" + varName[i]).c_str(), ("h_" + varName[i]).c_str(), binNum[i], bins[i].data());
        for(int j = 0; j < binNum[i]; j++) {
            h[i]->Fill(bins[i][j] + 0.01, xSec[i][j]);
            h[i]->SetBinError(j + 1, sta[i][j] + xSec[i][j] * combError(sys1, sys2, sys3[i][j], sys4[i][j]));
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
        // Construct template p.d.f.
        RooRealVar frac("frac", "frac", 0.5, 0, 1);
        RooRealVar n_sps_dps("n_sps_dps", "n_sps_dps", 10, 1, 1e4);
        RooAddPdf pdf_sps_dps("pdf_sps_dps", "pdf_sps_dps", RooArgList(*(hpdfSPS[i]), *(hpdfDPS[i])), frac);
        RooExtendPdf pdf_all("pdf_all", "pdf_all", pdf_sps_dps, n_sps_dps);
        // Fit and plot
        pdf_all.fitTo(*(dh[i]));
        plot_temp(varName[i], *(var[i]), dh[i], dhSPS[i], dhDPS[i], frac.getVal(), n_sps_dps.getVal());
        pdf_all.getVariables()->Print("v");
        fSPS[i] = 1.0 / (1.0 + (1.0 / frac.getVal() - 1.0) / (nSPSEvt / hSPS[i]->GetEntries()) * (nDPSEvt / hDPS[i]->GetEntries()));
        cout<<endl<<(nSPSEvt / hSPS[i]->GetEntries())<<' '<<(nDPSEvt / hDPS[i]->GetEntries())<<endl<<endl;
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
    // Simultaneous fit to specified variables
    //// Create category
    RooCategory cate("cate", "cate");
    RooArgList vars;
    map<string, RooDataHist*> hMap;
    for(int i = 0; i < varNum; i++) {
        if(!simul[i]) continue;
        cate.defineType(varName[i].c_str());
        hMap[varName[i]] = dh[i];
        vars.add(*(var[i]));
    }
    RooDataHist combData("combData", "combData", vars, cate, hMap);
    RooSimultaneous simPdf("simPdf", "simPdf", cate);
    //// Construct p.d.f.
    RooRealVar frac("frac", "frac", 0.5, 0, 1);
    RooRealVar **n_var = new RooRealVar*[varNum];
    RooAddPdf **pdf_var = new RooAddPdf*[varNum];
    RooExtendPdf **extpdf_var = new RooExtendPdf*[varNum];
    for(int i = 0; i < varNum; i++) {
        if(!simul[i]) continue;
        n_var[i] = new RooRealVar(("n_" + varName[i]).c_str(), ("n_" + varName[i]).c_str(), 100, 1, 1e4);
        pdf_var[i] = new RooAddPdf(("sim_" + varName[i]).c_str(), ("sim_" + varName[i]).c_str(), RooArgList(*(hpdfSPS[i]), *(hpdfDPS[i])), frac);
        extpdf_var[i] = new RooExtendPdf(("simex_" + varName[i]).c_str(), ("simex_" + varName[i]).c_str(), *(pdf_var[i]), *(n_var[i]));
        simPdf.addPdf(*(extpdf_var[i]), varName[i].c_str());
    }
    //// Template fit
    simPdf.fitTo(combData);
    for(int i = 0; i < varNum; i++) {
        if(!simul[i]) continue;
        plot_temp(varName[i] + "_simul", *(var[i]), dh[i], dhSPS[i], dhDPS[i], frac.getVal(), n_var[i]->getVal());
    }
    simPdf.getVariables()->Print("v");
    // Output to screen
    for(int i = 0; i < varNum; i++) {
        cout<<varName[i]<<": fSPS="<<fSPS[i]<<endl;
        for(int j = 0; j < binNum[i]; j++) cout<<xSec[i][j]<<"+/-"<<sta[i][j]<<"+/-"<<(xSec[i][j] * combError(sys1, sys2, sys3[i][j], sys4[i][j]))<<endl;
        cout<<endl;
    }
    cout<<endl;
    return;
}