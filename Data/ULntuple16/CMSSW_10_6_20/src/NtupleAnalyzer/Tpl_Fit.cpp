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
void plot_temp(string varName, RooRealVar &var, RooDataHist *dh, RooDataHist *dhSPS, RooDataHist *dhDPS, double frac, double n_sps_dps, string xTitle, string yTitle) {
    TCanvas *canvas = new TCanvas("canvas", "canvas", 1500, 1300);
    RooPlot *frame = var.frame(), *frame_pull = var.frame(RooFit::Title("Pull"));
    RooHistFunc func_SPS("func_SPS", "func_SPS", var, *dhSPS), func_DPS("func_DPS", "func_DPS", var, *dhDPS);
    RooRealVar coef_SPS("coef_SPS", "coef_SPS", frac * n_sps_dps);
    RooRealVar coef_DPS("coef_DPS", "coef_DPS", (1 - frac) * n_sps_dps);
    RooRealSumFunc func("func", "func", {func_SPS, func_DPS}, {coef_SPS, coef_DPS});
    dh->plotOn(frame, Name("Data"), LineColor(kBlack), LineWidth(2), MarkerSize(1), MarkerStyle(20));//, DrawOption("E1")
    func.plotOn(frame, LineColor(kBlue), Name("All")); 
    func_SPS.plotOn(frame, Normalization(coef_SPS.getVal()), LineColor(kRed), LineStyle(kDashed), LineWidth(2), Name("SPS"));
    func_DPS.plotOn(frame, Normalization(coef_DPS.getVal()), LineColor(40), DrawOption("F"), FillColor(40), MoveToBack(), Name("DPS"));
    RooHist* pull = frame->pullHist("Data", "All");
    frame_pull->addPlotable(pull, "P");
    TLegend *legend = new TLegend(.55, .60, .75, .85);
    legend->AddEntry(frame->findObject("Data"), "RunII 2016", "L");
    legend->AddEntry(frame->findObject("All"), "Total p.d.f.", "L");
    legend->AddEntry(frame->findObject("SPS"), "SPS p.d.f.", "L");
    legend->AddEntry(frame->findObject("DPS"), "DPS p.d.f.", "F");
    canvas->Divide(1, 2);
    canvas->cd(1)->SetPad(0.01, 0.20, 0.99, 0.99);
    gPad->SetLeftMargin(0.15);
    frame->SetTitle("");
    // frame->SetXTitle(xTitle.c_str());
    frame->SetYTitle(yTitle.c_str());
    frame->Draw();
    legend->DrawClone();
    TLatex latex;
    latex.SetNDC();
    add_latex(latex);
    canvas->cd(2)->SetPad(0.01, 0.03, 0.99, 0.25);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.25);
    frame_pull->SetTitle("");
    frame_pull->GetYaxis()->SetTitleSize(0.1);
    frame_pull->GetYaxis()->SetLabelSize(0.1);
    frame_pull->GetXaxis()->SetTitleSize(0.1);
    frame_pull->GetXaxis()->SetLabelSize(0.1);
    frame_pull->GetXaxis()->SetTitle(xTitle.c_str());
    frame_pull->Draw();
    TLine* l = new TLine(frame->GetXaxis()->GetXmin(), 0, frame->GetXaxis()->GetXmax(), 0);
    l->SetLineWidth(2);
    l->SetLineColor(kRed);
    l->SetLineStyle(kDashed);
    l->Draw("same");
    canvas->SaveAs(("fig/temp/Template_" + varName + ".png").c_str());
}

void Tpl_Fit() {
    // Names of kinematic variables
    const int varNum = 5;
    const string varName[] = {"delta_y", "delta_phi", "evt_mass", "evt_y", "evt_pt"};
    const string varLatx[] = {"|#Delta y(J/#psi_{1},J/#psi_{2})|", "|#Delta#phi(J/#psi_{1},J/#psi_{2})|", "M(J/#psi_{1},J/#psi_{2})", "|y(J/#psi_{1},J/#psi_{2})|", "p_{T}(J/#psi_{1},J/#psi_{2})"};
    const string varUnit[] = {"", "", "GeV", "", "GeV"};
    // Binning configurations
    const int binNum[] = {6, 8, 7, 5, 9};
    const vector<vector<double>> bins = {
        {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 4.0},
        {0.0_PI, 0.125_PI, 0.25_PI, 0.375_PI, 0.5_PI, 0.625_PI, 0.75_PI, 0.875_PI, 1.0_PI},
        {7.5, 17.5, 27.5, 37.5, 47.5, 57.5, 67.5, 107.5},
        // {7.5, 15.5, 23.5, 31.5, 39.5, 47.5, 57.5, 67.5, 107.5}
        {0.0, 0.4, 0.8, 1.2, 1.6, 2.0},
        {0, 5, 10, 15, 20, 25, 30, 35, 40, 80}
    };
    // Differential cross sections(event yield)
    vector<vector<double>> xSec = {
        {2585.85, 848.6, 514.317, 390.036, 264.352, 207.22},
        {1292.21, 767.777, 223.793, 336.095, 285.619, 329.894, 627.668, 993.237},
        {2113.59, 1038.69, 910.428, 397.103, 197.886, 47.2291, 106.221},
        // {916.264, 907.639, 1143.22, 836.383, 436.256, 175.983, 179.425},
        {1205.43, 1293.12, 1244.7, 785.147, 272.65},
        {654.707, 663.187, 537.083, 327.398, 1140.69, 782.025, 260.75, 233.562, 248.597}
    };
    // Statistical uncertainties(error of event yield)
    vector<vector<double>> sta = {
        {54.3052, 32.3728, 25.6471, 22.7905, 18.318, 16.6171},
        {37.797, 29.4854, 16.5064, 20.0209, 18.4893, 20.3366, 28.0881, 35.7158},
        {48.3568, 35.2975, 33.7636, 23.0093, 16.1421, 9.88096, 12.4474},
        // {33.2635, 32.3605, 37.5913, 31.8681, 22.7294, 15.3576, 14.8315},
        {38.806, 39.6443, 38.6797, 30.1464, 17.5307},
        {28.1011, 29.4673, 25.7194, 20.2752, 35.8797, 30.6585, 17.2883, 16.1139, 16.7508}
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
        {0.1495, 0.0976, 0.0721, 0.0839, 0.1499, 0.1265, 0.1399, 0.1049},
        {0.0671, 0.1035, 0.1459, 0.1701, 0.2337, 0.3067, 0.2037},
        // {0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0.0800, 0.1149, 0.1111, 0.1080, 0.0990},
        {0.0905, 0.1029, 0.1098, 0.1564, 0.0721, 0.1291, 0.1244, 0.1684, 0.3094}
    }, sys4 = {
        {0.0156, 0.0382, 0.0234, 0.0361, 0.0461, 0.0377},
        {0.0131, 0.0181, 0.0331, 0.0314, 0.0253, 0.0093, 0.0304, 0.0321},
        {0.0149, 0.0256, 0.0191, 0.0613, 0.0378, 0.1058, 0.0317},
        // {0.0671, 0.1035, 0.1459, 0.1701, 0.2337, 0.3067, 0.2037},
        {0.0372, 0.0231, 0.0241, 0.0173, 0.0001},
        {0.0281, 0.0300, 0.0190, 0.0305, 0.0169, 0.0311, 0.0137, 0.0161, 0.0205}
    };
    // Specify which variables participate in simultaneous fit
    const bool simul[] = {true, true, true, true, false};
    // Store results of f_SPS
    double *fSPS = new double[varNum];
    // Read SPS and DPS data for later use
    TFile spsFile("WeightNLO.root", "READ"), dpsFile("UnweightDPS.root", "READ");
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
        // Construct template p.d.f.
        RooRealVar frac("frac", "frac", 0.5, 0, 1);
        RooRealVar n_sps_dps("n_sps_dps", "n_sps_dps", 10, 1, 1e4);
        RooAddPdf pdf_sps_dps("pdf_sps_dps", "pdf_sps_dps", RooArgList(*(hpdfSPS[i]), *(hpdfDPS[i])), frac);
        RooExtendPdf pdf_all("pdf_all", "pdf_all", pdf_sps_dps, n_sps_dps);
        // Fit and plot
        pdf_all.fitTo(*(dh[i]));
        string xTitle = varLatx[i]+(varUnit[i].length() ? ("("+varUnit[i]+")") : ""), yTitle = "d#sigma/d"+varLatx[i]+" ("+(varUnit[i].length() ? ("pb/"+varUnit[i]) : "pb")+")";
        plot_temp(varName[i], *(var[i]), dh[i], dhSPS[i], dhDPS[i], frac.getVal(), n_sps_dps.getVal(), xTitle, yTitle);
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
            h2->SetBinError(j + 1, sta[i][j] + xSec[i][j] * combError(sys1, sys2, sys3[i][j], sys4[i][j], sys5));
        }
        TCanvas *canvas2 = new TCanvas("canvas2", "canvas2", 1000, 1000);
        gPad->SetLeftMargin(0.15);
        h1->SetTitle("");
        h1->SetLineWidth(4);
        h1->SetFillColor(40);
        h1->GetYaxis()->SetRangeUser(0, 1.3 * h1->GetMaximum());
        h1->GetXaxis()->SetTitle(xTitle.c_str());
        h1->GetYaxis()->SetTitle(yTitle.c_str());
        h1->SetStats(0);
        h1->Draw("E2");
        h2->SetLineColor(kBlack);
        h2->SetLineWidth(2);
        h2->SetMarkerStyle(20);
        h2->SetMarkerSize(1);
        h2->SetStats(0);
        h2->Draw("E same");
        TLatex latex;
        latex.SetNDC();
        add_latex(latex);
        canvas2->SaveAs(("fig/temp/Xsec_" + varName[i] + ".png").c_str());
    }
    // Plot systematic uncertainty summary
    for(int i = 0; i < varNum; i++) {
        TH1D *hColl = new TH1D[6];
        for(int j = 0; j <= 5; j++) {
            string s = "sys" + to_string(j) + "_" + varName[i];
            hColl[j] = TH1D(s.c_str(), s.c_str(), binNum[i], bins[i].data());
        }
        // TH1D *h0 = new TH1D(("sys_" + varName[i]).c_str(), ("sys_" + varName[i]).c_str(), binNum[i], bins[i].data());
        // TH1D *h1 = new TH1D(("sys1_" + varName[i]).c_str(), ("sys1_" + varName[i]).c_str(), binNum[i], bins[i].data());
        // TH1D *h2 = new TH1D(("sys2_" + varName[i]).c_str(), ("sys2_" + varName[i]).c_str(), binNum[i], bins[i].data());
        // TH1D *h3 = new TH1D(("sys3_" + varName[i]).c_str(), ("sys3_" + varName[i]).c_str(), binNum[i], bins[i].data());
        // TH1D *h4 = new TH1D(("sys4_" + varName[i]).c_str(), ("sys4_" + varName[i]).c_str(), binNum[i], bins[i].data());
        // TH1D *h5 = new TH1D(("sys5_" + varName[i]).c_str(), ("sys5_" + varName[i]).c_str(), binNum[i], bins[i].data());
        for(int j = 0; j < binNum[i]; j++) {
            hColl[0].SetBinContent(j + 1, 100 * combError(sys1, sys2, sys3[i][j], sys4[i][j], sys5));
            hColl[1].SetBinContent(j + 1, 100 * sys1);
            hColl[2].SetBinContent(j + 1, 100 * sys2);
            hColl[3].SetBinContent(j + 1, 100 * sys3[i][j]);
            hColl[4].SetBinContent(j + 1, 100 * sys4[i][j]);
            hColl[5].SetBinContent(j + 1, 100 * sys5);
        }
        TCanvas *canvas = new TCanvas("canvas2", "canvas2", 1000, 1000);
        hColl[0].SetTitle("");
        hColl[0].SetStats(0);
        hColl[0].SetLineWidth(2);
        hColl[0].SetLineColor(kBlack);
        hColl[1].SetLineColor(kRed);
        hColl[2].SetLineColor(kBlue);
        hColl[3].SetLineColor(kGreen);
        hColl[4].SetLineColor(kMagenta);
        hColl[5].SetLineColor(kYellow);
        string xTitle = varLatx[i]+(varUnit[i].length() ? ("("+varUnit[i]+")") : ""), yTitle = "d#sigma/d"+varLatx[i]+" ("+(varUnit[i].length() ? ("pb/"+varUnit[i]) : "pb")+")";
        hColl[0].GetXaxis()->SetTitle(xTitle.c_str());
        hColl[0].GetYaxis()->SetTitle(yTitle.c_str());
        hColl[0].GetYaxis()->SetRangeUser(0, 40);
        hColl[0].Draw();
        for(int j = 1; j <= 5; j++) hColl[j].Draw("same");
        TLatex latex;
        latex.SetNDC();
        add_latex(latex);
        TLegend *legend = new TLegend(.65, .50, .9, .9);
        legend->AddEntry(&hColl[0], "Total systematic", "L");
        legend->AddEntry(&hColl[1], "Branching ratio", "L");
        legend->AddEntry(&hColl[2], "Luminosity", "L");
        legend->AddEntry(&hColl[3], "Correction", "L");
        legend->AddEntry(&hColl[4], "Fitting stability", "L");
        legend->AddEntry(&hColl[5], "Time variable", "L");
        legend->DrawClone();
        canvas->SaveAs(("fig/temp/Sys_" + varName[i] + ".png").c_str());
        delete [] hColl;
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
        string xTitle = varLatx[i]+(varUnit[i].length() ? ("("+varUnit[i]+")") : ""), yTitle = "d#sigma/d"+varLatx[i]+" ("+(varUnit[i].length() ? ("pb/"+varUnit[i]) : "pb")+")";
        plot_temp(varName[i] + "_simul", *(var[i]), dh[i], dhSPS[i], dhDPS[i], frac.getVal(), n_var[i]->getVal(), xTitle, yTitle);
    }
    simPdf.getVariables()->Print("v");
    // Output to screen
    for(int i = 0; i < varNum; i++) {
        cout<<varName[i]<<": fSPS="<<setprecision(3)<<fSPS[i]<<endl;
        for(int j = 0; j < binNum[i]; j++) cout<<xSec[i][j]<<"\\pm"<<sta[i][j]<<"\\pm"<<(xSec[i][j] * combError(sys1, sys2, sys3[i][j], sys4[i][j], sys5))<<endl;
        cout<<endl;
    }
    cout<<endl;
    return;
}