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
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <string>
#include <vector>

using namespace std;
using namespace RooFit;

double combError(double sys1, double sys2, double sys3, double sys4) {
    return sqrt(sys1 * sys1 + sys2 * sys2 + sys3 * sys3 + sys4 * sys4);
}

bool parseRangeLine(const string &line, string &var, double &xmin, double &xmax) {
    size_t colon = line.find(':');
    size_t tilde = line.find('~');
    if(colon == string::npos || tilde == string::npos || tilde < colon) return false;
    var = line.substr(0, colon);
    xmin = atof(line.substr(colon + 1, tilde - colon - 1).c_str());
    xmax = atof(line.substr(tilde + 1).c_str());
    return true;
}

bool parseYieldLine(const string &line, double &xSec, double &sta, double &sys4MaxAbsFrac) {
    if(line.find("Event yield:") != 0) return false;
    size_t posYield = line.find(':');
    size_t posPlus = line.find("+/-");
    size_t posSta = line.find("(Sta.)");
    if(posYield == string::npos || posPlus == string::npos || posSta == string::npos) return false;
    xSec = atof(line.substr(posYield + 1, posPlus - posYield - 1).c_str());
    sta = atof(line.substr(posPlus + 3, posSta - (posPlus + 3)).c_str());
    double maxAbs = 0.0;
    size_t start = posSta;
    for(int k = 0; k < 3; k++) {
        size_t pct = line.find('%', start);
        if(pct == string::npos) return false;
        size_t comma = line.rfind(',', pct);
        size_t valStart = (comma == string::npos) ? start : comma + 1;
        double v = atof(line.substr(valStart, pct - valStart).c_str());
        maxAbs = max(maxAbs, fabs(v));
        start = pct + 1;
    }
    sys4MaxAbsFrac = maxAbs / 100.0;
    return true;
}

void plot_temp(const string &name, RooRealVar &var, RooDataHist *dh, RooDataHist *dhSPS, RooDataHist *dhDPS, double frac, double nTot) {
    TCanvas *canvas = new TCanvas(("c_" + name).c_str(), ("c_" + name).c_str(), 1500, 1300);
    RooPlot *frame = var.frame();
    RooHistFunc funcSPS(("funcSPS_" + name).c_str(), ("funcSPS_" + name).c_str(), var, *dhSPS);
    RooHistFunc funcDPS(("funcDPS_" + name).c_str(), ("funcDPS_" + name).c_str(), var, *dhDPS);
    RooRealVar coefSPS(("coefSPS_" + name).c_str(), ("coefSPS_" + name).c_str(), frac * nTot);
    RooRealVar coefDPS(("coefDPS_" + name).c_str(), ("coefDPS_" + name).c_str(), (1.0 - frac) * nTot);
    RooRealSumFunc func(("func_" + name).c_str(), ("func_" + name).c_str(), {funcSPS, funcDPS}, {coefSPS, coefDPS});
    dh->plotOn(frame, Name("Data"));
    func.plotOn(frame, LineColor(kBlack), Name("All"));
    funcSPS.plotOn(frame, Normalization(coefSPS.getVal()), LineColor(kBlue), LineStyle(kDashed), Name("SPS"));
    funcDPS.plotOn(frame, Normalization(coefDPS.getVal()), LineColor(kRed), LineStyle(kDotted), Name("DPS"));
    TLegend *legend = new TLegend(.65, .60, .85, .85);
    legend->AddEntry(frame->findObject("Data"), "Data xSec", "L");
    legend->AddEntry(frame->findObject("All"), "Total p.d.f.", "L");
    legend->AddEntry(frame->findObject("SPS"), "SPS p.d.f.", "L");
    legend->AddEntry(frame->findObject("DPS"), "DPS p.d.f.", "L");
    frame->SetXTitle(var.GetName());
    frame->SetYTitle((string("d#sigma/d(") + var.GetName() + ")").c_str());
    frame->Draw();
    legend->DrawClone();
    canvas->SaveAs(("fig/temp2/Template_" + name + ".png").c_str());
}

void Tpl_Fit() {
    const int varNum = 5;
    const string varName[] = {"delta_y", "delta_phi", "evt_mass", "evt_y", "evt_pt"};
    vector<vector<double>> bins(varNum), xSec(varNum), sta(varNum), sys4(varNum);

    ifstream inFile("log.txt");
    if(!inFile.is_open()) {
        cout<<"Cannot open log.txt"<<endl;
        return;
    }
    string line, curVar;
    double curMin = 0, curMax = 0;
    while(getline(inFile, line)) {
        if(line.empty()) continue;
        string parsedVar;
        double xmin = 0, xmax = 0;
        if(parseRangeLine(line, parsedVar, xmin, xmax)) {
            curVar = parsedVar;
            curMin = xmin;
            curMax = xmax;
            continue;
        }
        double y = 0, s = 0, sys4max = 0;
        if(!parseYieldLine(line, y, s, sys4max)) continue;
        int idx = -1;
        for(int i = 0; i < varNum; i++) if(curVar == varName[i]) idx = i;
        if(idx < 0) continue;
        if(bins[idx].empty()) bins[idx].push_back(curMin);
        bins[idx].push_back(curMax);
        xSec[idx].push_back(y);
        sta[idx].push_back(s);
        sys4[idx].push_back(sys4max);
    }
    inFile.close();

    vector<int> binNum(varNum, 0);
    for(int i = 0; i < varNum; i++) {
        binNum[i] = (int)xSec[i].size();
        if(binNum[i] == 0 || (int)bins[i].size() != binNum[i] + 1 || (int)sta[i].size() != binNum[i] || (int)sys4[i].size() != binNum[i]) {
            cout<<"Invalid or incomplete data for "<<varName[i]<<endl;
            return;
        }
    }

    const double sys1 = 0.015, sys2 = 0.026;
    vector<vector<double>> sys3(varNum);
    for(int i = 0; i < varNum; i++) sys3[i].assign(binNum[i], 0.0);
    for(int i = 0; i < varNum; i++) {
        for(int j = 0; j < binNum[i]; j++) {
            double scale = 1e-3 / (36.3 * 0.05961 * 0.05961) / (bins[i][j + 1] - bins[i][j]);
            xSec[i][j] *= scale;
            sta[i][j] *= scale;
        }
    }

    TFile spsFile("WeightSPS.root", "READ"), dpsFile("WeightDPS.root", "READ");
    TTree *spsTree = (TTree *)spsFile.Get("data"), *dpsTree = (TTree *)dpsFile.Get("data");
    if(!spsTree || !dpsTree) {
        cout<<"Cannot find tree data in WeightSPS.root or WeightDPS.root"<<endl;
        return;
    }

    const int fitIdx[2] = {0, 1}; // delta_y, delta_phi
    RooRealVar *var[2];
    TH1D *hData[2], *hSPS[2], *hDPS[2];
    RooDataHist *dhData[2], *dhSPS[2], *dhDPS[2];
    RooHistPdf *pdfSPS[2], *pdfDPS[2];
    double fSPS1D[2] = {0, 0};

    for(int k = 0; k < 2; k++) {
        int i = fitIdx[k];
        var[k] = new RooRealVar(varName[i].c_str(), varName[i].c_str(), bins[i][0], bins[i][binNum[i]]);
        hData[k] = new TH1D(("hData_" + varName[i]).c_str(), ("hData_" + varName[i]).c_str(), binNum[i], bins[i].data());
        hSPS[k] = new TH1D(("hSPS_" + varName[i]).c_str(), ("hSPS_" + varName[i]).c_str(), binNum[i], bins[i].data());
        hDPS[k] = new TH1D(("hDPS_" + varName[i]).c_str(), ("hDPS_" + varName[i]).c_str(), binNum[i], bins[i].data());

        for(int j = 0; j < binNum[i]; j++) {
            hData[k]->Fill(bins[i][j] + 0.01, xSec[i][j]);
            double sysAbs = xSec[i][j] * combError(sys1, sys2, sys3[i][j], sys4[i][j]);
            hData[k]->SetBinError(j + 1, sqrt(sta[i][j] * sta[i][j] + sysAbs * sysAbs));
        }
        dhData[k] = new RooDataHist(("dhData_" + varName[i]).c_str(), ("dhData_" + varName[i]).c_str(), *var[k], hData[k]);

        Double_t vSPS = 0, wSPS = 0, vDPS = 0, wDPS = 0;
        spsTree->SetBranchAddress(varName[i].c_str(), &vSPS);
        spsTree->SetBranchAddress("evt_weight", &wSPS);
        dpsTree->SetBranchAddress(varName[i].c_str(), &vDPS);
        dpsTree->SetBranchAddress("evt_weight", &wDPS);

        double sumSPS = 0, sumDPS = 0;
        Long64_t nSPS = spsTree->GetEntries(), nDPS = dpsTree->GetEntries();
        for(Long64_t j = 0; j < nSPS; j++) {
            spsTree->GetEntry(j);
            if(vSPS < bins[i][0] || vSPS > bins[i][binNum[i]]) continue;
            int pos = upper_bound(bins[i].begin(), bins[i].end(), vSPS) - bins[i].begin();
            double bwInv = 1.0 / (bins[i][pos] - bins[i][pos - 1]);
            hSPS[k]->Fill(vSPS, wSPS * bwInv);
            sumSPS += wSPS * bwInv;
        }
        for(Long64_t j = 0; j < nDPS; j++) {
            dpsTree->GetEntry(j);
            if(vDPS < bins[i][0] || vDPS > bins[i][binNum[i]]) continue;
            int pos = upper_bound(bins[i].begin(), bins[i].end(), vDPS) - bins[i].begin();
            double bwInv = 1.0 / (bins[i][pos] - bins[i][pos - 1]);
            hDPS[k]->Fill(vDPS, wDPS * bwInv);
            sumDPS += wDPS * bwInv;
        }
        if(sumSPS > 0) hSPS[k]->Scale(1.0 / sumSPS);
        if(sumDPS > 0) hDPS[k]->Scale(1.0 / sumDPS);

        dhSPS[k] = new RooDataHist(("dhSPS_" + varName[i]).c_str(), ("dhSPS_" + varName[i]).c_str(), *var[k], hSPS[k]);
        dhDPS[k] = new RooDataHist(("dhDPS_" + varName[i]).c_str(), ("dhDPS_" + varName[i]).c_str(), *var[k], hDPS[k]);
        pdfSPS[k] = new RooHistPdf(("pdfSPS_" + varName[i]).c_str(), ("pdfSPS_" + varName[i]).c_str(), *var[k], *dhSPS[k]);
        pdfDPS[k] = new RooHistPdf(("pdfDPS_" + varName[i]).c_str(), ("pdfDPS_" + varName[i]).c_str(), *var[k], *dhDPS[k]);

        RooRealVar frac(("frac_1D_" + varName[i]).c_str(), ("frac_1D_" + varName[i]).c_str(), 0.5, 0, 1);
        RooRealVar nAll(("n_1D_" + varName[i]).c_str(), ("n_1D_" + varName[i]).c_str(), 100, 1, 1e6);
        RooAddPdf mixPdf(("mix_1D_" + varName[i]).c_str(), ("mix_1D_" + varName[i]).c_str(), RooArgList(*pdfSPS[k], *pdfDPS[k]), frac);
        RooExtendPdf extPdf(("ext_1D_" + varName[i]).c_str(), ("ext_1D_" + varName[i]).c_str(), mixPdf, nAll);
        extPdf.fitTo(*dhData[k]);
        fSPS1D[k] = frac.getVal();
        plot_temp(varName[i] + "_1D", *var[k], dhData[k], dhSPS[k], dhDPS[k], frac.getVal(), nAll.getVal());
    }

    RooCategory cate("cate", "cate");
    cate.defineType("delta_y");
    cate.defineType("delta_phi");
    RooArgList vars;
    vars.add(*var[0]);
    vars.add(*var[1]);
    map<string, RooDataHist*> hMap;
    hMap["delta_y"] = dhData[0];
    hMap["delta_phi"] = dhData[1];
    RooDataHist combData("combData", "combData", vars, cate, hMap);

    RooRealVar frac2D("frac2D", "frac2D", 0.5, 0, 1);
    RooRealVar nDY("nDY", "nDY", 100, 1, 1e6), nDPhi("nDPhi", "nDPhi", 100, 1, 1e6);
    RooAddPdf pdfDY("pdfDY", "pdfDY", RooArgList(*pdfSPS[0], *pdfDPS[0]), frac2D);
    RooAddPdf pdfDPhi("pdfDPhi", "pdfDPhi", RooArgList(*pdfSPS[1], *pdfDPS[1]), frac2D);
    RooExtendPdf extDY("extDY", "extDY", pdfDY, nDY), extDPhi("extDPhi", "extDPhi", pdfDPhi, nDPhi);
    RooSimultaneous simPdf("simPdf", "simPdf", cate);
    simPdf.addPdf(extDY, "delta_y");
    simPdf.addPdf(extDPhi, "delta_phi");
    simPdf.fitTo(combData);

    const double fSPS2D = frac2D.getVal();
    plot_temp("delta_y_2D", *var[0], dhData[0], dhSPS[0], dhDPS[0], fSPS2D, nDY.getVal());
    plot_temp("delta_phi_2D", *var[1], dhData[1], dhSPS[1], dhDPS[1], fSPS2D, nDPhi.getVal());

    double relSys = 0.0;
    if(fSPS2D != 0.0) {
        relSys = max(fabs(fSPS2D - fSPS1D[0]), fabs(fSPS2D - fSPS1D[1])) / fabs(fSPS2D);
    }

    ofstream out("fig/temp2/fSPS_summary.txt");
    out<<"f_SPS_2D = "<<fSPS2D<<"\n";
    out<<"f_SPS_1D_delta_y = "<<fSPS1D[0]<<"\n";
    out<<"f_SPS_1D_delta_phi = "<<fSPS1D[1]<<"\n";
    out<<"relative_sys = "<<relSys<<"\n";
    out.close();

    cout<<"f_SPS_2D = "<<fSPS2D<<endl;
    cout<<"f_SPS_1D_delta_y = "<<fSPS1D[0]<<endl;
    cout<<"f_SPS_1D_delta_phi = "<<fSPS1D[1]<<endl;
    cout<<"relative_sys = "<<relSys<<endl;
}
