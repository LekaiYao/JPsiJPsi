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
#include "RooFitResult.h"
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

void plot_temp3(const string &name, RooRealVar &var, RooDataHist *dh, RooDataHist *dhLO, RooDataHist *dhNLO, RooDataHist *dhDPS, double fLO, double fNLO, double nTot) {
    TCanvas *canvas = new TCanvas(("c_" + name).c_str(), ("c_" + name).c_str(), 1500, 1300);
    RooPlot *frame = var.frame();
    RooHistFunc funcLO(("funcLO_" + name).c_str(), ("funcLO_" + name).c_str(), var, *dhLO);
    RooHistFunc funcNLO(("funcNLO_" + name).c_str(), ("funcNLO_" + name).c_str(), var, *dhNLO);
    RooHistFunc funcDPS(("funcDPS_" + name).c_str(), ("funcDPS_" + name).c_str(), var, *dhDPS);
    double fDPS = 1.0 - fLO - fNLO;
    RooRealVar coefLO(("coefLO_" + name).c_str(), ("coefLO_" + name).c_str(), fLO * nTot);
    RooRealVar coefNLO(("coefNLO_" + name).c_str(), ("coefNLO_" + name).c_str(), fNLO * nTot);
    RooRealVar coefDPS(("coefDPS_" + name).c_str(), ("coefDPS_" + name).c_str(), fDPS * nTot);
    RooRealSumFunc func(("func_" + name).c_str(), ("func_" + name).c_str(), {funcLO, funcNLO, funcDPS}, {coefLO, coefNLO, coefDPS});
    dh->plotOn(frame, Name("Data"));
    func.plotOn(frame, LineColor(kBlack), Name("All"));
    funcLO.plotOn(frame, Normalization(coefLO.getVal()), LineColor(kBlue), LineStyle(kDashed), Name("LO"));
    funcNLO.plotOn(frame, Normalization(coefNLO.getVal()), LineColor(kGreen + 2), LineStyle(kDashDotted), Name("NLO"));
    funcDPS.plotOn(frame, Normalization(coefDPS.getVal()), LineColor(kRed), LineStyle(kDotted), Name("DPS"));
    TLegend *legend = new TLegend(.62, .56, .88, .85);
    legend->AddEntry(frame->findObject("Data"), "Data xSec", "L");
    legend->AddEntry(frame->findObject("All"), "Total p.d.f.", "L");
    legend->AddEntry(frame->findObject("LO"), "LO p.d.f.", "L");
    legend->AddEntry(frame->findObject("NLO"), "NLO p.d.f.", "L");
    legend->AddEntry(frame->findObject("DPS"), "DPS p.d.f.", "L");
    frame->SetXTitle(var.GetName());
    frame->SetYTitle((string("d#sigma/d(") + var.GetName() + ")").c_str());
    frame->Draw();
    legend->DrawClone();
    canvas->SaveAs(("fig/temp3/Template3_" + name + ".png").c_str());
}

void Tpl3_Fit() {
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

    TFile loFile("WeightSPS_pythia.root", "READ");
    TFile nloFile("WeightSPSNLO.root", "READ");
    TFile dpsFile("WeightDPS.root", "READ");
    TTree *loTree = (TTree *)loFile.Get("data");
    TTree *nloTree = (TTree *)nloFile.Get("data");
    TTree *dpsTree = (TTree *)dpsFile.Get("data");
    if(!loTree || !nloTree || !dpsTree) {
        cout<<"Cannot find tree data in WeightSPS_pythia.root / WeightSPSNLO.root / WeightDPS.root"<<endl;
        return;
    }

    const int fitIdx[2] = {0, 1}; // delta_y, delta_phi
    RooRealVar *var[2];
    TH1D *hData[2], *hLO[2], *hNLO[2], *hDPS[2];
    RooDataHist *dhData[2], *dhLO[2], *dhNLO[2], *dhDPS[2];
    RooHistPdf *pdfLO[2], *pdfNLO[2], *pdfDPS[2];
    double fLO1D[2] = {0, 0}, fNLO1D[2] = {0, 0}, fDPS1D[2] = {0, 0};

    for(int k = 0; k < 2; k++) {
        int i = fitIdx[k];
        var[k] = new RooRealVar(varName[i].c_str(), varName[i].c_str(), bins[i][0], bins[i][binNum[i]]);
        hData[k] = new TH1D(("hData_" + varName[i]).c_str(), ("hData_" + varName[i]).c_str(), binNum[i], bins[i].data());
        hLO[k] = new TH1D(("hLO_" + varName[i]).c_str(), ("hLO_" + varName[i]).c_str(), binNum[i], bins[i].data());
        hNLO[k] = new TH1D(("hNLO_" + varName[i]).c_str(), ("hNLO_" + varName[i]).c_str(), binNum[i], bins[i].data());
        hDPS[k] = new TH1D(("hDPS_" + varName[i]).c_str(), ("hDPS_" + varName[i]).c_str(), binNum[i], bins[i].data());

        for(int j = 0; j < binNum[i]; j++) {
            hData[k]->Fill(bins[i][j] + 0.01, xSec[i][j]);
            double sysAbs = xSec[i][j] * combError(sys1, sys2, sys3[i][j], sys4[i][j]);
            hData[k]->SetBinError(j + 1, sqrt(sta[i][j] * sta[i][j] + sysAbs * sysAbs));
        }
        dhData[k] = new RooDataHist(("dhData_" + varName[i]).c_str(), ("dhData_" + varName[i]).c_str(), *var[k], hData[k]);

        Double_t vLO = 0, wLO = 0, vNLO = 0, wNLO = 0, vD = 0, wD = 0;
        loTree->SetBranchAddress(varName[i].c_str(), &vLO);
        loTree->SetBranchAddress("evt_weight", &wLO);
        nloTree->SetBranchAddress(varName[i].c_str(), &vNLO);
        nloTree->SetBranchAddress("evt_weight", &wNLO);
        dpsTree->SetBranchAddress(varName[i].c_str(), &vD);
        dpsTree->SetBranchAddress("evt_weight", &wD);

        double sumLO = 0, sumNLO = 0, sumDPS = 0;
        Long64_t nLO = loTree->GetEntries(), nNLO = nloTree->GetEntries(), nDPS = dpsTree->GetEntries();
        for(Long64_t j = 0; j < nLO; j++) {
            loTree->GetEntry(j);
            if(vLO < bins[i][0] || vLO > bins[i][binNum[i]]) continue;
            int pos = upper_bound(bins[i].begin(), bins[i].end(), vLO) - bins[i].begin();
            double bwInv = 1.0 / (bins[i][pos] - bins[i][pos - 1]);
            hLO[k]->Fill(vLO, wLO * bwInv);
            sumLO += wLO * bwInv;
        }
        for(Long64_t j = 0; j < nNLO; j++) {
            nloTree->GetEntry(j);
            if(vNLO < bins[i][0] || vNLO > bins[i][binNum[i]]) continue;
            int pos = upper_bound(bins[i].begin(), bins[i].end(), vNLO) - bins[i].begin();
            double bwInv = 1.0 / (bins[i][pos] - bins[i][pos - 1]);
            hNLO[k]->Fill(vNLO, wNLO * bwInv);
            sumNLO += wNLO * bwInv;
        }
        for(Long64_t j = 0; j < nDPS; j++) {
            dpsTree->GetEntry(j);
            if(vD < bins[i][0] || vD > bins[i][binNum[i]]) continue;
            int pos = upper_bound(bins[i].begin(), bins[i].end(), vD) - bins[i].begin();
            double bwInv = 1.0 / (bins[i][pos] - bins[i][pos - 1]);
            hDPS[k]->Fill(vD, wD * bwInv);
            sumDPS += wD * bwInv;
        }
        if(sumLO > 0) hLO[k]->Scale(1.0 / sumLO);
        if(sumNLO > 0) hNLO[k]->Scale(1.0 / sumNLO);
        if(sumDPS > 0) hDPS[k]->Scale(1.0 / sumDPS);

        dhLO[k] = new RooDataHist(("dhLO_" + varName[i]).c_str(), ("dhLO_" + varName[i]).c_str(), *var[k], hLO[k]);
        dhNLO[k] = new RooDataHist(("dhNLO_" + varName[i]).c_str(), ("dhNLO_" + varName[i]).c_str(), *var[k], hNLO[k]);
        dhDPS[k] = new RooDataHist(("dhDPS_" + varName[i]).c_str(), ("dhDPS_" + varName[i]).c_str(), *var[k], hDPS[k]);
        pdfLO[k] = new RooHistPdf(("pdfLO_" + varName[i]).c_str(), ("pdfLO_" + varName[i]).c_str(), *var[k], *dhLO[k]);
        pdfNLO[k] = new RooHistPdf(("pdfNLO_" + varName[i]).c_str(), ("pdfNLO_" + varName[i]).c_str(), *var[k], *dhNLO[k]);
        pdfDPS[k] = new RooHistPdf(("pdfDPS_" + varName[i]).c_str(), ("pdfDPS_" + varName[i]).c_str(), *var[k], *dhDPS[k]);

        RooRealVar fLO(("fLO_1D_" + varName[i]).c_str(), ("fLO_1D_" + varName[i]).c_str(), 0.3, 0, 1);
        RooRealVar fNLO(("fNLO_1D_" + varName[i]).c_str(), ("fNLO_1D_" + varName[i]).c_str(), 0.3, 0, 1);
        RooFormulaVar fDPS(("fDPS_1D_" + varName[i]).c_str(), "1.0-@0-@1", RooArgList(fLO, fNLO));
        RooRealVar nAll(("n_1D_" + varName[i]).c_str(), ("n_1D_" + varName[i]).c_str(), 100, 1, 1e6);
        RooAddPdf mixPdf(("mix_1D_" + varName[i]).c_str(), ("mix_1D_" + varName[i]).c_str(), RooArgList(*pdfLO[k], *pdfNLO[k], *pdfDPS[k]), RooArgList(fLO, fNLO), true);
        RooExtendPdf extPdf(("ext_1D_" + varName[i]).c_str(), ("ext_1D_" + varName[i]).c_str(), mixPdf, nAll);
        extPdf.fitTo(*dhData[k]);
        fLO1D[k] = fLO.getVal();
        fNLO1D[k] = fNLO.getVal();
        fDPS1D[k] = fDPS.getVal();
        plot_temp3(varName[i] + "_1D", *var[k], dhData[k], dhLO[k], dhNLO[k], dhDPS[k], fLO.getVal(), fNLO.getVal(), nAll.getVal());
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

    RooRealVar fLO2D("fLO2D", "fLO2D", 0.3, 0, 1);
    RooRealVar fNLO2D("fNLO2D", "fNLO2D", 0.3, 0, 1);
    RooFormulaVar fDPS2D("fDPS2D", "1.0-@0-@1", RooArgList(fLO2D, fNLO2D));
    RooRealVar nDY("nDY", "nDY", 100, 1, 1e6), nDPhi("nDPhi", "nDPhi", 100, 1, 1e6);
    RooAddPdf pdfDY("pdfDY", "pdfDY", RooArgList(*pdfLO[0], *pdfNLO[0], *pdfDPS[0]), RooArgList(fLO2D, fNLO2D), true);
    RooAddPdf pdfDPhi("pdfDPhi", "pdfDPhi", RooArgList(*pdfLO[1], *pdfNLO[1], *pdfDPS[1]), RooArgList(fLO2D, fNLO2D), true);
    RooExtendPdf extDY("extDY", "extDY", pdfDY, nDY), extDPhi("extDPhi", "extDPhi", pdfDPhi, nDPhi);
    RooSimultaneous simPdf("simPdf", "simPdf", cate);
    simPdf.addPdf(extDY, "delta_y");
    simPdf.addPdf(extDPhi, "delta_phi");
    RooFitResult *fitRes2D = simPdf.fitTo(combData, Save(true));

    double fLO2DVal = fLO2D.getVal();
    double fNLO2DVal = fNLO2D.getVal();
    double fDPS2DVal = fDPS2D.getVal();
    double fDPS2DStatErr = 0.0;
    if(fitRes2D) {
        double errLO = fLO2D.getError();
        double errNLO = fNLO2D.getError();
        double cov = fitRes2D->correlation("fLO2D", "fNLO2D") * errLO * errNLO;
        double varDPS = errLO * errLO + errNLO * errNLO + 2.0 * cov;
        if(varDPS > 0.0) fDPS2DStatErr = sqrt(varDPS);
    }
    plot_temp3("delta_y_2D", *var[0], dhData[0], dhLO[0], dhNLO[0], dhDPS[0], fLO2DVal, fNLO2DVal, nDY.getVal());
    plot_temp3("delta_phi_2D", *var[1], dhData[1], dhLO[1], dhNLO[1], dhDPS[1], fLO2DVal, fNLO2DVal, nDPhi.getVal());

    double fDPS2DSysErr = max(fabs(fDPS2DVal - fDPS1D[0]), fabs(fDPS2DVal - fDPS1D[1]));
    double statPct = 0.0, sysPct = 0.0;
    if(fDPS2DVal != 0.0) {
        statPct = 100.0 * fDPS2DStatErr / fabs(fDPS2DVal);
        sysPct = 100.0 * fDPS2DSysErr / fabs(fDPS2DVal);
    }

    ofstream out("fig/temp3/fractions3_summary.txt");
    out<<"1D_delta_y f_LO = "<<fLO1D[0]<<"\n";
    out<<"1D_delta_y f_NLO = "<<fNLO1D[0]<<"\n";
    out<<"1D_delta_y f_DPS = "<<fDPS1D[0]<<"\n";
    out<<"1D_delta_phi f_LO = "<<fLO1D[1]<<"\n";
    out<<"1D_delta_phi f_NLO = "<<fNLO1D[1]<<"\n";
    out<<"1D_delta_phi f_DPS = "<<fDPS1D[1]<<"\n";
    out<<"2D f_LO = "<<fLO2DVal<<"\n";
    out<<"2D f_NLO = "<<fNLO2DVal<<"\n";
    out<<"2D f_DPS = "<<fDPS2DVal<<"\n";
    out<<"2D f_DPS_stat_err = "<<fDPS2DStatErr<<"\n";
    out<<"2D f_DPS_sys_err = "<<fDPS2DSysErr<<"\n";
    out<<"2D f_DPS_stat_pct = "<<statPct<<"\n";
    out<<"2D f_DPS_sys_pct = "<<sysPct<<"\n";
    out.close();

    cout<<"1D delta_y: f_LO="<<fLO1D[0]<<", f_NLO="<<fNLO1D[0]<<", f_DPS="<<fDPS1D[0]<<endl;
    cout<<"1D delta_phi: f_LO="<<fLO1D[1]<<", f_NLO="<<fNLO1D[1]<<", f_DPS="<<fDPS1D[1]<<endl;
    cout<<"2D: f_LO="<<fLO2DVal<<", f_NLO="<<fNLO2DVal<<", f_DPS="<<fDPS2DVal<<endl;
    cout<<"2D f_DPS_stat_err = "<<fDPS2DStatErr<<endl;
    cout<<"2D f_DPS_sys_err = "<<fDPS2DSysErr<<endl;
    cout<<"2D f_DPS_stat_pct = "<<statPct<<endl;
    cout<<"2D f_DPS_sys_pct = "<<sysPct<<endl;
}
