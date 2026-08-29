#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TTree.h"
#include "RooAbsPdf.h"
#include "RooExtendPdf.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include <iostream>
#include <fstream>

using namespace std;
using namespace RooFit;

void Fit_pull() {
    ifstream inFile("summary.txt");
    double sig, sigErr, x;
    TTree *t = new TTree("tree", "tree");
    t->Branch("x", &x);
    if(!inFile.is_open()) {
        cerr << "Cannot open summary.txt" << endl;
        return;
    }
    for(int i = 0; i < 5000 && (inFile >> sig >> sigErr); i++) {
        if(sig < 0 || sigErr <= 0) continue;
        x = (sig - 4797.8161) / sigErr;
        t->Fill();
    }
    inFile.close();
    if(t->GetEntries() == 0) {
        cerr << "No valid pull entries found" << endl;
        return;
    }

    RooRealVar n_P_P("x", "x", -4, 4);
    RooDataSet *data = new RooDataSet("data", "data", t, RooArgSet(n_P_P));

    RooRealVar mean("mean", "mean", 0, -1, 1);
    RooRealVar sigma("sigma", "sigma", 1, 0.2, 5);
    RooGaussian gauss("gauss", "gauss", n_P_P, mean, sigma);

    RooRealVar n("n", "n", 1e3, 1, 1e4);
    RooExtendPdf pdf("pdf", "pdf", gauss, n);
    // RooDataSet *data = pdf.generate(x, Extended(kTRUE));
    RooFitResult *res = 0;
    for(int attempt = 0; attempt < 5; ++attempt) {
        delete res;
        res = pdf.fitTo(*data, Save());
        if(res && !res->status()) break;
    }
    if(!res || res->status()) {
        cerr << "Pull fit failed after five attempts" << endl;
        delete res;
        return;
    }

    TCanvas *c = new TCanvas("c", "c", 2000, 1800);
    RooPlot *f = n_P_P.frame(RooFit::Title("n_P_P Pull Distribution;N_{corr}"), Bins(40));
    data->plotOn(f, DataError(RooAbsData::SumW2), Name("Data"));//
    pdf.plotOn(f, LineColor(kBlue), LineWidth(2), Name("GS"));
    TLegend *l = new TLegend(.65, .60, .85, .85);
    l->AddEntry(f->findObject("Data"), "N_{corr} Pull", "LEP");
    l->AddEntry(f->findObject("GS"), "Gaussian", "L");
    f->Draw();
    l->DrawClone();
    gSystem->mkdir("fig", kTRUE);
    c->SaveAs("fig/pull.png");
    pdf.getVariables()->Print("v");
    delete res;
    return;
}
