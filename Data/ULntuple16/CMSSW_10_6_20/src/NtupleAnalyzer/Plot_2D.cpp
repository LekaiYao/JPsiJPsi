#include "Plot_4D.hpp"
#include "TH2D.h"
#include "TSystem.h"

void Plot_2D(string var, double vmin, double vmax) {
    // Define variables
    RooRealVar Jpsi_mass1("Jpsi_mass1", "Jpsi_mass1", 2.95, 3.25);
    RooRealVar Jpsi_mass2("Jpsi_mass2", "Jpsi_mass2", 2.95, 3.25);
    // RooRealVar Jpsi_ctau1("Jpsi_ctau1", "Jpsi_ctau1", -0.03, 0.16);
    // RooRealVar Jpsi_ctau2("Jpsi_ctau2", "Jpsi_ctau2", -0.03, 0.16);
    RooRealVar evt_weight("evt_weight", "evt_weight", 0, 1000);
    RooRealVar kine(var.c_str(), var.c_str(), vmin, vmax);
    RooArgSet variables;
    variables.add(Jpsi_mass1);
    variables.add(Jpsi_mass2);
    // variables.add(Jpsi_ctau1);
    // variables.add(Jpsi_ctau2);
    variables.add(evt_weight);
    variables.add(kine);
    // Read file and apply cut
    string sel = var + " > " + to_string(vmin) + " && " + var + " < " + to_string(vmax);
    TFile *dataFile = new TFile("WeightData.root", "READ");
    if(!dataFile || dataFile->IsZombie()) {
        cerr << "Cannot open WeightData.root" << endl;
        return;
    }
    TTree *dataTree = (TTree*)dataFile->Get("data");
    if(!dataTree) {
        cerr << "Cannot find tree 'data' in WeightData.root" << endl;
        return;
    }
    RooDataSet *data = new RooDataSet("data", "data", dataTree, variables, sel.c_str(), "evt_weight");

    // Parameters below are imported from total cross section 4D fit
    TFile *f = new TFile("Model_4D_tot.root");
    if(!f || f->IsZombie()) {
        cerr << "Cannot open Model_4D_tot.root" << endl;
        return;
    }
    RooWorkspace *wsp = (RooWorkspace *)f->Get("wsp");
    RooAddPdf *pdf_all = wsp ? dynamic_cast<RooAddPdf *>(wsp->pdf("pdf_all")) : 0;
    if(!pdf_all) {
        cerr << "Cannot find RooAddPdf 'pdf_all'" << endl;
        return;
    }
    // RooAbsPdf &pdf_P_P = *(wsp->pdf("pdf_P_P"));
    // RooAbsPdf &pdf_P_NP = *(wsp->pdf("pdf_P_NP"));
    // RooAbsPdf &pdf_NP_P = *(wsp->pdf("pdf_NP_P"));
    // RooAbsPdf &pdf_NP_NP = *(wsp->pdf("pdf_NP_NP"));
    // RooAbsPdf &pdf_Sig_Comb = *(wsp->pdf("pdf_Sig_Comb"));
    // RooAbsPdf &pdf_Comb_Sig = *(wsp->pdf("pdf_Comb_Sig"));
    // RooAbsPdf &pdf_Comb_Comb = *(wsp->pdf("pdf_Comb_Comb"));

    // Draw data point and p.d.f. curve
    string prefix = "fig/2d/";
    TCanvas *canvas = new TCanvas("canvas", "canvas", 2000, 1000);
    TH2D *hPdf = (TH2D *)pdf_all->createHistogram("hPdf", Jpsi_mass1, RooFit::YVar(Jpsi_mass2));
    TH2D *hData = (TH2D *)data->createHistogram(Jpsi_mass1, Jpsi_mass2, 100, 100, "", "hData");
    canvas->Divide(2, 1);
    canvas->cd(1);
    hPdf->Draw("surf");
    canvas->cd(2);
    hData->Draw("lego");
    gSystem->mkdir(prefix.c_str(), kTRUE);
    canvas->SaveAs((prefix + "2D_JpsiMass12.png").c_str());
    return;
}
