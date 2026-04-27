#include "Plot_4D.hpp"

void Fit_4D_diff(string var, double vmin, double vmax) {
    // Define variables
    RooRealVar Jpsi_mass1("Jpsi_mass1", "Jpsi_mass1", 2.95, 3.25);
    RooRealVar Jpsi_mass2("Jpsi_mass2", "Jpsi_mass2", 2.95, 3.25);
    RooRealVar Jpsi_ctau1("Jpsi_ctau1", "Jpsi_ctau1", -0.03, 0.16);
    RooRealVar Jpsi_ctau2("Jpsi_ctau2", "Jpsi_ctau2", -0.03, 0.16);
    RooRealVar evt_weight("evt_weight", "evt_weight", 0, 1000);
    RooRealVar kine(var.c_str(), var.c_str(), vmin, vmax);
    RooArgSet variables;
    variables.add(Jpsi_mass1);
    variables.add(Jpsi_mass2);
    variables.add(Jpsi_ctau1);
    variables.add(Jpsi_ctau2);
    variables.add(evt_weight);
    variables.add(kine);
    // Read file and apply cut
    string sel = var + " > " + to_string(vmin) + " && " + var + " < " + to_string(vmax);
    TFile *dataFile = new TFile("WeightData.root", "READ");
    TTree *dataTree = (TTree*)dataFile->Get("data");
    RooDataSet *data = new RooDataSet("data", "data", dataTree, variables, sel.c_str(), "evt_weight");
    
    // Parameters below are imported from total cross section 4D fit
    // All parameters but psi2S_a, psi2S_b are fixed already
    // Event yields are reconfigured
    TFile *f = new TFile("Model_4D_tot.root");
    RooWorkspace *wsp = (RooWorkspace *)f->Get("wsp");
    RooAddPdf &pdf_all = dynamic_cast<RooAddPdf &>(wsp->allPdfs()["pdf_all"]);
    RooAbsPdf &pdf_P_P = *(wsp->pdf("pdf_P_P"));
    RooAbsPdf &pdf_P_NP = *(wsp->pdf("pdf_P_NP"));
    RooAbsPdf &pdf_NP_P = *(wsp->pdf("pdf_NP_P"));
    RooAbsPdf &pdf_NP_NP = *(wsp->pdf("pdf_NP_NP"));
    RooAbsPdf &pdf_Sig_Comb = *(wsp->pdf("pdf_Sig_Comb"));
    RooAbsPdf &pdf_Comb_Sig = *(wsp->pdf("pdf_Comb_Sig"));
    RooAbsPdf &pdf_Comb_Comb = *(wsp->pdf("pdf_Comb_Comb"));
    wsp->var("n_P_P")->setVal(1e3);
    wsp->var("n_P_NP")->setVal(1e2);
    wsp->var("n_NP_NP")->setVal(1e2);
    wsp->var("n_Sig_Comb")->setVal(1e3);
    wsp->var("n_Comb_Comb")->setVal(1e2);
    wsp->var("n_P_P")->setMin(1);
    wsp->var("n_P_NP")->setMin(0);
    wsp->var("n_NP_NP")->setMin(0);
    wsp->var("n_Sig_Comb")->setMin(1);
    wsp->var("n_Comb_Comb")->setMin(1);

    // 4D fitting
    RooFitResult *res;
    while(true) {
        res = pdf_all.fitTo(*data->reduce(RooArgSet(Jpsi_mass1, Jpsi_mass2, Jpsi_ctau1, Jpsi_ctau2)), Save());
        if(!res->status() && res->edm()<0.01) break;
    }
    // Draw data point and p.d.f. curve
    string prefix = "fig/diff/";
    Plot_4D(data, pdf_all, pdf_P_P, pdf_P_NP, pdf_NP_P, pdf_NP_NP, pdf_Sig_Comb, pdf_Comb_Sig, pdf_Comb_Comb,
        prefix, Jpsi_mass1, Jpsi_mass2, Jpsi_ctau1, Jpsi_ctau2);
    pdf_all.getVariables()->Print("v");

    // Save parameters to file
    wsp->writeToFile("Model_4D_diff.root");
    return;
}