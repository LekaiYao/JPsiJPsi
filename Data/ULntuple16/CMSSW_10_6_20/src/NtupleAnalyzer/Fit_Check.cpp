#include "Plot_4D.hpp"

void Fit_Check(bool verbose) {
    // Read pseudo-data from file
    RooRealVar Jpsi_mass1("Jpsi_mass1", "Jpsi_mass1", 2.95, 3.25);
    RooRealVar Jpsi_mass2("Jpsi_mass2", "Jpsi_mass2", 2.95, 3.25);
    RooRealVar Jpsi_ctau1("Jpsi_ctau1", "Jpsi_ctau1", -0.03, 0.16);
    RooRealVar Jpsi_ctau2("Jpsi_ctau2", "Jpsi_ctau2", -0.03, 0.16);
    TFile *fData = new TFile("Data_4D.root");
    RooWorkspace *wk = (RooWorkspace *)fData->Get("wk");
    RooAbsPdf *model_all = wk->pdf("model_all");
    srand(time(0));
    RooRandom::randomGenerator()->SetSeed(rand());
    RooDataSet *data = model_all->generate(RooArgSet(Jpsi_mass1, Jpsi_mass2, Jpsi_ctau1, Jpsi_ctau2), Extended(kTRUE));
    // Read fit model from file
    TFile *fModel = new TFile("Model_4D_diff.root");
    RooWorkspace *wsp = (RooWorkspace *)fModel->Get("wsp");
    RooAddPdf &pdf_all = dynamic_cast<RooAddPdf &>(wsp->allPdfs()["pdf_all"]);
    RooAbsPdf &pdf_P_P = *(wsp->pdf("pdf_P_P"));
    RooAbsPdf &pdf_P_NP = *(wsp->pdf("pdf_P_NP"));
    RooAbsPdf &pdf_NP_P = *(wsp->pdf("pdf_NP_P"));
    RooAbsPdf &pdf_NP_NP = *(wsp->pdf("pdf_NP_NP"));
    RooAbsPdf &pdf_Sig_Comb = *(wsp->pdf("pdf_Sig_Comb"));
    RooAbsPdf &pdf_Comb_Sig = *(wsp->pdf("pdf_Comb_Sig"));
    RooAbsPdf &pdf_Comb_Comb = *(wsp->pdf("pdf_Comb_Comb"));
    double n_P_P = wsp->var("n_P_P")->getVal(), sta = wsp->var("n_P_P")->getError();
    // 4D fitting
    RooFitResult *res;
    while(true) {
        res = pdf_all.fitTo(*data, Save());
        if(!res->status()) break;
    }
    // Draw data point and p.d.f. curve
    string prefix = "fig/check/";
    Plot_4D(data, pdf_all, pdf_P_P, pdf_P_NP, pdf_NP_P, pdf_NP_NP, pdf_Sig_Comb, pdf_Comb_Sig, pdf_Comb_Comb,
        prefix, Jpsi_mass1, Jpsi_mass2, Jpsi_ctau1, Jpsi_ctau2);
    cout<<"Status: "<<res->status()<<endl;
    double sys4 = (wsp->var("n_P_P")->getVal() - n_P_P) / n_P_P * 100;// / wk->var("vScale")->getVal()
    if(verbose) cout<<"Event yield: "<<n_P_P<<" +/- "<<sta<<"(Sta.), "<<sys4<<"\%(Sys.4)"<<endl;
    else cout<<", "<<sys4<<"\%(Sys.4)"<<endl;
    return;
}