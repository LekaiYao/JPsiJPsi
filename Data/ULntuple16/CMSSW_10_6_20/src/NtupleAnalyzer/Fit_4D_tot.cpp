#include "Plot_4D.hpp"

void Fit_4D_tot() {
    // Define variables and read input file
    RooRealVar Jpsi_mass1("Jpsi_mass1", "Jpsi_mass1", 2.95, 3.25);
    RooRealVar Jpsi_mass2("Jpsi_mass2", "Jpsi_mass2", 2.95, 3.25);
    RooRealVar Jpsi_ctau1("Jpsi_ctau1", "Jpsi_ctau1", -0.03, 0.16);
    RooRealVar Jpsi_ctau2("Jpsi_ctau2", "Jpsi_ctau2", -0.03, 0.16);
    RooRealVar evt_weight("evt_weight", "evt_weight", 0, 1000);
    RooArgSet variables;
    variables.add(Jpsi_mass1);
    variables.add(Jpsi_mass2);
    variables.add(Jpsi_ctau1);
    variables.add(Jpsi_ctau2);
    variables.add(evt_weight);
    TFile *dataFile = new TFile("WeightData.root", "READ");
    TTree *dataTree = (TTree*)dataFile->Get("data");
    RooDataSet *data = new RooDataSet("data", "data", dataTree, variables, "", "evt_weight");
    
    // Define J/psi Mass p.d.f.
    // Signal p.d.f.
    RooRealVar Jpsi_mean("Jpsi_mean", "Jpsi_mean", 3.0969, 3.05, 3.15);
    RooRealVar Jpsi_devia1("Jpsi_devia1", "Jpsi_devia1", 0.01, 0, 0.03);
    RooRealVar Jpsi_devia2("Jpsi_devia2", "Jpsi_devia2", 0.05, 0.02, 0.08);
    RooRealVar Jpsi_alpha2("Jpsi_alpha2", "Jpsi_alpha2", 1.5, 0.1, 3.5);
    RooRealVar Jpsi_nx2("Jpsi_nx2", "Jpsi_nx2", 1, 0, 500);
    RooRealVar Jpsi_ratio("Jpsi_ratio", "Jpsi_ratio", 0.6, 0, 1);
    RooGaussian Jpsi_gaussian_1("Jpsi_gaussian_1", "Jpsi_gaussian_1", Jpsi_mass1, Jpsi_mean, Jpsi_devia1);
    RooCBShape Jpsi_crysBall_1("Jpsi_crysBall_1", "Jpsi_crysBall_1", Jpsi_mass1, Jpsi_mean, Jpsi_devia2, Jpsi_alpha2, Jpsi_nx2);
    RooGaussian Jpsi_gaussian_2("Jpsi_gaussian_2", "Jpsi_gaussian_2", Jpsi_mass2, Jpsi_mean, Jpsi_devia1);
    RooCBShape Jpsi_crysBall_2("Jpsi_crysBall_2", "Jpsi_crysBall_2", Jpsi_mass2, Jpsi_mean, Jpsi_devia2, Jpsi_alpha2, Jpsi_nx2);
    RooAddPdf JpsiMassSig1("JpsiMassSig1", "JpsiMassSig1", RooArgList(Jpsi_gaussian_1, Jpsi_crysBall_1), Jpsi_ratio);
    RooAddPdf JpsiMassSig2("JpsiMassSig2", "JpsiMassSig2", RooArgList(Jpsi_gaussian_2, Jpsi_crysBall_2), Jpsi_ratio);
    // Background p.d.f.
    RooChebychev JpsiMassComb1("JpsiMassComb1", "JpsiMassComb1", Jpsi_mass1, RooArgList());
    RooChebychev JpsiMassComb2("JpsiMassComb2", "JpsiMassComb2", Jpsi_mass2, RooArgList());

    // Define J/psi Ctau p.d.f.
    // Signal p.d.f.
    RooRealVar Jpsi_mu1("Jpsi_mu1", "Jpsi_mu1", 0, -0.005, 0.005);
    RooRealVar Jpsi_sigma1("Jpsi_sigma1", "Jpsi_sigma1", 0.001, 0, 0.003);
    RooRealVar Jpsi_sigma2("Jpsi_sigma2", "Jpsi_sigma2", 0.004, 0.002, 0.008);
    RooRealVar Jpsi_prop1("Jpsi_prop1", "Jpsi_prop1", 0.5, 0, 1);
    RooGaussian Jpsi_gauss1_1("Jpsi_gauss1_1", "Jpsi_gauss1_1", Jpsi_ctau1, Jpsi_mu1, Jpsi_sigma1);
    RooGaussian Jpsi_gauss2_1("Jpsi_gauss2_1", "Jpsi_gauss2_1", Jpsi_ctau1, Jpsi_mu1, Jpsi_sigma2);
    RooGaussian Jpsi_gauss1_2("Jpsi_gauss1_2", "Jpsi_gauss1_2", Jpsi_ctau2, Jpsi_mu1, Jpsi_sigma1);
    RooGaussian Jpsi_gauss2_2("Jpsi_gauss2_2", "Jpsi_gauss2_2", Jpsi_ctau2, Jpsi_mu1, Jpsi_sigma2);
    RooAddPdf JpsiCtauSig1("JpsiCtauSig1", "JpsiCtauSig1", RooArgList(Jpsi_gauss1_1, Jpsi_gauss2_1), Jpsi_prop1);
    RooAddPdf JpsiCtauSig2("JpsiCtauSig2", "JpsiCtauSig2", RooArgList(Jpsi_gauss1_2, Jpsi_gauss2_2), Jpsi_prop1);
    // Background p.d.f.
    RooRealVar Jpsi_sigma3("Jpsi_sigma3", "Jpsi_sigma3", 0.001, 0, 0.01);
    RooRealVar Jpsi_coef1("Jpsi_coef1", "Jpsi_coef1", 0.06, 0.01, 0.1);
    RooGExpModel JpsiCtauBkg1("JpsiCtauBkg1", "JpsiCtauBkg1", Jpsi_ctau1, Jpsi_sigma3, Jpsi_coef1, false, RooGExpModel::Type::Flipped);
    RooGExpModel JpsiCtauBkg2("JpsiCtauBkg2", "JpsiCtauBkg2", Jpsi_ctau2, Jpsi_sigma3, Jpsi_coef1, false, RooGExpModel::Type::Flipped);
    // Combinatorial signal p.d.f.
    RooRealVar Jpsi_prop2("Jpsi_prop2", "Jpsi_prop2", 0.5, 0, 1);
    RooAddPdf JpsiCtauCombSig1("JpsiCtauCombSig1", "JpsiCtauCombSig1", RooArgList(JpsiCtauSig1, JpsiCtauBkg1), Jpsi_prop2);
    RooAddPdf JpsiCtauCombSig2("JpsiCtauCombSig2", "JpsiCtauCombSig2", RooArgList(JpsiCtauSig2, JpsiCtauBkg2), Jpsi_prop2);
    // Combinatorial background p.d.f.
    RooRealVar Jpsi_sigma7("Jpsi_sigma7", "Jpsi_sigma7", 0.003, 0.001, 0.02);
    RooRealVar Jpsi_coef3("Jpsi_coef3", "Jpsi_coef3", 0.06, 0.01, 0.2);
    RooGExpModel JpsiCtauCombBkg1("JpsiCtauCombBkg1", "JJpsiCtauCombBkg1", Jpsi_ctau1, Jpsi_sigma7, Jpsi_coef3, false, RooGExpModel::Type::Flipped);
    RooGExpModel JpsiCtauCombBkg2("JpsiCtauCombBkg2", "JJpsiCtauCombBkg2", Jpsi_ctau2, Jpsi_sigma7, Jpsi_coef3, false, RooGExpModel::Type::Flipped);
    
    // Form 4-dim p.d.f.
    Int_t N = 60000;
    // J/psi1 Mass + J/psi2 Mass dimension
    RooProdPdf pdf_mass_SigSig("pdf_mass_SigSig", "pdf_mass_SigSig", JpsiMassSig1, JpsiMassSig2);
    RooProdPdf pdf_mass_SigComb("pdf_mass_SigComb", "pdf_mass_SigComb", JpsiMassSig1, JpsiMassComb2);
    RooProdPdf pdf_mass_CombSig("pdf_mass_CombSig", "pdf_mass_CombSig", JpsiMassComb1, JpsiMassSig2);
    RooProdPdf pdf_mass_CombComb("pdf_mass_CombComb", "pdf_mass_CombComb", JpsiMassComb1, JpsiMassComb2);
    // J/psi1 Ctau + J/psi2 Ctau dimension
    RooProdPdf pdf_ctau_PP("pdf_ctau_PP", "pdf_ctau_PP", JpsiCtauSig1, JpsiCtauSig2);
    RooProdPdf pdf_ctau_PNP("pdf_ctau_PNP", "pdf_ctau_PNP", JpsiCtauSig1, JpsiCtauBkg2);
    RooProdPdf pdf_ctau_NPP("pdf_ctau_NPP", "pdf_ctau_NPP", JpsiCtauBkg1, JpsiCtauSig2);
    RooProdPdf pdf_ctau_NPNP("pdf_ctau_NPNP", "pdf_ctau_NPNP", JpsiCtauBkg1, JpsiCtauBkg2);
    RooProdPdf pdf_ctau_SigBkg("pdf_ctau_SigBkg", "pdf_ctau_SigBkg", JpsiCtauCombSig1, JpsiCtauCombBkg2);
    RooProdPdf pdf_ctau_BkgSig("pdf_ctau_BkgSig", "pdf_ctau_BkgSig", JpsiCtauCombBkg1, JpsiCtauCombSig2);
    RooProdPdf pdf_ctau_BkgBkg("pdf_ctau_BkgBkg", "pdf_ctau_BkgBkg", JpsiCtauCombBkg1, JpsiCtauCombBkg2);
    // Combine Mass and Ctau dimensions
    RooProdPdf pdf_P_P("pdf_P_P", "pdf_P_P", pdf_mass_SigSig, pdf_ctau_PP);
    RooProdPdf pdf_P_NP("pdf_P_NP", "pdf_P_NP", pdf_mass_SigSig, pdf_ctau_PNP);
    RooProdPdf pdf_NP_P("pdf_NP_P", "pdf_NP_P", pdf_mass_SigSig, pdf_ctau_NPP);
    RooProdPdf pdf_NP_NP("pdf_NP_NP", "pdf_NP_NP", pdf_mass_SigSig, pdf_ctau_NPNP);
    RooProdPdf pdf_Sig_Comb("pdf_Sig_Comb", "pdf_Sig_Comb", pdf_mass_SigComb, pdf_ctau_SigBkg);
    RooProdPdf pdf_Comb_Sig("pdf_Comb_Sig", "pdf_Comb_Sig", pdf_mass_CombSig, pdf_ctau_BkgSig);
    RooProdPdf pdf_Comb_Comb("pdf_Comb_Comb", "pdf_Comb_Comb", pdf_mass_CombComb, pdf_ctau_BkgBkg);
    RooRealVar n_P_P("n_P_P", "n_P_P", 1e3, 1e1, N);
    RooRealVar n_P_NP("n_P_NP", "n_P_NP", 1e3, 1, N);
    RooRealVar n_NP_NP("n_NP_NP", "n_NP_NP", 1e3, 1, N);
    RooRealVar n_Sig_Comb("n_Sig_Comb", "n_Sig_Comb", 1e4, 1, N);
    RooRealVar n_Comb_Comb("n_Comb_Comb", "n_Comb_Comb", 1e2, 1, N);
    RooAddPdf pdf_all("pdf_all", "pdf_all",
        RooArgList(pdf_P_P, pdf_P_NP, pdf_NP_P, pdf_NP_NP, pdf_Sig_Comb, pdf_Comb_Sig, pdf_Comb_Comb),
        RooArgList(n_P_P, n_P_NP, n_P_NP, n_NP_NP, n_Sig_Comb, n_Sig_Comb, n_Comb_Comb)
    );
    RooFitResult *res;
    while(true) {
        res = pdf_all.fitTo(*data->reduce(RooArgSet(Jpsi_mass1, Jpsi_mass2, Jpsi_ctau1, Jpsi_ctau2)), Save());
        if(!res->status() && res->edm()<0.01) break;
    }

    // Draw data point and p.d.f. curve
    string prefix = "fig/";
    Plot_4D(data, pdf_all, pdf_P_P, pdf_P_NP, pdf_NP_P, pdf_NP_NP, pdf_Sig_Comb, pdf_Comb_Sig, pdf_Comb_Comb,
        prefix, Jpsi_mass1, Jpsi_mass2, Jpsi_ctau1, Jpsi_ctau2);
    pdf_all.getVariables()->Print("v");

    // Set parameters to constant save model to file
    Jpsi_mean.setConstant(kTRUE);
    Jpsi_devia1.setConstant(kTRUE);
    Jpsi_devia2.setConstant(kTRUE);
    Jpsi_alpha2.setConstant(kTRUE);
    Jpsi_nx2.setConstant(kTRUE);
    Jpsi_ratio.setConstant(kTRUE);

    Jpsi_mu1.setConstant(kTRUE);
    Jpsi_sigma1.setConstant(kTRUE);
    Jpsi_sigma2.setConstant(kTRUE);
    Jpsi_prop1.setConstant(kTRUE);
    Jpsi_sigma3.setConstant(kTRUE);
    Jpsi_coef1.setConstant(kTRUE);
    Jpsi_prop2.setConstant(kTRUE);
    Jpsi_sigma7.setConstant(kTRUE);
    Jpsi_coef3.setConstant(kTRUE);

    RooWorkspace *wsp = new RooWorkspace("wsp", "wsp");
    wsp->import(pdf_all);
    wsp->writeToFile("Model_4D_tot.root");
    // Log message to screen
    cout<<"Status: "<<res->status()<<endl;
    cout<<"Event yield: "<<n_P_P.getVal()<<" +/- "<<n_P_P.getError()<<endl;
    return;
}