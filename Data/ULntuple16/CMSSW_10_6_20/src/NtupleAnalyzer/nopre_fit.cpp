#include "TCanvas.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TLatex.h"
#include "RooAbsPdf.h"
#include "RooExtendPdf.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooGExpModel.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooProdPdf.h"
#include "RooFitResult.h"

using namespace std;
using namespace RooFit;

void nopre_fit() {
    Int_t N = 60000, BinNum = 100;
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
    // RooRealVar Jpsi_alpha1("Jpsi_alpha1", "Jpsi_alpha1", 1, 0.1, 5);
    // RooRealVar Jpsi_nx1("Jpsi_nx1", "Jpsi_nx1", 1, 0, 50);
    RooRealVar Jpsi_devia2("Jpsi_devia2", "Jpsi_devia2", 0.05, 0.02, 0.08);
    RooRealVar Jpsi_alpha2("Jpsi_alpha2", "Jpsi_alpha2", 1.5, 0.1, 3.5);
    RooRealVar Jpsi_nx2("Jpsi_nx2", "Jpsi_nx2", 1, 0, 50);
    RooRealVar Jpsi_ratio("Jpsi_ratio", "Jpsi_ratio", 0.6, 0, 1);
    // RooCBShape Jpsi_crysBall1_1("Jpsi_crysBall1_1", "Jpsi_crysBall1_1", Jpsi_mass1, Jpsi_mean, Jpsi_devia1, Jpsi_alpha1, Jpsi_nx1);
    RooGaussian Jpsi_crysBall1_1("Jpsi_crysBall1_1", "Jpsi_crysBall1_1", Jpsi_mass1, Jpsi_mean, Jpsi_devia1);
    RooCBShape Jpsi_crysBall2_1("Jpsi_crysBall2_1", "Jpsi_crysBall2_1", Jpsi_mass1, Jpsi_mean, Jpsi_devia2, Jpsi_alpha2, Jpsi_nx2);
    // RooCBShape Jpsi_crysBall1_2("Jpsi_crysBall1_2", "Jpsi_crysBall1_2", Jpsi_mass2, Jpsi_mean, Jpsi_devia1, Jpsi_alpha1, Jpsi_nx1);
    RooGaussian Jpsi_crysBall1_2("Jpsi_crysBall1_2", "Jpsi_crysBall1_2", Jpsi_mass2, Jpsi_mean, Jpsi_devia1);
    RooCBShape Jpsi_crysBall2_2("Jpsi_crysBall2_2", "Jpsi_crysBall2_2", Jpsi_mass2, Jpsi_mean, Jpsi_devia2, Jpsi_alpha2, Jpsi_nx2);
    RooAddPdf JpsiMassSig1("JpsiMassSig1", "JpsiMassSig1", RooArgList(Jpsi_crysBall1_1, Jpsi_crysBall2_1), Jpsi_ratio);
    RooAddPdf JpsiMassSig2("JpsiMassSig2", "JpsiMassSig2", RooArgList(Jpsi_crysBall1_2, Jpsi_crysBall2_2), Jpsi_ratio);
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
    // RooRealVar Jpsi_mu3("Jpsi_mu3", "Jpsi_mu3", 0, -0.005, 0.005);
    // RooRealVar Jpsi_sigma6("Jpsi_sigma6", "Jpsi_sigma6", 0.003, 0.001, 0.01);//0.02
    RooRealVar Jpsi_sigma7("Jpsi_sigma7", "Jpsi_sigma7", 0.003, 0.001, 0.02);//0.01
    RooRealVar Jpsi_coef3("Jpsi_coef3", "Jpsi_coef3", 0.06, 0.01, 0.2);
    // RooRealVar Jpsi_prop3("Jpsi_prop3", "Jpsi_prop3", 0.5, 0, 1);
    // RooRealVar Jpsi_alpha3("Jpsi_alpha3", "Jpsi_alpha3", 1.5, 0.1, 10);
    // RooRealVar Jpsi_nx3("Jpsi_nx3", "Jpsi_nx3", 3, 0.5, 20);
    // RooCBShape Jpsi_gauss4("Jpsi_gauss4", "Jpsi_gauss4", Jpsi_ctau1, Jpsi_mu3, Jpsi_sigma6, Jpsi_alpha3, Jpsi_nx3);
    // RooGaussian Jpsi_gauss4("Jpsi_gauss4", "Jpsi_gauss4", Jpsi_ctau1, Jpsi_mu3, Jpsi_sigma6);
    // RooGExpModel Jpsi_expGau3("Jpsi_expGau3", "Jpsi_expGau3", Jpsi_ctau1, Jpsi_sigma7, Jpsi_coef3, false, RooGExpModel::Type::Flipped);
    // RooAddPdf JpsiCtauCombBkg("JpsiCtauCombBkg", "JpsiCtauCombBkg", RooArgList(Jpsi_gauss4, Jpsi_expGau3), Jpsi_prop3);
    RooGExpModel JpsiCtauCombBkg1("JpsiCtauCombBkg1", "JJpsiCtauCombBkg1", Jpsi_ctau1, Jpsi_sigma7, Jpsi_coef3, false, RooGExpModel::Type::Flipped);
    RooGExpModel JpsiCtauCombBkg2("JpsiCtauCombBkg2", "JJpsiCtauCombBkg2", Jpsi_ctau2, Jpsi_sigma7, Jpsi_coef3, false, RooGExpModel::Type::Flipped);
    
    // Form 4-dim p.d.f.
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
    RooRealVar n_P_P("n_P_P", "n_P_P", 1e3, 1e1, N);//1e2
    RooRealVar n_P_NP("n_P_NP", "n_P_NP", 1e3, 1, N);//1e1
    // RooRealVar n_NP_P("n_NP_P", "n_NP_P", 1e3, 0, N);//1e1
    RooRealVar n_NP_NP("n_NP_NP", "n_NP_NP", 1e3, 1, N);
    RooRealVar n_Sig_Comb("n_Sig_Comb", "n_Sig_Comb", 1e4, 1, N);
    // RooRealVar n_Comb_Sig("n_Comb_Sig", "n_Comb_Sig", 1e2, 0, N);
    RooRealVar n_Comb_Comb("n_Comb_Comb", "n_Comb_Comb", 1e2, 1, N);
    RooAddPdf pdf_all("pdf_all", "pdf_all",
        RooArgList(pdf_P_P, pdf_P_NP, pdf_NP_P, pdf_NP_NP, pdf_Sig_Comb, pdf_Comb_Sig, pdf_Comb_Comb),//
        RooArgList(n_P_P, n_P_NP, n_P_NP, n_NP_NP, n_Sig_Comb, n_Sig_Comb, n_Comb_Comb)//
    );
    RooFitResult *res;
    while(true) {
        res = pdf_all.fitTo(*data->reduce(RooArgSet(Jpsi_mass1, Jpsi_mass2, Jpsi_ctau1, Jpsi_ctau2)), Save());
        if(!res->status() && res->edm()<0.01) break;//
    }
    // res = pdf_all.fitTo(*data->reduce(RooArgSet(Jpsi_mass1, Jpsi_mass2, Jpsi_ctau1, Jpsi_ctau2)), Save());
    // Draw data point and p.d.f. curve
    TCanvas *canvas7 = new TCanvas("canvas7", "canvas7", 1500, 1000);
    RooPlot *frame7 = Jpsi_mass1.frame(RooFit::Title("Jpsi-1 Mass 4D"), Bins(BinNum));
    data->plotOn(frame7, DataError(RooAbsData::SumW2), Name("Data"));
    pdf_all.plotOn(frame7, LineColor(kBlack), LineWidth(2), Name("All"));
    pdf_all.plotOn(frame7, Components(pdf_P_P), LineColor(kBlue), LineStyle(kSolid), LineWidth(2), Name("P_P"));
    pdf_all.plotOn(frame7, Components(pdf_P_NP), LineColor(kBlue), LineStyle(kDashed), LineWidth(2), Name("P_NP"));
    pdf_all.plotOn(frame7, Components(pdf_NP_P), LineColor(kBlue), LineStyle(kDotted), LineWidth(2), Name("NP_P"));
    pdf_all.plotOn(frame7, Components(pdf_NP_NP), LineColor(kBlue), LineStyle(kDashDotted), LineWidth(2), Name("NP_NP"));
    pdf_all.plotOn(frame7, Components(pdf_Sig_Comb), LineColor(kRed), LineWidth(1), Name("Sig_Comb"));
    pdf_all.plotOn(frame7, Components(pdf_Comb_Sig), LineColor(kGreen), LineWidth(1), Name("Comb_Sig"));
    pdf_all.plotOn(frame7, Components(pdf_Comb_Comb), LineColor(kMagenta), LineWidth(1), Name("Comb_Comb"));
    TLegend *legend7 = new TLegend(.65, .60, .85, .85);
    legend7->AddEntry(frame7->findObject("Data"), "RunII 2018", "L");
    legend7->AddEntry(frame7->findObject("All"), "Total p.d.f.", "L");
    legend7->AddEntry(frame7->findObject("P_P"), "prompt, prompt", "L");
    legend7->AddEntry(frame7->findObject("P_NP"), "prompt, non-prompt", "L");
    legend7->AddEntry(frame7->findObject("NP_P"), "non-prompt, prompt", "L");
    legend7->AddEntry(frame7->findObject("NP_NP"), "non-prompt, non-prompt", "L");
    legend7->AddEntry(frame7->findObject("Sig_Comb"), "J/#psi, #mu^{+}#mu^{-}", "L");
    legend7->AddEntry(frame7->findObject("Comb_Sig"), "#mu^{+}#mu^{-}, J/#psi", "L");
    legend7->AddEntry(frame7->findObject("Comb_Comb"), "#mu^{+}#mu^{-}, #mu^{+}#mu^{-}", "L");
    frame7->Draw();
    legend7->DrawClone();
    // canvas7->SaveAs("fits/4D_JpsiMass.pdf");
    canvas7->SaveAs("fits/4D_JpsiMass1.png");
    TCanvas *canvas8 = new TCanvas("canvas8", "canvas8", 1500, 1000);
    RooPlot *frame8 = Jpsi_mass2.frame(RooFit::Title("Jpsi-2 Mass 4D"), Bins(BinNum));
    data->plotOn(frame8, DataError(RooAbsData::SumW2), Name("Data"));
    pdf_all.plotOn(frame8, LineColor(kBlack), LineWidth(2), Name("All"));
    pdf_all.plotOn(frame8, Components(pdf_P_P), LineColor(kBlue), LineStyle(kSolid), LineWidth(2), Name("P_P"));
    pdf_all.plotOn(frame8, Components(pdf_P_NP), LineColor(kBlue), LineStyle(kDashed), LineWidth(2), Name("P_NP"));
    pdf_all.plotOn(frame8, Components(pdf_NP_P), LineColor(kBlue), LineStyle(kDotted), LineWidth(2), Name("NP_P"));
    pdf_all.plotOn(frame8, Components(pdf_NP_NP), LineColor(kBlue), LineStyle(kDashDotted), LineWidth(2), Name("NP_NP"));
    pdf_all.plotOn(frame8, Components(pdf_Sig_Comb), LineColor(kRed), LineWidth(1), Name("Sig_Comb"));
    pdf_all.plotOn(frame8, Components(pdf_Comb_Sig), LineColor(kGreen), LineWidth(1), Name("Comb_Sig"));
    pdf_all.plotOn(frame8, Components(pdf_Comb_Comb), LineColor(kMagenta), LineWidth(1), Name("Comb_Comb"));
    TLegend *legend8 = new TLegend(.65, .60, .85, .85);
    legend8->AddEntry(frame8->findObject("Data"), "RunII 2018", "L");
    legend8->AddEntry(frame8->findObject("All"), "Total p.d.f.", "L");
    legend8->AddEntry(frame8->findObject("P_P"), "prompt, prompt", "L");
    legend8->AddEntry(frame8->findObject("P_NP"), "prompt, non-prompt", "L");
    legend8->AddEntry(frame8->findObject("NP_P"), "non-prompt, prompt", "L");
    legend8->AddEntry(frame8->findObject("NP_NP"), "non-prompt, non-prompt", "L");
    legend8->AddEntry(frame8->findObject("Sig_Comb"), "J/#psi, #mu^{+}#mu^{-}", "L");
    legend8->AddEntry(frame8->findObject("Comb_Sig"), "#mu^{+}#mu^{-}, J/#psi", "L");
    legend8->AddEntry(frame8->findObject("Comb_Comb"), "#mu^{+}#mu^{-}, #mu^{+}#mu^{-}", "L");
    frame8->Draw();
    legend8->DrawClone();
    // canvas8->SaveAs("fits/4D_psi2SMass.pdf");
    canvas8->SaveAs("fits/4D_JpsiMass2.png");
    TCanvas *canvas9 = new TCanvas("canvas9", "canvas9", 1500, 1000);
    RooPlot *frame9 = Jpsi_ctau1.frame(RooFit::Title("Jpsi-1 Ctau 4D"), Bins(BinNum));
    data->plotOn(frame9, DataError(RooAbsData::SumW2), Name("Data"));
    pdf_all.plotOn(frame9, LineColor(kBlack), LineWidth(2), Name("All"));
    pdf_all.plotOn(frame9, Components(pdf_P_P), LineColor(kBlue), LineStyle(kSolid), LineWidth(2), Name("P_P"));
    pdf_all.plotOn(frame9, Components(pdf_P_NP), LineColor(kBlue), LineStyle(kDashed), LineWidth(2), Name("P_NP"));
    pdf_all.plotOn(frame9, Components(pdf_NP_P), LineColor(kBlue), LineStyle(kDotted), LineWidth(2), Name("NP_P"));
    pdf_all.plotOn(frame9, Components(pdf_NP_NP), LineColor(kBlue), LineStyle(kDashDotted), LineWidth(2), Name("NP_NP"));
    pdf_all.plotOn(frame9, Components(pdf_Sig_Comb), LineColor(kRed), LineWidth(1), Name("Sig_Comb"));
    pdf_all.plotOn(frame9, Components(pdf_Comb_Sig), LineColor(kGreen), LineWidth(1), Name("Comb_Sig"));
    pdf_all.plotOn(frame9, Components(pdf_Comb_Comb), LineColor(kMagenta), LineWidth(1), Name("Comb_Comb"));
    TLegend *legend9 = new TLegend(.65, .60, .85, .85);
    legend9->AddEntry(frame9->findObject("Data"), "RunII 2018", "L");
    legend9->AddEntry(frame9->findObject("All"), "Total p.d.f.", "L");
    legend9->AddEntry(frame9->findObject("P_P"), "prompt, prompt", "L");
    legend9->AddEntry(frame9->findObject("P_NP"), "prompt, non-prompt", "L");
    legend9->AddEntry(frame9->findObject("NP_P"), "non-prompt, prompt", "L");
    legend9->AddEntry(frame9->findObject("NP_NP"), "non-prompt, non-prompt", "L");
    legend9->AddEntry(frame9->findObject("Sig_Comb"), "J/#psi, #mu^{+}#mu^{-}", "L");
    legend9->AddEntry(frame9->findObject("Comb_Sig"), "#mu^{+}#mu^{-}, J/#psi", "L");
    legend9->AddEntry(frame9->findObject("Comb_Comb"), "#mu^{+}#mu^{-}, #mu^{+}#mu^{-}", "L");
    frame9->SetAxisRange(1, 5e3, "Y");
    gPad->SetLogy();
    frame9->Draw();
    legend9->DrawClone();
    // canvas9->SaveAs("fits/4D_JpsiCtau.pdf");
    canvas9->SaveAs("fits/4D_JpsiCtau1.png");
    TCanvas *canvas0 = new TCanvas("canvas0", "canvas0", 1500, 1000);
    RooPlot *frame0 = Jpsi_ctau2.frame(RooFit::Title("Jpsi-2 Ctau 4D"), Bins(BinNum));
    data->plotOn(frame0, DataError(RooAbsData::SumW2), Name("Data"));
    pdf_all.plotOn(frame0, LineColor(kBlack), LineWidth(2), Name("All"));
    pdf_all.plotOn(frame0, Components(pdf_P_P), LineColor(kBlue), LineStyle(kSolid), LineWidth(2), Name("P_P"));
    pdf_all.plotOn(frame0, Components(pdf_P_NP), LineColor(kBlue), LineStyle(kDashed), LineWidth(2), Name("P_NP"));
    pdf_all.plotOn(frame0, Components(pdf_NP_P), LineColor(kBlue), LineStyle(kDotted), LineWidth(2), Name("NP_P"));
    pdf_all.plotOn(frame0, Components(pdf_NP_NP), LineColor(kBlue), LineStyle(kDashDotted), LineWidth(2), Name("NP_NP"));
    pdf_all.plotOn(frame0, Components(pdf_Sig_Comb), LineColor(kRed), LineWidth(1), Name("Sig_Comb"));
    pdf_all.plotOn(frame0, Components(pdf_Comb_Sig), LineColor(kGreen), LineWidth(1), Name("Comb_Sig"));
    pdf_all.plotOn(frame0, Components(pdf_Comb_Comb), LineColor(kMagenta), LineWidth(1), Name("Comb_Comb"));
    TLegend *legend0 = new TLegend(.65, .60, .85, .85);
    legend0->AddEntry(frame0->findObject("Data"), "RunII 2018", "L");
    legend0->AddEntry(frame0->findObject("All"), "Total p.d.f.", "L");
    legend0->AddEntry(frame0->findObject("P_P"), "prompt, prompt", "L");
    legend0->AddEntry(frame0->findObject("P_NP"), "prompt, non-prompt", "L");
    legend0->AddEntry(frame0->findObject("NP_P"), "non-prompt, prompt", "L");
    legend0->AddEntry(frame0->findObject("NP_NP"), "non-prompt, non-prompt", "L");
    legend0->AddEntry(frame0->findObject("Sig_Comb"), "J/#psi, #mu^{+}#mu^{-}", "L");
    legend0->AddEntry(frame0->findObject("Comb_Sig"), "#mu^{+}#mu^{-}, J/#psi", "L");
    legend0->AddEntry(frame0->findObject("Comb_Comb"), "#mu^{+}#mu^{-}, #mu^{+}#mu^{-}", "L");
    frame0->SetAxisRange(1, 5e3, "Y");
    gPad->SetLogy();
    frame0->Draw();
    legend0->DrawClone();
    // canvas0->SaveAs("fits/4D_psi2SCtau.pdf");
    canvas0->SaveAs("fits/4D_JpsiCtau2.png");
    pdf_all.getVariables()->Print("v");
    cout<<"Status: "<<res->status()<<endl;
    return;
}