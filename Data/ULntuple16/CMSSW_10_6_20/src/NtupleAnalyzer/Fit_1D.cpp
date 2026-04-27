#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
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
#include "RooWorkspace.h"

using namespace std;
using namespace RooFit;

void Fit_1D(string var, double vmin, double vmax) {
    // Define variables
    RooRealVar Jpsi_mass1("Jpsi_mass1", "Jpsi_mass1", 2.95, 3.25);
    RooRealVar Jpsi_mass2("Jpsi_mass2", "Jpsi_mass2", 2.95, 3.25);
    RooRealVar Jpsi_ctau1("Jpsi_ctau1", "Jpsi_ctau1", -0.03, 0.16);
    RooRealVar Jpsi_ctau2("Jpsi_ctau2", "Jpsi_ctau2", -0.03, 0.16);
    RooRealVar evt_weight("evt_weight", "evt_weight", 0, 100);
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
    // TFile *BdecayFile = new TFile("WeightBdecaySmall.root", "READ");
    // TFile *MixFile = new TFile("WeightMix.root", "READ");
    // TTree *BdecayTree = (TTree*)BdecayFile->Get("data");
    // TTree *MixTree = (TTree*)MixFile->Get("data");
    // RooDataSet *Bdecay = new RooDataSet("Bdecay", "Bdecay", BdecayTree, variables, sel.c_str(), "evt_weight");
    // RooDataSet *Mix = new RooDataSet("Mix", "Mix", MixTree, variables, sel.c_str(), "evt_weight");
    // Read parameters from differential 4D fit
    TFile *f = new TFile("Model_4D_diff.root");
    RooWorkspace *wsp = (RooWorkspace *)f->Get("wsp");
    
    // Define p.d.f.
    //// J/psi Mass
    ////// Mixture of SPS/DPS=2:1
    RooRealVar Jpsi_m1("Jpsi_m1", "Jpsi_m1", wsp->var("Jpsi_mean")->getVal());
    RooRealVar Jpsi_d1("Jpsi_d1", "Jpsi_d1", wsp->var("Jpsi_devia1")->getVal());
    RooRealVar Jpsi_d2("Jpsi_d2", "Jpsi_d2", wsp->var("Jpsi_devia2")->getVal());
    RooRealVar Jpsi_a2("Jpsi_a2", "Jpsi_a2", wsp->var("Jpsi_alpha2")->getVal());
    RooRealVar Jpsi_n2("Jpsi_n2", "Jpsi_n2", wsp->var("Jpsi_nx2")->getVal());
    RooRealVar Jpsi_r1("Jpsi_r1", "Jpsi_r1", wsp->var("Jpsi_ratio")->getVal());
    RooGaussian Jpsi_g1_1("Jpsi_g1_1", "Jpsi_g1_1", Jpsi_mass1, Jpsi_m1, Jpsi_d1);
    RooCBShape Jpsi_c1_1("Jpsi_c1_1", "Jpsi_c1_1", Jpsi_mass1, Jpsi_m1, Jpsi_d2, Jpsi_a2, Jpsi_n2);
    RooGaussian Jpsi_g1_2("Jpsi_g1_2", "Jpsi_g1_2", Jpsi_mass2, Jpsi_m1, Jpsi_d1);
    RooCBShape Jpsi_c1_2("Jpsi_c1_2", "Jpsi_c1_2", Jpsi_mass2, Jpsi_m1, Jpsi_d2, Jpsi_a2, Jpsi_n2);
    RooAddPdf JpsiMassMix1("JpsiMassMix1", "JpsiMassMix1", RooArgList(Jpsi_g1_1, Jpsi_c1_1), Jpsi_r1);
    RooAddPdf JpsiMassMix2("JpsiMassMix2", "JpsiMassMix2", RooArgList(Jpsi_g1_2, Jpsi_c1_2), Jpsi_r1);
    ////// Bdecay
    RooRealVar Jpsi_m2("Jpsi_m2", "Jpsi_m2", wsp->var("Jpsi_mean")->getVal());
    RooRealVar Jpsi_d3("Jpsi_d3", "Jpsi_d3", wsp->var("Jpsi_devia1")->getVal());
    RooRealVar Jpsi_d4("Jpsi_d4", "Jpsi_d4", wsp->var("Jpsi_devia2")->getVal());
    RooRealVar Jpsi_a4("Jpsi_a4", "Jpsi_a4", wsp->var("Jpsi_alpha2")->getVal());
    RooRealVar Jpsi_n4("Jpsi_n4", "Jpsi_n4", wsp->var("Jpsi_nx2")->getVal());
    RooRealVar Jpsi_r2("Jpsi_r2", "Jpsi_r2", wsp->var("Jpsi_ratio")->getVal());
    RooGaussian Jpsi_g2_1("Jpsi_g2_1", "Jpsi_g2_1", Jpsi_mass1, Jpsi_m2, Jpsi_d3);
    RooCBShape Jpsi_c2_1("Jpsi_c2_1", "Jpsi_c2_1", Jpsi_mass1, Jpsi_m2, Jpsi_d4, Jpsi_a4, Jpsi_n4);
    RooGaussian Jpsi_g2_2("Jpsi_g2_2", "Jpsi_g2_2", Jpsi_mass2, Jpsi_m2, Jpsi_d3);
    RooCBShape Jpsi_c2_2("Jpsi_c2_2", "Jpsi_c2_2", Jpsi_mass2, Jpsi_m2, Jpsi_d4, Jpsi_a4, Jpsi_n4);
    RooAddPdf JpsiMassB1("JpsiMassB1", "JpsiMassB1", RooArgList(Jpsi_g2_1, Jpsi_c2_1), Jpsi_r2);
    RooAddPdf JpsiMassB2("JpsiMassB2", "JpsiMassB2", RooArgList(Jpsi_g2_2, Jpsi_c2_2), Jpsi_r2);
    ////// Combinatorial signal
    RooRealVar Jpsi_r3("Jpsi_r3", "Jpsi_r3", wsp->var("Jpsi_prop2")->getVal());
    RooAddPdf JpsiMassCSig1("JpsiMassCSig1", "JpsiMassCSig1", RooArgList(JpsiMassMix1, JpsiMassB1), Jpsi_r3);
    RooAddPdf JpsiMassCSig2("JpsiMassCSig2", "JpsiMassCSig2", RooArgList(JpsiMassMix2, JpsiMassB2), Jpsi_r3);
    ////// Combinatorial background
    RooChebychev JpsiMassCBkg1("JpsiMassCBkg1", "JpsiMassCBkg1", Jpsi_mass1, RooArgList());
    RooChebychev JpsiMassCBkg2("JpsiMassCBkg2", "JpsiMassCBkg2", Jpsi_mass2, RooArgList());

    //// J/psi Ctau
    ////// Mixture of SPS/DPS=2:1
    RooRealVar Jpsi_u1("Jpsi_u1", "Jpsi_u1", wsp->var("Jpsi_mu1")->getVal());
    RooRealVar Jpsi_s1("Jpsi_s1", "Jpsi_s1", wsp->var("Jpsi_sigma1")->getVal());
    RooRealVar Jpsi_s2("Jpsi_s2", "Jpsi_s2", wsp->var("Jpsi_sigma2")->getVal());
    RooRealVar Jpsi_p1("Jpsi_p1", "Jpsi_p1", wsp->var("Jpsi_prop1")->getVal());
    RooGaussian Jpsi_g3_1("Jpsi_g3_1", "Jpsi_g3_1", Jpsi_ctau1, Jpsi_u1, Jpsi_s1);
    RooGaussian Jpsi_g4_1("Jpsi_g4_1", "Jpsi_g4_1", Jpsi_ctau1, Jpsi_u1, Jpsi_s2);
    RooGaussian Jpsi_g3_2("Jpsi_g3_2", "Jpsi_g3_2", Jpsi_ctau2, Jpsi_u1, Jpsi_s1);
    RooGaussian Jpsi_g4_2("Jpsi_g4_2", "Jpsi_g4_2", Jpsi_ctau2, Jpsi_u1, Jpsi_s2);
    RooAddPdf JpsiCtauMix1("JpsiCtauMix1", "JpsiCtauMix1", RooArgList(Jpsi_g3_1, Jpsi_g4_1), Jpsi_p1);
    RooAddPdf JpsiCtauMix2("JpsiCtauMix2", "JpsiCtauMix2", RooArgList(Jpsi_g3_2, Jpsi_g4_2), Jpsi_p1);
    ////// Bdecay
    RooRealVar Jpsi_s3("Jpsi_s3", "Jpsi_s3", wsp->var("Jpsi_sigma3")->getVal());
    // RooRealVar Jpsi_s4("Jpsi_s4", "Jpsi_s4", 0.001, 0.002, 0.008);
    RooRealVar Jpsi_e1("Jpsi_e1", "Jpsi_e1", wsp->var("Jpsi_coef1")->getVal());
    // RooRealVar Jpsi_e2("Jpsi_e2", "Jpsi_e2", 0.03, 0.01, 0.1);
    // RooRealVar Jpsi_p2("Jpsi_p2", "Jpsi_p2", 0.5, 0, 1);
    // RooGExpModel Jpsi_ge1("Jpsi_ge1", "Jpsi_ge1", Jpsi_ctau, Jpsi_s3, Jpsi_e1, false, RooGExpModel::Type::Flipped);
    // RooGExpModel Jpsi_ge2("Jpsi_ge2", "Jpsi_ge2", Jpsi_ctau, Jpsi_s4, Jpsi_e2, false, RooGExpModel::Type::Flipped);
    // RooAddPdf JpsiCtauB("JpsiCtauB", "JpsiCtauB", RooArgList(Jpsi_ge1, Jpsi_ge2), Jpsi_p2);
    RooGExpModel JpsiCtauB1("JpsiCtauB1", "JpsiCtauB1", Jpsi_ctau1, Jpsi_s3, Jpsi_e1, false, RooGExpModel::Type::Flipped);
    RooGExpModel JpsiCtauB2("JpsiCtauB2", "JpsiCtauB2", Jpsi_ctau2, Jpsi_s3, Jpsi_e1, false, RooGExpModel::Type::Flipped);
    ////// Combinatorial signal
    RooRealVar Jpsi_p3("Jpsi_p3", "Jpsi_p3", wsp->var("Jpsi_prop2")->getVal());
    RooAddPdf JpsiCtauCSig1("JpsiCtauCSig1", "JpsiCtauCSig1", RooArgList(JpsiCtauMix1, JpsiCtauB1), Jpsi_p3);
    RooAddPdf JpsiCtauCSig2("JpsiCtauCSig2", "JpsiCtauCSig2", RooArgList(JpsiCtauMix2, JpsiCtauB2), Jpsi_p3);
    ////// Combinatorial background
    RooRealVar Jpsi_s7("Jpsi_s7", "Jpsi_s7", wsp->var("Jpsi_sigma7")->getVal());
    RooRealVar Jpsi_e3("Jpsi_e3", "Jpsi_e3", wsp->var("Jpsi_coef3")->getVal());
    RooGExpModel JpsiCtauCBkg1("JpsiCtauCBkg1", "JpsiCtauCBkg1", Jpsi_ctau1, Jpsi_s7, Jpsi_e3, false, RooGExpModel::Type::Flipped);
    RooGExpModel JpsiCtauCBkg2("JpsiCtauCBkg2", "JpsiCtauCBkg2", Jpsi_ctau2, Jpsi_s7, Jpsi_e3, false, RooGExpModel::Type::Flipped);
    
    // Form 4-dim p.d.f.
    // double scale = 5056 / wsp->var("n_P_P")->getVal();
    // RooRealVar vScale("vScale", "vScale", scale);
    RooProdPdf model_P_P("model_P_P", "model_P_P", RooArgList(JpsiMassMix1, JpsiMassMix2, JpsiCtauMix1, JpsiCtauMix2));
    RooProdPdf model_P_NP("model_P_NP", "model_P_NP", RooArgList(JpsiMassMix1, JpsiMassB2, JpsiCtauMix1, JpsiCtauB2));
    RooProdPdf model_NP_P("model_NP_P", "model_NP_P", RooArgList(JpsiMassB1, JpsiMassMix2, JpsiCtauB1, JpsiCtauMix2));
    RooProdPdf model_NP_NP("model_NP_NP", "model_NP_NP", RooArgList(JpsiMassB1, JpsiMassB2, JpsiCtauB1, JpsiCtauB2));
    RooProdPdf model_Sig_Comb("model_Sig_Comb", "model_Sig_Comb", RooArgList(JpsiMassCSig1, JpsiMassCBkg2, JpsiCtauCSig1, JpsiCtauCBkg2));
    RooProdPdf model_Comb_Sig("model_Comb_Sig", "model_Comb_Sig", RooArgList(JpsiMassCBkg1, JpsiMassCSig2, JpsiCtauCBkg1, JpsiCtauCSig2));
    RooProdPdf model_Comb_Comb("model_Comb_Comb", "model_Comb_Comb", RooArgList(JpsiMassCBkg1, JpsiMassCBkg2, JpsiCtauCBkg1, JpsiCtauCBkg2));
    // RooRealVar y_P_P("y_P_P", "y_P_P", 5056);
    // RooRealVar y_P_NP("y_P_NP", "y_P_NP", wsp->var("n_P_NP")->getVal() * scale);
    // RooRealVar y_NP_NP("y_NP_NP", "y_NP_NP", wsp->var("n_NP_NP")->getVal() * scale);
    // RooRealVar y_Sig_Comb("y_Sig_Comb", "y_Sig_Comb", wsp->var("n_Sig_Comb")->getVal() * scale);
    // RooRealVar y_Comb_Comb("y_Comb_Comb", "y_Comb_Comb", wsp->var("n_Comb_Comb")->getVal() * scale);
    RooRealVar y_P_P("y_P_P", "y_P_P", wsp->var("n_P_P")->getVal());
    RooRealVar y_P_NP("y_P_NP", "y_P_NP", wsp->var("n_P_NP")->getVal());
    RooRealVar y_NP_NP("y_NP_NP", "y_NP_NP", wsp->var("n_NP_NP")->getVal());
    RooRealVar y_Sig_Comb("y_Sig_Comb", "y_Sig_Comb", wsp->var("n_Sig_Comb")->getVal());
    RooRealVar y_Comb_Comb("y_Comb_Comb", "y_Comb_Comb", wsp->var("n_Comb_Comb")->getVal());
    RooAddPdf model_all("model_all", "model_all",
        RooArgList(model_P_P, model_P_NP, model_NP_P, model_NP_NP, model_Sig_Comb, model_Comb_Sig, model_Comb_Comb),
        RooArgList(y_P_P, y_P_NP, y_P_NP, y_NP_NP, y_Sig_Comb, y_Sig_Comb, y_Comb_Comb)
    );
    // Generate data
    Jpsi_m1.setConstant(kTRUE);
    Jpsi_d1.setConstant(kTRUE);
    Jpsi_d2.setConstant(kTRUE);
    Jpsi_a2.setConstant(kTRUE);
    Jpsi_n2.setConstant(kTRUE);
    Jpsi_r1.setConstant(kTRUE);

    Jpsi_m2.setConstant(kTRUE);
    Jpsi_d3.setConstant(kTRUE);
    Jpsi_d4.setConstant(kTRUE);
    Jpsi_a4.setConstant(kTRUE);
    Jpsi_n4.setConstant(kTRUE);
    Jpsi_r2.setConstant(kTRUE);

    Jpsi_r3.setConstant(kTRUE);

    Jpsi_u1.setConstant(kTRUE);
    Jpsi_s1.setConstant(kTRUE);
    Jpsi_s2.setConstant(kTRUE);
    Jpsi_p1.setConstant(kTRUE);

    Jpsi_s3.setConstant(kTRUE);
    Jpsi_e1.setConstant(kTRUE);

    Jpsi_p3.setConstant(kTRUE);

    Jpsi_s7.setConstant(kTRUE);
    Jpsi_e3.setConstant(kTRUE);

    RooWorkspace *wk = new RooWorkspace("wk", "wk");
    // wk->import(vScale);
    wk->import(model_all);
    wk->writeToFile("Data_4D.root");
    // Display fitting status
    // cout<<"Status of all 1D fits:"<<endl;
    // cout<<res1->status()<<' '<<res2->status()<<' '<<res3->status()<<' '<<res4->status()<<endl;
    // cout<<res5->status()<<' '<<res6->status()<<' '<<res7->status()<<' '<<res8->status()<<endl;
    return;
}