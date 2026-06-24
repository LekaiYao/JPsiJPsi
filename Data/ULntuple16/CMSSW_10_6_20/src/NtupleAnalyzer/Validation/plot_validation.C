// Appendix C MC-vs-data validation plot (sideband-subtracted data vs SPS/DPS MC).
// Templated: @DATAVAR@ = data branch base, @MCVAR@ = MC branch, @TITLE@ = axis title,
// @NBINS@/@XMIN@/@XMAX@ = histogram range. Produced by run_validation.sh.
//
// SPS MC = MC_SPS_NLOstar.root (nominal NLO*, made by make_validation_mc.cpp).
// DPS MC and the sWeighted data file are carried over unchanged from the original
// validation (only the SPS generator switched Pythia8 -> NLO*).
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TColor.h"

void plot_@DATAVAR@() {
    TString varName  = "@DATAVAR@_";   // data branch
    TString mcName   = "@MCVAR@";      // MC branch (WeightedTree)
    TString massVar1 = "J1Mass_", massVar2 = "J2Mass_";
    TString ctauVar1 = "x1_", ctauVar2 = "x2_";

    TFile *f_data = TFile::Open("Sample/Data_MixWeighted_cut0_sw.root");
    TFile *f_mc   = TFile::Open("MC_SPS_NLOstar.root");
    TFile *f_mc2  = TFile::Open("Sample/MC_DPS.root");

    double sig_low = 3.049, sig_high = 3.148;
    double sb1_low = 2.90, sb1_high = 3.00;
    double sb2_low = 3.20, sb2_high = 3.30;
    double ctau_low = -0.005, ctau_high = 0.004;

    TTree *tree_data = (TTree*)f_data->Get("tree");
    TTree *tree_mc1  = (TTree*)f_mc->Get("WeightedTree");
    TTree *tree_mc2  = (TTree*)f_mc2->Get("WeightedTree");

    double var_min = @XMIN@, var_max = @XMAX@;
    int nBins = @NBINS@;

    TH1F *h_mc  = new TH1F("h_mc",  mcName,  nBins, var_min, var_max);
    TH1F *h_mc2 = new TH1F("h_mc2", mcName,  nBins, var_min, var_max);
    TH1F *h_sig = new TH1F("h_sig", "sig",   nBins, var_min, var_max);
    TH1F *h_sb  = new TH1F("h_sb",  "sb",    nBins, var_min, var_max);
    TH1F *h_sbs = new TH1F("h_sbs", "sbs",   nBins, var_min, var_max);

    // MC: reco-level shape (unweighted, normalized) of SPS and DPS
    tree_mc1->Draw(mcName + ">>h_mc",  "", "goff");
    tree_mc2->Draw(mcName + ">>h_mc2", "", "goff");
    // Data: sideband-subtracted signal, prompt (ctau-cut)
    tree_data->Draw(varName + ">>h_sig", Form("%s>%.3f && %s<%.3f && %s>%.3f && %s<%.3f && (%s>%.3f && %s<%.3f && %s>%.3f && %s<%.3f)",
                    massVar1.Data(), sig_low, massVar1.Data(), sig_high,
                    massVar2.Data(), sig_low, massVar2.Data(), sig_high,
                    ctauVar1.Data(), ctau_low, ctauVar1.Data(), ctau_high,
                    ctauVar2.Data(), ctau_low, ctauVar2.Data(), ctau_high), "goff");
    tree_data->Draw(varName + ">>h_sb", Form("((%s>%.3f && %s<%.3f)||(%s>%.3f && %s<%.3f))&&((%s>%.3f && %s<%.3f)||(%s>%.3f && %s<%.3f))&& (%s>%.3f && %s<%.3f && %s>%.3f && %s<%.3f)",
                    massVar1.Data(), sb1_low, massVar1.Data(), sb1_high,
                    massVar1.Data(), sb2_low, massVar1.Data(), sb2_high,
                    massVar2.Data(), sb1_low, massVar2.Data(), sb1_high,
                    massVar2.Data(), sb2_low, massVar2.Data(), sb2_high,
                    ctauVar1.Data(), ctau_low, ctauVar1.Data(), ctau_high,
                    ctauVar2.Data(), ctau_low, ctauVar2.Data(), ctau_high), "goff");
    double scale = (sig_high - sig_low) / (sb2_high - sb2_low + sb1_high - sb1_low);
    h_sb->Scale(scale);
    h_sbs->Add(h_sig);
    h_sbs->Add(h_sb, -1);

    if (h_mc->Integral() > 0)  h_mc->Scale(1.0 / h_mc->Integral());
    if (h_mc2->Integral() > 0) h_mc2->Scale(1.0 / h_mc2->Integral());
    if (h_sbs->Integral() > 0) h_sbs->Scale(1.0 / h_sbs->Integral());

    int azure = TColor::GetColor(51, 153, 255);
    h_mc->SetLineColor(azure);  h_mc->SetLineWidth(2);
    h_mc->SetFillColor(azure);  h_mc->SetFillStyle(3345);
    h_mc2->SetLineColor(kOrange + 7); h_mc2->SetLineWidth(2);
    h_mc2->SetFillColorAlpha(kOrange + 7, 0.5); h_mc2->SetFillStyle(1001);
    int yellow = TColor::GetColor(255, 235, 0);
    h_sbs->SetLineColor(yellow); h_sbs->SetLineWidth(2);
    h_sbs->SetFillColorAlpha(yellow, 0.3); h_sbs->SetFillStyle(1001);

    float max_val = std::max({h_mc->GetMaximum(), h_mc2->GetMaximum(), h_sbs->GetMaximum()});
    h_mc->SetMaximum(max_val * 1.25);

    TCanvas *c1 = new TCanvas("c1", "Comparison", 800, 600);
    gStyle->SetOptStat(0);
    h_mc->SetTitle("Validation of @TITLE@; @TITLE@ ; Normalized Entries");
    h_mc->Draw("hist");        h_mc->Draw("E1 same");
    h_mc2->Draw("hist same");  h_mc2->Draw("E1 same");
    h_sbs->Draw("hist same");  h_sbs->Draw("E1 same");

    TLegend *leg = new TLegend(0.65, 0.68, 0.88, 0.88);
    leg->AddEntry(h_sbs, "Data (SBS)", "lf");
    leg->AddEntry(h_mc,  "SPS MC (NLO*)", "lf");
    leg->AddEntry(h_mc2, "DPS MC", "lf");
    leg->Draw();

    c1->SaveAs("output_sbs/@DATAVAR@_compare_SBS_MC.pdf");
}
