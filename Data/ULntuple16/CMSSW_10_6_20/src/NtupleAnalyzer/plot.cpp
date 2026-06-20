// Author: Shiyang CHEN
// Description: Draw weighted distributions from Weight*.root.
#include <iostream>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TSystem.h"
using namespace std;
#define PI 3.14159265359

// I/O and behavior switches.
// Change inputFileName to WeightDPS.root or WeightData.root when needed.
const string inputFileName = "WeightSPSstar.root";
const string outPlotDir = "plots/SPSstar/";
const bool drawWeightComparison = false;  // false: weighted only; true: unweighted vs weighted.

void drawWeightedOnly(TH1D &weighted, const string &outName) {
    TCanvas canvas(("c_" + outName).c_str(), outName.c_str(), 1000, 800);
    weighted.SetStats(0);
    weighted.Draw("hist");
    canvas.SaveAs((outPlotDir + outName + ".png").c_str());
}

void drawComparison(TH1D &unweighted, TH1D &weighted, const string &outName) {
    TCanvas canvas(("c_" + outName).c_str(), outName.c_str(), 1000, 500);
    canvas.Divide(2, 1);
    canvas.cd(1);
    unweighted.SetStats(0);
    unweighted.Draw("hist");
    canvas.cd(2);
    weighted.SetStats(0);
    weighted.Draw("hist");
    canvas.SaveAs((outPlotDir + outName + ".png").c_str());
}

void savePlot(TH1D &unweighted, TH1D &weighted, const string &outName) {
    if(drawWeightComparison) drawComparison(unweighted, weighted, outName);
    else drawWeightedOnly(weighted, outName);
}

void plot() {
    gSystem->mkdir(outPlotDir.c_str(), kTRUE);

    TFile inFile(inputFileName.c_str(), "READ");
    TTree *inTree = (TTree *)inFile.Get("data");
    if(!inTree) {
        cout<<"Cannot find tree data in "<<inputFileName<<endl;
        return;
    }
    // Define variables
    Double_t Jpsi_mass1, Jpsi_ctau1, Jpsi_pt1, Jpsi_Z1;
    Double_t Jpsi_mass2, Jpsi_ctau2, Jpsi_pt2;
    Double_t evt_weight, evt_mass, delta_y, delta_phi, evt_Z, vtx_Z;
    inTree->SetBranchAddress("Jpsi_mass1", &Jpsi_mass1);
    inTree->SetBranchAddress("Jpsi_ctau1", &Jpsi_ctau1);
    inTree->SetBranchAddress("Jpsi_mass2", &Jpsi_mass2);
    inTree->SetBranchAddress("Jpsi_ctau2", &Jpsi_ctau2);
    inTree->SetBranchAddress("evt_mass", &evt_mass);
    inTree->SetBranchAddress("delta_y", &delta_y);
    inTree->SetBranchAddress("delta_phi", &delta_phi);
    inTree->SetBranchAddress("evt_weight", &evt_weight);
    inTree->SetBranchAddress("Jpsi_pt1", &Jpsi_pt1);
    inTree->SetBranchAddress("Jpsi_pt2", &Jpsi_pt2);
    inTree->SetBranchAddress("Jpsi_Z1", &Jpsi_Z1);
    inTree->SetBranchAddress("evt_Z", &evt_Z);
    inTree->SetBranchAddress("vtx_Z", &vtx_Z);
    // Define histograms
    TH1D h_Jpsi_mass1("h_Jpsi_mass1", "Jpsi_mass_unweighted1", 50, 2.95, 3.25), wh_Jpsi_mass1("wh_Jpsi_mass1", "Jpsi_mass_weighted1", 50, 2.95, 3.25);
    TH1D h_Jpsi_mass2("h_Jpsi_mass2", "Jpsi_mass_unweighted2", 50, 2.95, 3.25), wh_Jpsi_mass2("wh_Jpsi_mass2", "Jpsi_mass_weighted2", 50, 2.95, 3.25);
    TH1D h_Jpsi_ctau1("h_Jpsi_ctau1", "Jpsi_ctau_unweighted1", 50, -0.02, 0.15), wh_Jpsi_ctau1("wh_Jpsi_ctau1", "Jpsi_ctau_weighted1", 50, -0.02, 0.15);
    TH1D h_Jpsi_ctau2("h_Jpsi_ctau2", "Jpsi_ctau_unweighted2", 50, -0.02, 0.15), wh_Jpsi_ctau2("wh_Jpsi_ctau2", "Jpsi_ctau_weighted2", 50, -0.02, 0.15);
    wh_Jpsi_mass1.Sumw2();
    wh_Jpsi_mass2.Sumw2();
    wh_Jpsi_ctau1.Sumw2();
    wh_Jpsi_ctau2.Sumw2();
    TH1D h_evt_mass("h_evt_mass", "evt_mass_unweighted", 100, 7.5, 107.5), wh_evt_mass("wh_evt_mass", "evt_mass_weighted", 100, 7.5, 107.5);
    TH1D h_delta_y("h_delta_y", "delta_y_unweighted", 80, 0, 4.0), wh_delta_y("wh_delta_y", "delta_y_weighted", 80, 0, 4.0);
    TH1D h_delta_phi("h_delta_phi", "delta_phi_unweighted", 40, 0, PI), wh_delta_phi("wh_delta_phi", "delta_phi_weighted", 40, 0, PI);//40 DPS,80 SPS
    wh_evt_mass.Sumw2();
    wh_delta_y.Sumw2();
    wh_delta_phi.Sumw2();
    TH1D h_Jpsi_pt1("h_Jpsi_pt1", "Jpsi_pt_unweighted1", 50, 10, 40), wh_Jpsi_pt1("wh_Jpsi_pt1", "Jpsi_pt_weighted1", 50, 10, 40);
    TH1D h_Jpsi_pt2("h_Jpsi_pt2", "Jpsi_pt_unweighted2", 50, 10, 40), wh_Jpsi_pt2("wh_Jpsi_pt2", "Jpsi_pt_weighted2", 50, 10, 40);
    TH1D h_Jpsi_Z1("h_Jpsi_Z1", "Jpsi_Z_unweighted1", 50, -10, 10), wh_Jpsi_Z1("wh_Jpsi_Z1", "Jpsi_Z_weighted1", 50, -10, 10);
    TH1D h_evt_Z("h_evt_Z", "evt_Z_unweighted", 50, -10, 10), wh_evt_Z("wh_evt_Z", "evt_Z_weighted", 50, -10, 10);
    TH1D h_vtx_Z("h_vtx_Z", "vtx_Z_unweighted", 50, -10, 10), wh_vtx_Z("wh_vtx_Z", "vtx_Z_weighted", 50, -10, 10);
    wh_Jpsi_pt1.Sumw2();
    wh_Jpsi_pt2.Sumw2();
    wh_Jpsi_Z1.Sumw2();
    wh_evt_Z.Sumw2();
    wh_vtx_Z.Sumw2();
    int nEvent = inTree->GetEntries();
    for(int i = 0; i < nEvent; i++) {
        inTree->GetEntry(i);
        h_Jpsi_mass1.Fill(Jpsi_mass1);
        wh_Jpsi_mass1.Fill(Jpsi_mass1, evt_weight);
        h_Jpsi_mass2.Fill(Jpsi_mass2);
        wh_Jpsi_mass2.Fill(Jpsi_mass2, evt_weight);
        h_Jpsi_ctau1.Fill(Jpsi_ctau1);
        wh_Jpsi_ctau1.Fill(Jpsi_ctau1, evt_weight);
        h_Jpsi_ctau2.Fill(Jpsi_ctau2);
        wh_Jpsi_ctau2.Fill(Jpsi_ctau2, evt_weight);
        h_delta_y.Fill(delta_y);
        wh_delta_y.Fill(delta_y, evt_weight);
        h_delta_phi.Fill(delta_phi);
        wh_delta_phi.Fill(delta_phi, evt_weight);
        h_evt_mass.Fill(evt_mass);
        wh_evt_mass.Fill(evt_mass, evt_weight);
        h_Jpsi_pt1.Fill(Jpsi_pt1);
        wh_Jpsi_pt1.Fill(Jpsi_pt1, evt_weight);
        h_Jpsi_pt2.Fill(Jpsi_pt2);
        wh_Jpsi_pt2.Fill(Jpsi_pt2, evt_weight);
        h_Jpsi_Z1.Fill(Jpsi_Z1);
        wh_Jpsi_Z1.Fill(Jpsi_Z1, evt_weight);
        h_evt_Z.Fill(evt_Z);
        wh_evt_Z.Fill(evt_Z, evt_weight);
        h_vtx_Z.Fill(vtx_Z);
        wh_vtx_Z.Fill(vtx_Z, evt_weight);
    }
    savePlot(h_Jpsi_mass1, wh_Jpsi_mass1, "JpsiMass1");
    savePlot(h_Jpsi_mass2, wh_Jpsi_mass2, "JpsiMass2");
    savePlot(h_Jpsi_ctau1, wh_Jpsi_ctau1, "JpsiCtau1");
    savePlot(h_Jpsi_ctau2, wh_Jpsi_ctau2, "JpsiCtau2");
    savePlot(h_evt_mass, wh_evt_mass, "evtMass");
    savePlot(h_delta_y, wh_delta_y, "deltaY");
    savePlot(h_delta_phi, wh_delta_phi, "deltaPhi");
    savePlot(h_Jpsi_pt1, wh_Jpsi_pt1, "JpsiPt1");
    savePlot(h_Jpsi_pt2, wh_Jpsi_pt2, "JpsiPt2");
    savePlot(h_Jpsi_Z1, wh_Jpsi_Z1, "JpsiZ1");
    savePlot(h_evt_Z, wh_evt_Z, "evtZ");
    savePlot(h_vtx_Z, wh_vtx_Z, "vtxZ");

    cout<<"Input file: "<<inputFileName<<endl;
    cout<<"Total events: "<<nEvent<<endl;
    cout<<"drawWeightComparison = "<<drawWeightComparison<<endl;
    return;
}
