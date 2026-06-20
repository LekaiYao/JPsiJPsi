// Author: Shiyang CHEN
// Description: for 2016 NLO Monte-Carlo
// Implementation: 
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
using namespace std;
#define PI 3.14159265359

class Process {
    private:
    void readTree(string& fileName) {
        // Handle root file
        TFile file(fileName.c_str(), "READ");
        TTree *inTree = (TTree *)file.Get("rootuple/oniaTree");
        if(!inTree) return;
        int nEvent = inTree->GetEntries();
        // Define input tree variables
        Double_t GEJpsi1_pt, GEJpsi1_eta, GEJpsi1_y, GEJpsi1_phi, GEJpsi1_mass;
        Double_t GEJpsi2_pt, GEJpsi2_eta, GEJpsi2_y, GEJpsi2_phi, GEJpsi2_mass;
        // Set input tree SetBranchAddress address
        inTree->SetBranchAddress("GEJpsi1_pt", &GEJpsi1_pt);
        inTree->SetBranchAddress("GEJpsi1_eta", &GEJpsi1_eta);
        inTree->SetBranchAddress("GEJpsi1_y", &GEJpsi1_y);
        inTree->SetBranchAddress("GEJpsi1_phi", &GEJpsi1_phi);
        inTree->SetBranchAddress("GEJpsi1_mass", &GEJpsi1_mass);

        inTree->SetBranchAddress("GEJpsi2_pt", &GEJpsi2_pt);
        inTree->SetBranchAddress("GEJpsi2_eta", &GEJpsi2_eta);
        inTree->SetBranchAddress("GEJpsi2_y", &GEJpsi2_y);
        inTree->SetBranchAddress("GEJpsi2_phi", &GEJpsi2_phi);
        inTree->SetBranchAddress("GEJpsi2_mass", &GEJpsi2_mass);
        // Loop on input tree entries
        cout<<fileName<<" has events: "<<nEvent<<endl;
        for(int i = 0; i < nEvent; i++) {
            cout<<"Processing No."<<i<<'\r';
            inTree->GetEntry(i);
            // Apply fiducial cut
            if(GEJpsi1_pt > 40 || GEJpsi1_pt < 10 || GEJpsi2_pt > 40 || GEJpsi2_pt < 10) continue;
            if(fabs(GEJpsi1_y) > 2 || fabs(GEJpsi2_y) > 2) continue;

            TLorentzVector JpsiLV1, JpsiLV2;
            JpsiLV1.SetPtEtaPhiM(GEJpsi1_pt, GEJpsi1_eta, GEJpsi1_phi, GEJpsi1_mass);
            JpsiLV2.SetPtEtaPhiM(GEJpsi2_pt, GEJpsi2_eta, GEJpsi2_phi, GEJpsi2_mass);
            if((JpsiLV1 + JpsiLV2).M() < 7.5) continue;
            evt_mass.push_back((JpsiLV1 + JpsiLV2).M());
            delta_phi.push_back(PI - fabs(fabs(GEJpsi1_phi - GEJpsi2_phi) - PI));
            delta_y.push_back(fabs(GEJpsi1_y - GEJpsi2_y));
            evt_pt.push_back((JpsiLV1 + JpsiLV2).Pt());
            evt_y.push_back(fabs((JpsiLV1 + JpsiLV2).Rapidity()));

            Jpsi_mass1.push_back(GEJpsi1_mass);
            Jpsi_pt1.push_back(GEJpsi1_pt);
            Jpsi_y1.push_back(GEJpsi1_y);

            Jpsi_mass2.push_back(GEJpsi2_mass);
            Jpsi_pt2.push_back(GEJpsi2_pt);
            Jpsi_y2.push_back(GEJpsi2_y);
            totEvent++;
        }
        cout<<"Finished.  >> "<<endl;
        return;
    }

    public:
    int totEvent = 0;
    vector<Double_t> Jpsi_mass1, Jpsi_pt1, Jpsi_y1;
    vector<Double_t> Jpsi_mass2, Jpsi_pt2, Jpsi_y2;
    vector<Double_t> evt_mass, evt_y, evt_pt, delta_y, delta_phi;
    
    void loopOn() {
        string prefix = "/eos/home-c/chensh/JPsiPsi2s/HELAC_Onia/condorIO/JJSPS/ntuple/HO2016_Ntuple_";
        for(int i = 1; i <= 10; i++) {
            string filename = prefix + to_string(i) + ".root";
            readTree(filename);
        }
    }
};

string outFile = "UnweightNLO.root";
void slim() {
    Process process;
    process.loopOn();
    TFile file(outFile.c_str(), "RECREATE");
    TTree *outTree = new TTree("data", "data");
    // Define output tree variables
    Double_t Jpsi_mass1, Jpsi_pt1, Jpsi_y1;
    Double_t Jpsi_mass2, Jpsi_pt2, Jpsi_y2;
    Double_t evt_mass, evt_y, evt_pt, delta_y, delta_phi;
    Double_t evt_weight;
    // Set output tree SetBranchAddress address
    outTree->Branch("Jpsi_mass1", &Jpsi_mass1, "Jpsi_mass1/D");
    outTree->Branch("Jpsi_pt1", &Jpsi_pt1, "Jpsi_pt1/D");
    outTree->Branch("Jpsi_y1", &Jpsi_y1, "Jpsi_y1/D");

    outTree->Branch("Jpsi_mass2", &Jpsi_mass2, "Jpsi_mass2/D");
    outTree->Branch("Jpsi_pt2", &Jpsi_pt2, "Jpsi_pt2/D");
    outTree->Branch("Jpsi_y2", &Jpsi_y2, "Jpsi_y2/D");

    outTree->Branch("evt_mass", &evt_mass, "evt_mass/D");
    outTree->Branch("evt_y", &evt_y, "evt_y/D");
    outTree->Branch("evt_pt", &evt_pt, "evt_pt/D");
    outTree->Branch("delta_y", &delta_y, "delta_y/D");
    outTree->Branch("delta_phi", &delta_phi, "delta_phi/D");

    outTree->Branch("evt_weight", &evt_weight, "evt_weight/D");
    for(int i = 0; i < process.totEvent; i++) {
        Jpsi_mass1 = process.Jpsi_mass1[i];
        Jpsi_pt1 = process.Jpsi_pt1[i];
        Jpsi_y1 = process.Jpsi_y1[i];

        Jpsi_mass2 = process.Jpsi_mass2[i];
        Jpsi_pt2 = process.Jpsi_pt2[i];
        Jpsi_y2 = process.Jpsi_y2[i];
        
        evt_mass = process.evt_mass[i];
        evt_y = process.evt_y[i];
        evt_pt = process.evt_pt[i];
        delta_y = process.delta_y[i];
        delta_phi = process.delta_phi[i];

        evt_weight = 1;
        outTree->Fill();
    }
    outTree->Write();
    cout<<"Done! Total event number: "<<process.totEvent<<endl;
    return;
}
