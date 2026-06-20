// Author: Shiyang CHEN
// Description: NtupleAnalyzer for 2016 MINIAOD/Monte-Carlo
// Implementation: Take Ntuple as input, apply physical seletion and acceptance&efficiency correction.
// Output root file contains following information: event(mass, weight), dimuon(mass, pt, eta, phi, ctau), muon(pt)
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
        TTree *inTree = (TTree *)file.Get("GenAnalyzer/gen_tree");
        int nEvent = inTree->GetEntries();
        // Define input tree variables
        vector<Double_t> *GENjpsi_pt = 0, *GENjpsi_eta = 0, *GENjpsi_phi = 0, *GENjpsi_y = 0, *GENjpsi_mass = 0;
        // Set input tree SetBranchAddress address
        inTree->SetBranchAddress("GENjpsi_pt", &GENjpsi_pt);
        inTree->SetBranchAddress("GENjpsi_eta", &GENjpsi_eta);
        inTree->SetBranchAddress("GENjpsi_phi", &GENjpsi_phi);
        inTree->SetBranchAddress("GENjpsi_y", &GENjpsi_y);
        inTree->SetBranchAddress("GENjpsi_mass", &GENjpsi_mass);
        // Loop on input tree entries
        cout<<fileName<<" has events: "<<nEvent<<endl;
        for(int i = 0; i < nEvent; i++) {
            cout<<"Processing No."<<i<<'\r';
            inTree->GetEntry(i);
            // Apply fiducial cut
            for(int j = 0; j < (int)GENjpsi_pt->size(); j++) {
                if(GENjpsi_pt->at(j) > 40 || GENjpsi_pt->at(j) < 10 || fabs(GENjpsi_y->at(j)) > 2) continue;
                Jpsi_pt.push_back(GENjpsi_pt->at(j));
                Jpsi_eta.push_back(GENjpsi_eta->at(j));
                Jpsi_phi.push_back(GENjpsi_phi->at(j));
                Jpsi_y.push_back(GENjpsi_y->at(j));
                Jpsi_mass.push_back(GENjpsi_mass->at(j));
                totJspi++;
            }
        }
        cout<<"Finished.  >> "<<endl;
        return;
    }

    public:
    int totJspi = 0;
    vector<Double_t> Jpsi_pt, Jpsi_eta, Jpsi_phi, Jpsi_y, Jpsi_mass;
    
    void loopOn() {
        string prefix = "direct/DPS_2016_JJ_";
        for(int i = 1; i <= 10; i++) {
            string filename = prefix + to_string(i) + ".root";
            readTree(filename);
        }
    }
};

string outFile = "UnweightDPS.root";
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
    int nEvent = 10000;
    srand(time(0));
    delete gRandom;
    gRandom = new TRandom3(rand());
    for(int i = 0; i < nEvent; i++) {
        int i1 = (int)(gRandom->Rndm() * process.totJspi), i2 = (int)(gRandom->Rndm() * process.totJspi);
        if(i1 == i2) continue;

        Jpsi_mass1 = process.Jpsi_mass[i1];
        Jpsi_pt1 = process.Jpsi_pt[i1];
        Jpsi_y1 = process.Jpsi_y[i1];

        Jpsi_mass2 = process.Jpsi_mass[i2];
        Jpsi_pt2 = process.Jpsi_pt[i2];
        Jpsi_y2 = process.Jpsi_y[i2];

        TLorentzVector JpsiLV1, JpsiLV2;
        JpsiLV1.SetPtEtaPhiM(Jpsi_pt1, process.Jpsi_eta[i1], process.Jpsi_phi[i1], Jpsi_mass1);
        JpsiLV2.SetPtEtaPhiM(Jpsi_pt2, process.Jpsi_eta[i2], process.Jpsi_phi[i2], Jpsi_mass2);
        if((JpsiLV1 + JpsiLV2).M() < 7.5) continue;
        
        evt_mass = (JpsiLV1 + JpsiLV2).M();
        evt_y = fabs((JpsiLV1 + JpsiLV2).Rapidity());
        evt_pt = (JpsiLV1 + JpsiLV2).Pt();
        delta_y = fabs(Jpsi_y1 - Jpsi_y2);
        delta_phi = PI - fabs(fabs(process.Jpsi_phi[i1] - process.Jpsi_phi[i2]) - PI);

        evt_weight = 1;
        outTree->Fill();
    }
    outTree->Write();
    cout<<"Done! Total event number: "<<outTree->GetEntries()<<endl;
    return;
}
