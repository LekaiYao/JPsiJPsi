// MC validation sample maker for Appendix C (kinematics MC-vs-data comparison).
// Reads a raw oniaTree ntuple, applies the SAME selection as rephrase.cpp, and
// writes a "WeightedTree" in the rich format expected by the Validation plotting
// scripts (J/psi-pair-level + per-muon RECO branches + acc&eff Weight_sum).
//
// Purpose: regenerate the SPS MC sample with the nominal NLO* (HELAC-Onia) ntuple
// to replace the legacy Pythia8(LO) MC_SPS.root used in the original Appendix C.
//
// Usage:  root -l -q -b make_validation_mc.cpp
// Output: MC_SPS_NLOstar.root  (tree "WeightedTree")
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
using namespace std;
#define PI 3.14159265359

// ---------------- I/O settings ----------------
// Nominal NLO* (HELAC-Onia, gluon pT>0.8 GeV IR cut) sample, same input as the
// active SPSstar block of rephrase.cpp.
#define N_DIR 1
string prefix[N_DIR] = {""};
string infix = "/eos/home-l/leyao/26JJ/MC_Maker/HelacOnia2016/CMSSW_10_6_20/src/NTUPLE/NLO_gpt0p8/Ntuple_2016_SPSstar";
int suffix[N_DIR] = {1};                       // process files _1.._suffix
string outFile = "MC_SPS_NLOstar.root";
string accPath = "../acceptance.txt";
string effPath = "../efficiency_0_0.6.txt";

class Process {
    private:
    vector<Double_t> acc_pt, acc_y, eff_pt, eff_y;
    Double_t **nGen_Jpsi, **nAcc_Jpsi;
    Double_t **nBin_Jpsi, **nVtx_Jpsi, **nVtx_evt, **nTrg_evt;
    // identical acceptance&efficiency weight as rephrase.cpp::calWeight
    Double_t calWeight(Double_t pt1, Double_t y1, Double_t pt2, Double_t y2) {
        int i = upper_bound(acc_y.begin(), acc_y.end(), y1) - acc_y.begin() - 1, j = upper_bound(acc_pt.begin(), acc_pt.end(), pt1) - acc_pt.begin() - 1;
        int k = upper_bound(acc_y.begin(), acc_y.end(), y2) - acc_y.begin() - 1, l = upper_bound(acc_pt.begin(), acc_pt.end(), pt2) - acc_pt.begin() - 1;
        Double_t w = nGen_Jpsi[i][j] / nAcc_Jpsi[i][j] * nGen_Jpsi[k][l] / nAcc_Jpsi[k][l];
        i = upper_bound(eff_y.begin(), eff_y.end(), y1) - eff_y.begin() - 1;
        j = upper_bound(eff_pt.begin(), eff_pt.end(), pt1) - eff_pt.begin() - 1;
        k = upper_bound(eff_y.begin(), eff_y.end(), y2) - eff_y.begin() - 1;
        l = upper_bound(eff_pt.begin(), eff_pt.end(), pt2) - eff_pt.begin() - 1;
        w *= nBin_Jpsi[i][j] / nVtx_Jpsi[i][j] * nBin_Jpsi[k][l] / nVtx_Jpsi[k][l] * nVtx_evt[l][j] / nTrg_evt[l][j];
        return w;
    }
    void readTree(string& fileName) {
        TFile file(fileName.c_str(), "READ");
        TTree *inTree = (TTree *)file.Get("rootuple/oniaTree");
        if(!inTree) { cout<<"  skip (no tree): "<<fileName<<endl; return; }
        int nEvent = inTree->GetEntries();
        vector<Double_t> *REmu_pt = 0, *REmu_eta = 0, *REmu_phi = 0;
        vector<Double_t> *REJpsi_pt = 0, *REJpsi_eta = 0, *REJpsi_y = 0, *REJpsi_phi = 0, *REJpsi_mass = 0, *REJpsi_ctau = 0;
        vector<int> *REJpsi_muId1 = 0, *REJpsi_muId2 = 0;
        vector<Double_t> *REevt_fourMuMass = 0;
        vector<int> *REevt_JpsiId1 = 0, *REevt_JpsiId2 = 0;
        vector<bool> *REevt_matchTrg = 0, *REevt_passHLT = 0, *REevt_samePV = 0;
        inTree->SetBranchAddress("REmu_pt", &REmu_pt);
        inTree->SetBranchAddress("REmu_eta", &REmu_eta);
        inTree->SetBranchAddress("REmu_phi", &REmu_phi);
        inTree->SetBranchAddress("REJpsi_pt", &REJpsi_pt);
        inTree->SetBranchAddress("REJpsi_eta", &REJpsi_eta);
        inTree->SetBranchAddress("REJpsi_y", &REJpsi_y);
        inTree->SetBranchAddress("REJpsi_phi", &REJpsi_phi);
        inTree->SetBranchAddress("REJpsi_mass", &REJpsi_mass);
        inTree->SetBranchAddress("REJpsi_ctau", &REJpsi_ctau);
        inTree->SetBranchAddress("REJpsi_muId1", &REJpsi_muId1);
        inTree->SetBranchAddress("REJpsi_muId2", &REJpsi_muId2);
        inTree->SetBranchAddress("REevt_fourMuMass", &REevt_fourMuMass);
        inTree->SetBranchAddress("REevt_JpsiId1", &REevt_JpsiId1);
        inTree->SetBranchAddress("REevt_JpsiId2", &REevt_JpsiId2);
        inTree->SetBranchAddress("REevt_matchTrg", &REevt_matchTrg);
        inTree->SetBranchAddress("REevt_passHLT", &REevt_passHLT);
        inTree->SetBranchAddress("REevt_samePV", &REevt_samePV);
        cout<<fileName<<" has events: "<<nEvent<<endl;
        for(int i = 0; i < nEvent; i++) {
            inTree->GetEntry(i);
            for(int j = (int)REevt_matchTrg->size() - 1; j >= 0; j--) {
                if(!REevt_passHLT->at(j)) continue;
                if(!REevt_matchTrg->at(j)) continue;
                int Id1 = REevt_JpsiId1->at(j), Id2 = REevt_JpsiId2->at(j);
                if(REJpsi_pt->at(Id1) > 40 || REJpsi_pt->at(Id1) < 10) continue;
                if(REJpsi_pt->at(Id2) > 40 || REJpsi_pt->at(Id2) < 10) continue;
                if(!REevt_samePV->at(j)) continue;
                TLorentzVector J1, J2;
                J1.SetPtEtaPhiM(REJpsi_pt->at(Id1), REJpsi_eta->at(Id1), REJpsi_phi->at(Id1), REJpsi_mass->at(Id1));
                J2.SetPtEtaPhiM(REJpsi_pt->at(Id2), REJpsi_eta->at(Id2), REJpsi_phi->at(Id2), REJpsi_mass->at(Id2));
                if((J1 + J2).M() < 7.5) continue;
                // muon indices: J1 -> mu1,mu2 ; J2 -> mu3,mu4
                int m1 = REJpsi_muId1->at(Id1), m2 = REJpsi_muId2->at(Id1);
                int m3 = REJpsi_muId1->at(Id2), m4 = REJpsi_muId2->at(Id2);
                // J/psi-pair-level (same naming as data branches)
                v_FourMuonMass.push_back((J1 + J2).M());
                v_FourMuonPt.push_back((J1 + J2).Pt());
                v_FourMuonY.push_back(fabs((J1 + J2).Rapidity()));
                v_DeltaY.push_back(fabs(REJpsi_y->at(Id1) - REJpsi_y->at(Id2)));
                v_DeltaPhi.push_back(PI - fabs(fabs(REJpsi_phi->at(Id1) - REJpsi_phi->at(Id2)) - PI));
                // per-J RECO
                v_J1PtReco.push_back(REJpsi_pt->at(Id1));   v_J2PtReco.push_back(REJpsi_pt->at(Id2));
                v_J1yReco.push_back(REJpsi_y->at(Id1));     v_J2yReco.push_back(REJpsi_y->at(Id2));
                v_J1PhiReco.push_back(REJpsi_phi->at(Id1)); v_J2PhiReco.push_back(REJpsi_phi->at(Id2));
                v_J1Mass.push_back(REJpsi_mass->at(Id1));   v_J2Mass.push_back(REJpsi_mass->at(Id2));
                v_x1.push_back(REJpsi_ctau->at(Id1));       v_x2.push_back(REJpsi_ctau->at(Id2));
                // per-muon RECO
                v_Mu1PtReco.push_back(REmu_pt->at(m1)); v_Mu2PtReco.push_back(REmu_pt->at(m2));
                v_Mu3PtReco.push_back(REmu_pt->at(m3)); v_Mu4PtReco.push_back(REmu_pt->at(m4));
                v_Mu1EtaReco.push_back(REmu_eta->at(m1)); v_Mu2EtaReco.push_back(REmu_eta->at(m2));
                v_Mu3EtaReco.push_back(REmu_eta->at(m3)); v_Mu4EtaReco.push_back(REmu_eta->at(m4));
                v_Weight.push_back(calWeight(REJpsi_pt->at(Id1), REJpsi_y->at(Id1), REJpsi_pt->at(Id2), REJpsi_y->at(Id2)));
                totEvent++;
                break;
            }
        }
        cout<<"Finished.  >> kept so far: "<<totEvent<<endl;
        return;
    }

    public:
    int totEvent = 0;
    vector<Double_t> v_FourMuonMass, v_FourMuonPt, v_FourMuonY, v_DeltaY, v_DeltaPhi;
    vector<Double_t> v_J1PtReco, v_J2PtReco, v_J1yReco, v_J2yReco, v_J1PhiReco, v_J2PhiReco;
    vector<Double_t> v_J1Mass, v_J2Mass, v_x1, v_x2;
    vector<Double_t> v_Mu1PtReco, v_Mu2PtReco, v_Mu3PtReco, v_Mu4PtReco;
    vector<Double_t> v_Mu1EtaReco, v_Mu2EtaReco, v_Mu3EtaReco, v_Mu4EtaReco;
    vector<Double_t> v_Weight;
    void readMatrix() {
        string line;
        ifstream accFile(accPath.c_str());
        if(!accFile.is_open()) { cout<<"ERROR: cannot open "<<accPath<<endl; return; }
        int acc_ptBin = 0, acc_yBin = 0, lineCnt = 0;
        while(getline(accFile, line)) {
            istringstream iss(line);
            if(lineCnt == 0) { iss>>acc_ptBin>>acc_yBin; nGen_Jpsi = new Double_t*[acc_yBin]; nAcc_Jpsi = new Double_t*[acc_yBin]; }
            else if(lineCnt == 1) { Double_t pt; for(int i = 0; i <= acc_ptBin; i++) { iss>>pt; acc_pt.push_back(pt); } }
            else if(lineCnt == 2) { Double_t y; for(int i = 0; i <= acc_yBin; i++) { iss>>y; acc_y.push_back(y); } }
            else if(lineCnt <= 2 + acc_yBin) { int i = lineCnt - 3; nAcc_Jpsi[i] = new Double_t[acc_ptBin]; nGen_Jpsi[i] = new Double_t[acc_ptBin]; for(int j = 0; j < acc_ptBin; j++) iss>>nAcc_Jpsi[i][j]>>nGen_Jpsi[i][j]; }
            else break;
            lineCnt++;
        }
        accFile.close();
        ifstream effFile(effPath.c_str());
        if(!effFile.is_open()) { cout<<"ERROR: cannot open "<<effPath<<endl; return; }
        int eff_ptBin = 0, eff_yBin = 0; lineCnt = 0;
        while(getline(effFile, line)) {
            istringstream iss(line);
            if(lineCnt == 0) { iss>>eff_ptBin>>eff_yBin; nBin_Jpsi = new Double_t*[eff_yBin]; nVtx_Jpsi = new Double_t*[eff_yBin]; nVtx_evt = new Double_t*[eff_ptBin]; nTrg_evt = new Double_t*[eff_ptBin]; }
            else if(lineCnt == 1) { Double_t pt; for(int i = 0; i <= eff_ptBin; i++) { iss>>pt; eff_pt.push_back(pt); } }
            else if(lineCnt == 2) { Double_t y; for(int i = 0; i <= eff_yBin; i++) { iss>>y; eff_y.push_back(y); } }
            else if(lineCnt <= 2 + eff_yBin) { int i = lineCnt - 3; nVtx_Jpsi[i] = new Double_t[eff_ptBin]; nBin_Jpsi[i] = new Double_t[eff_ptBin]; for(int j = 0; j < eff_ptBin; j++) iss>>nVtx_Jpsi[i][j]>>nBin_Jpsi[i][j]; }
            else if(lineCnt <= 2 + eff_yBin + eff_ptBin) { int i = lineCnt - 3 - eff_yBin; nVtx_evt[i] = new Double_t[eff_ptBin]; nTrg_evt[i] = new Double_t[eff_ptBin]; for(int j = 0; j < eff_ptBin; j++) iss>>nTrg_evt[i][j]>>nVtx_evt[i][j]; }
            else break;
            lineCnt++;
        }
        effFile.close();
        return;
    }
    void loopOn() {
        for(int i = 0; i < N_DIR; i++)
            for(int j = 1; j <= suffix[i]; j++) {
                string fileName = infix + prefix[i] + "_" + to_string(j) + ".root";
                readTree(fileName);
            }
    }
};

void make_validation_mc() {
    Process p;
    p.readMatrix();
    p.loopOn();
    TFile file(outFile.c_str(), "RECREATE");
    TTree *t = new TTree("WeightedTree", "WeightedTree");
    Float_t FourMuonMass_, FourMuonPt_, FourMuonY_, DeltaY_, DeltaPhi_;
    Float_t J1PtReco_, J2PtReco_, J1yReco_, J2yReco_, J1PhiReco_, J2PhiReco_;
    Float_t J1Mass_, J2Mass_, x1_, x2_;
    Float_t Mu1PtReco_, Mu2PtReco_, Mu3PtReco_, Mu4PtReco_;
    Float_t Mu1EtaReco_, Mu2EtaReco_, Mu3EtaReco_, Mu4EtaReco_;
    Float_t Weight_sum;
    t->Branch("FourMuonMass_", &FourMuonMass_, "FourMuonMass_/F");
    t->Branch("FourMuonPt_", &FourMuonPt_, "FourMuonPt_/F");
    t->Branch("FourMuonY_", &FourMuonY_, "FourMuonY_/F");
    t->Branch("DeltaY_", &DeltaY_, "DeltaY_/F");
    t->Branch("DeltaPhi_", &DeltaPhi_, "DeltaPhi_/F");
    t->Branch("J1PtReco_", &J1PtReco_, "J1PtReco_/F");
    t->Branch("J2PtReco_", &J2PtReco_, "J2PtReco_/F");
    t->Branch("J1yReco_", &J1yReco_, "J1yReco_/F");
    t->Branch("J2yReco_", &J2yReco_, "J2yReco_/F");
    t->Branch("J1PhiReco_", &J1PhiReco_, "J1PhiReco_/F");
    t->Branch("J2PhiReco_", &J2PhiReco_, "J2PhiReco_/F");
    t->Branch("J1Mass_", &J1Mass_, "J1Mass_/F");
    t->Branch("J2Mass_", &J2Mass_, "J2Mass_/F");
    t->Branch("x1_", &x1_, "x1_/F");
    t->Branch("x2_", &x2_, "x2_/F");
    t->Branch("Mu1PtReco_", &Mu1PtReco_, "Mu1PtReco_/F");
    t->Branch("Mu2PtReco_", &Mu2PtReco_, "Mu2PtReco_/F");
    t->Branch("Mu3PtReco_", &Mu3PtReco_, "Mu3PtReco_/F");
    t->Branch("Mu4PtReco_", &Mu4PtReco_, "Mu4PtReco_/F");
    t->Branch("Mu1EtaReco_", &Mu1EtaReco_, "Mu1EtaReco_/F");
    t->Branch("Mu2EtaReco_", &Mu2EtaReco_, "Mu2EtaReco_/F");
    t->Branch("Mu3EtaReco_", &Mu3EtaReco_, "Mu3EtaReco_/F");
    t->Branch("Mu4EtaReco_", &Mu4EtaReco_, "Mu4EtaReco_/F");
    t->Branch("Weight_sum", &Weight_sum, "Weight_sum/F");
    for(int i = 0; i < p.totEvent; i++) {
        FourMuonMass_ = p.v_FourMuonMass[i]; FourMuonPt_ = p.v_FourMuonPt[i]; FourMuonY_ = p.v_FourMuonY[i];
        DeltaY_ = p.v_DeltaY[i]; DeltaPhi_ = p.v_DeltaPhi[i];
        J1PtReco_ = p.v_J1PtReco[i]; J2PtReco_ = p.v_J2PtReco[i];
        J1yReco_ = p.v_J1yReco[i]; J2yReco_ = p.v_J2yReco[i];
        J1PhiReco_ = p.v_J1PhiReco[i]; J2PhiReco_ = p.v_J2PhiReco[i];
        J1Mass_ = p.v_J1Mass[i]; J2Mass_ = p.v_J2Mass[i]; x1_ = p.v_x1[i]; x2_ = p.v_x2[i];
        Mu1PtReco_ = p.v_Mu1PtReco[i]; Mu2PtReco_ = p.v_Mu2PtReco[i]; Mu3PtReco_ = p.v_Mu3PtReco[i]; Mu4PtReco_ = p.v_Mu4PtReco[i];
        Mu1EtaReco_ = p.v_Mu1EtaReco[i]; Mu2EtaReco_ = p.v_Mu2EtaReco[i]; Mu3EtaReco_ = p.v_Mu3EtaReco[i]; Mu4EtaReco_ = p.v_Mu4EtaReco[i];
        Weight_sum = p.v_Weight[i];
        t->Fill();
    }
    t->Write();
    cout<<"Done! Output "<<outFile<<" with "<<p.totEvent<<" selected events."<<endl;
    return;
}
