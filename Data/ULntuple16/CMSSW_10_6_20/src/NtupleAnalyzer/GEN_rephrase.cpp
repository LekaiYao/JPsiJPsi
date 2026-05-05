// Author: Shiyang CHEN
// Description: NtupleAnalyzer for 2018 MINIAOD/Monte-Carlo
// Implementation: Take Ntuple as input, apply physical seletion and acceptance&efficiency correction.
// Output root file contains following information: event(mass, weight), dimuon(mass, pt, eta, phi, ctau), muon(pt)
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
using namespace std;
#define PI 3.14159265359
// I/O settings area
// Data
// #define N_DIR 7
// string prefix[N_DIR] = {"B", "C", "D", "E", "F", "G", "H"};
// string infix = "/Ntuple_2016_";
// int suffix[N_DIR] = {20, 9, 14, 3, 8, 29, 36};
// string outFile = "WeightData.root";
// SPS LO
// #define N_DIR 1
// string prefix[N_DIR] = {""};
// string infix = "/eos/home-l/leyao/26JJ/MC_Maker/HelacOnia2016/CMSSW_10_6_20/src/NTUPLE/LO/Ntuple_2016_SPSLO";
// int suffix[N_DIR] = {1};
// string outFile = "WeightSPSLO_GEN.root";
// SPS NLO
// #define N_DIR 1
// string prefix[N_DIR] = {""};
// string infix = "/eos/home-l/leyao/26JJ/MC_Maker/HelacOnia2016/CMSSW_10_6_20/src/NTUPLE/NLO/Ntuple_2016_SPSNLO";
// int suffix[N_DIR] = {1};
// string outFile = "WeightSPSNLO_GEN.root";
// SPS
#define N_DIR 1
string prefix[N_DIR] = {""};
string infix = "/eos/home-c/chensh/JPsiJPsi/SKIM_tightfilter/SPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc_eff/Ntuple_2016_SPS";
int suffix[N_DIR] = {1};
string outFile = "WeightSPS_GEN.root";
// DPS
// #define N_DIR 1
// string prefix[N_DIR] = {""};
// string infix = "/eos/home-c/chensh/JPsiPsi2s/SKIM_tightfilter/DPS/ULPythia2018/CMSSW_10_2_5/src/4mu_acc_eff/pTHat4/Ntuple_2018_DPS";//Ntuple/Ntuple_2018_DPS";
// int suffix[N_DIR] = {40};//{25};
// string outFile = "WeightDPS.root";
// Bdecay
// #define N_DIR 2
// string prefix[N_DIR] = {"", ""};
// string infix = "/eos/home-c/chensh/JPsiPsi2s/SKIM_tightfilter/Bdecay/ULPythia2018/CMSSW_10_2_5/src/";
// string subinfix[N_DIR] = {"JpsiPsi2S/Ntuple_2018_Bdecay", "Psi2SJpsi/Ntuple_2018_Bdecay"};
// int suffix[N_DIR] = {9, 9};
// string outFile = "WeightBdecay.root";

class Process {
    private:
    vector<Double_t> acc_pt, acc_y, eff_pt, eff_y;
    Double_t **nGen_Jpsi, **nAcc_Jpsi;
    Double_t **nBin_Jpsi, **nVtx_Jpsi, **nVtx_evt, **nTrg_evt;
    Double_t calWeight(Double_t Jpsi_pt1, Double_t Jpsi_y1, Double_t Jpsi_pt2, Double_t Jpsi_y2) {
        int i = upper_bound(acc_y.begin(), acc_y.end(), Jpsi_y1) - acc_y.begin() - 1, j = upper_bound(acc_pt.begin(), acc_pt.end(), Jpsi_pt1) - acc_pt.begin() - 1;// i-J/psi y, j-J/psi pT
        int k = upper_bound(acc_y.begin(), acc_y.end(), Jpsi_y2) - acc_y.begin() - 1, l = upper_bound(acc_pt.begin(), acc_pt.end(), Jpsi_pt2) - acc_pt.begin() - 1;// k-psi(2S) y, l-psi(2S) pT
        Double_t w = nGen_Jpsi[i][j] / nAcc_Jpsi[i][j] * nGen_Jpsi[k][l] / nAcc_Jpsi[k][l];
        i = upper_bound(eff_y.begin(), eff_y.end(), Jpsi_y1) - eff_y.begin() - 1;
        j = upper_bound(eff_pt.begin(), eff_pt.end(), Jpsi_pt1) - eff_pt.begin() - 1;
        k = upper_bound(eff_y.begin(), eff_y.end(), Jpsi_y2) - eff_y.begin() - 1;
        l = upper_bound(eff_pt.begin(), eff_pt.end(), Jpsi_pt2) - eff_pt.begin() - 1;
        // w = nBin_Jpsi[i][j] / nVtx_Jpsi[i][j] * nBin_psi2S[k][l] / nVtx_psi2S[k][l] * nVtx_evt[l][j] / nTrg_evt[l][j];
        w *= nBin_Jpsi[i][j] / nVtx_Jpsi[i][j] * nBin_Jpsi[k][l] / nVtx_Jpsi[k][l] * nVtx_evt[l][j] / nTrg_evt[l][j];
        return w;
    }
    void readTree(string& fileName) {
        // Handle root file
        TFile file(fileName.c_str(), "READ");
        TTree *inTree = (TTree *)file.Get("rootuple/oniaTree");
        int nEvent = inTree->GetEntries();
        // Define input tree variables
        vector<Double_t> *REmu_pt = 0;
        vector<Double_t> *REJpsi_pt = 0, *REJpsi_eta = 0, *REJpsi_y = 0, *REJpsi_phi = 0, *REJpsi_mass = 0, *REJpsi_ctau = 0;//, *REJpsi_sigLxy = 0;
        Double_t GEJpsi1_phi = 0, GEJpsi2_phi = 0;
        vector<int> *REJpsi_muId1 = 0;//, *REJpsi_muId2 = 0;
        vector<Double_t> *REevt_fourMuMass = 0, *REevt_vtxProb = 0;//, *REevt_d = 0, *REevt_L1muPtMax = 0;
        vector<int> *REevt_JpsiId1 = 0, *REevt_JpsiId2 = 0;
        vector<bool> *REevt_matchTrg = 0, *REevt_passHLT = 0, *REevt_fourMuFit = 0, *REevt_samePV = 0;
        // Set input tree SetBranchAddress address
        inTree->SetBranchAddress("REmu_pt", &REmu_pt);
        inTree->SetBranchAddress("REJpsi_pt", &REJpsi_pt);
        inTree->SetBranchAddress("REJpsi_eta", &REJpsi_eta);
        inTree->SetBranchAddress("REJpsi_y", &REJpsi_y);
        inTree->SetBranchAddress("REJpsi_phi", &REJpsi_phi);
        inTree->SetBranchAddress("REJpsi_mass", &REJpsi_mass);
        inTree->SetBranchAddress("REJpsi_ctau", &REJpsi_ctau);
        inTree->SetBranchAddress("GEJpsi1_phi", &GEJpsi1_phi);
        inTree->SetBranchAddress("GEJpsi2_phi", &GEJpsi2_phi);
        // inTree->SetBranchAddress("REJpsi_sigLxy", &REJpsi_sigLxy);
        inTree->SetBranchAddress("REJpsi_muId1", &REJpsi_muId1);
        // inTree->SetBranchAddress("REJpsi_muId2", &REJpsi_muId2);
        inTree->SetBranchAddress("REevt_fourMuMass", &REevt_fourMuMass);
        inTree->SetBranchAddress("REevt_vtxProb", &REevt_vtxProb);
        // inTree->SetBranchAddress("REevt_d", &REevt_d);
        inTree->SetBranchAddress("REevt_JpsiId1", &REevt_JpsiId1);
        inTree->SetBranchAddress("REevt_JpsiId2", &REevt_JpsiId2);
        inTree->SetBranchAddress("REevt_matchTrg", &REevt_matchTrg);
        inTree->SetBranchAddress("REevt_passHLT", &REevt_passHLT);
        inTree->SetBranchAddress("REevt_fourMuFit", &REevt_fourMuFit);
        inTree->SetBranchAddress("REevt_samePV", &REevt_samePV);
        // inTree->SetBranchAddress("REevt_L1muPtMax", &REevt_L1muPtMax);
        // Loop on input tree entries
        cout<<fileName<<" has events: "<<nEvent<<endl;
        for(int i = 0; i < nEvent; i++) {// && i < 100
            cout<<"Processing No."<<i<<'\r';
            inTree->GetEntry(i);
            // Apply cut: pass trigger matching && dimuon pT within [10, 40]
            for(int j = (int)REevt_matchTrg->size() - 1; j >= 0; j--) {
                if(!REevt_passHLT->at(j)) continue;
                if(!REevt_matchTrg->at(j)) continue;
                // if(REevt_vtxProb->at(j) == 0) continue;
                int JpsiId1 = REevt_JpsiId1->at(j), JpsiId2 = REevt_JpsiId2->at(j);
                if(REJpsi_pt->at(JpsiId1) > 40 || REJpsi_pt->at(JpsiId1) < 10) continue;
                if(REJpsi_pt->at(JpsiId2) > 40 || REJpsi_pt->at(JpsiId2) < 10) continue;
                if(!REevt_samePV->at(j)) continue;
                //// Choose one from following cut:
                // Additional cut for No.1 L1 muon pT
                // if(REevt_L1muPtMax->at(j) < 5.5) continue;
                // Additional cut for No.1 muon pT
                // if(REmu_pt->at(REJpsi_muId1->at(JpsiId)) < 5.5 && REmu_pt->at(REpsi2S_muId1->at(psi2SId)) < 5.5) continue;
                ////
                // Differential cross-section bin selection
                TLorentzVector JpsiLV1, JpsiLV2;
                JpsiLV1.SetPtEtaPhiM(REJpsi_pt->at(JpsiId1), REJpsi_eta->at(JpsiId1), REJpsi_phi->at(JpsiId1), REJpsi_mass->at(JpsiId1));
                JpsiLV2.SetPtEtaPhiM(REJpsi_pt->at(JpsiId2), REJpsi_eta->at(JpsiId2), REJpsi_phi->at(JpsiId2), REJpsi_mass->at(JpsiId2));

                // if(REevt_fourMuMass->at(j) <= 57.5 || REevt_fourMuMass->at(j) > 67.5) continue;
                evt_mass.push_back((JpsiLV1 + JpsiLV2).M());
                if((JpsiLV1 + JpsiLV2).M() < 7.5) continue;
                evt_mass2.push_back(REevt_fourMuMass->at(j));

                // double deltaPhi = PI - fabs(fabs(REJpsi_phi->at(JpsiId) - REpsi2S_phi->at(psi2SId)) - PI);
                // if(deltaPhi < 2 || deltaPhi >= 3) continue;
                delta_phi.push_back(PI - fabs(fabs(REJpsi_phi->at(JpsiId1) - REJpsi_phi->at(JpsiId2)) - PI));
                delta_phi_GEN.push_back(PI - fabs(fabs(GEJpsi1_phi - GEJpsi2_phi) - PI));

                // double deltaY = fabs(REJpsi_y->at(JpsiId) - REpsi2S_y->at(psi2SId));
                // if(deltaY < 2.0 || deltaY >= 2.4) continue;
                delta_y.push_back(fabs(REJpsi_y->at(JpsiId1) - REJpsi_y->at(JpsiId2)));

                // if((JpsiLV + psi2SLV).Pt() < 24 || (JpsiLV + psi2SLV).Pt() > 32) continue;
                evt_pt.push_back((JpsiLV1 + JpsiLV2).Pt());

                // if(fabs((JpsiLV + psi2SLV).Rapidity()) < 2.4 || fabs((JpsiLV + psi2SLV).Rapidity()) > 3) continue;
                evt_y.push_back(fabs((JpsiLV1 + JpsiLV2).Rapidity()));

                ////
                // vector<Double_t> temp_mu_pt;
                // temp_mu_pt.push_back(REmu_pt->at(REJpsi_muId1->at(JpsiId)));
                // temp_mu_pt.push_back(REmu_pt->at(REJpsi_muId2->at(JpsiId)));
                // temp_mu_pt.push_back(REmu_pt->at(REpsi2S_muId1->at(psi2SId)));
                // temp_mu_pt.push_back(REmu_pt->at(REpsi2S_muId2->at(psi2SId)));
                // mu_pt.push_back(temp_mu_pt);
                Jpsi_mass1.push_back(REJpsi_mass->at(JpsiId1));
                Jpsi_ctau1.push_back(REJpsi_ctau->at(JpsiId1));
                // Jpsi_sigLxy.push_back(REJpsi_sigLxy->at(JpsiId));
                Jpsi_pt1.push_back(REJpsi_pt->at(JpsiId1));
                Jpsi_y1.push_back(REJpsi_y->at(JpsiId1));
                // Jpsi_eta.push_back(REJpsi_eta->at(JpsiId));
                // Jpsi_phi.push_back(REJpsi_phi->at(JpsiId));
                Jpsi_mass2.push_back(REJpsi_mass->at(JpsiId2));
                Jpsi_ctau2.push_back(REJpsi_ctau->at(JpsiId2));
                // psi2S_sigLxy.push_back(REpsi2S_sigLxy->at(psi2SId));
                Jpsi_pt2.push_back(REJpsi_pt->at(JpsiId2));
                Jpsi_y2.push_back(REJpsi_y->at(JpsiId2));
                // psi2S_eta.push_back(REpsi2S_eta->at(psi2SId));
                // psi2S_phi.push_back(REpsi2S_phi->at(psi2SId));
                // // evt_d.push_back(REevt_d->at(j));
                evt_vtxProb.push_back(REevt_vtxProb->at(j));
                evt_weight.push_back(calWeight(REJpsi_pt->at(JpsiId1), REJpsi_y->at(JpsiId1), REJpsi_pt->at(JpsiId2), REJpsi_y->at(JpsiId2)));
                // // evt_weight.push_back(calWeight(REJpsi_pt->at(JpsiId), REJpsi_eta->at(JpsiId), REpsi2S_pt->at(psi2SId), REpsi2S_eta->at(psi2SId)));
                totEvent++;
                if(REevt_passHLT->at(j)) hltEvent++;
                // if(REevt_fourMuFit->at(j)) vtxEvent++;
                if(REevt_samePV->at(j)) vtxEvent++;
		break;
            }
        }
        totEntry += nEvent;
        // file.Close();
        cout<<"Finished.  >> "<<endl;
        return;
    }

    public:
    int totEvent = 0, hltEvent = 0, vtxEvent = 0, totEntry = 0;
    // vector<vector<Double_t>> mu_pt;
    vector<Double_t> Jpsi_mass1, Jpsi_ctau1, Jpsi_pt1, Jpsi_y1;//, Jpsi_sigLxy, Jpsi_eta, Jpsi_phi;
    vector<Double_t> Jpsi_mass2, Jpsi_ctau2, Jpsi_pt2, Jpsi_y2;//, psi2S_sigLxy, psi2S_eta, psi2S_phi;
    vector<Double_t> evt_weight, evt_vtxProb;//, evt_d;
    vector<Double_t> evt_mass, evt_y, evt_pt, delta_y, delta_phi, delta_phi_GEN, evt_mass2;
    void readMatrix() {
        string line;
        // Save acc&eff in arrays
        ifstream accFile("acceptance.txt");
        if(!accFile.is_open()) return;
        int acc_ptBin = 0, acc_yBin = 0, lineCnt = 0;
        while(getline(accFile, line)) {
            istringstream iss(line);
            if(lineCnt == 0) {
                iss>>acc_ptBin>>acc_yBin;
                nGen_Jpsi = new Double_t*[acc_yBin];
                nAcc_Jpsi = new Double_t*[acc_yBin];
            }else if(lineCnt == 1) {
                Double_t pt;
                for(int i = 0; i <= acc_ptBin; i++) {
                    iss>>pt;
                    acc_pt.push_back(pt);
                }
            }else if(lineCnt == 2) {
                Double_t y;
                for(int i = 0; i <= acc_yBin; i++) {
                    iss>>y;
                    acc_y.push_back(y);
                }
            }else if(lineCnt <= 2 + acc_yBin) {
                int i = lineCnt - 3;
                nAcc_Jpsi[i] = new Double_t[acc_ptBin];
                nGen_Jpsi[i] = new Double_t[acc_ptBin];
                for(int j = 0; j < acc_ptBin; j++) iss>>nAcc_Jpsi[i][j]>>nGen_Jpsi[i][j];
            }else break;
            lineCnt++;
        }
        accFile.close();
        ifstream effFile("efficiency_2SPS+DPS.txt");
        if(!effFile.is_open()) return;
        int eff_ptBin = 0, eff_yBin = 0;
        lineCnt = 0;
        while(getline(effFile, line)) {
            istringstream iss(line);
            if(lineCnt == 0) {
                iss>>eff_ptBin>>eff_yBin;
                nBin_Jpsi = new Double_t*[eff_yBin];
                nVtx_Jpsi = new Double_t*[eff_yBin];
                nVtx_evt = new Double_t*[eff_ptBin];
                nTrg_evt = new Double_t*[eff_ptBin];
            }else if(lineCnt == 1) {
                Double_t pt;
                for(int i = 0; i <= eff_ptBin; i++) {
                    iss>>pt;
                    eff_pt.push_back(pt);
                }
            }else if(lineCnt == 2) {
                Double_t y;
                for(int i = 0; i <= eff_yBin; i++) {
                    iss>>y;
                    eff_y.push_back(y);
                }
            }else if(lineCnt <= 2 + eff_yBin) {
                int i = lineCnt - 3;
                nVtx_Jpsi[i] = new Double_t[eff_ptBin];
                nBin_Jpsi[i] = new Double_t[eff_ptBin];
                for(int j = 0; j < eff_ptBin; j++) iss>>nVtx_Jpsi[i][j]>>nBin_Jpsi[i][j];
            }else if(lineCnt <= 2 + eff_yBin + eff_ptBin) {
                int i = lineCnt - 3 - eff_yBin;
                nVtx_evt[i] = new Double_t[eff_ptBin];
                nTrg_evt[i] = new Double_t[eff_ptBin];
                for(int j = 0; j < eff_ptBin; j++) iss>>nTrg_evt[i][j]>>nVtx_evt[i][j];
            }else break;
            lineCnt++;
        }
        effFile.close();
        return;
    }
    void loopOn() {
        for(int i = 0; i < N_DIR; i++) {
            for(int j = 1; j <= suffix[i]; j++) {
                string fileName = prefix[i] + infix + prefix[i] + "_" + to_string(j) + ".root";
                // string dir = "/eos/home-c/chensh/Data2018/A/Charmonium/2018A_Ntuple_chensh_v1/251210_052421/0000";
                // string fileName = dir + infix + prefix[i] + "_" + to_string(j) + ".root";
                // string fileName = prefix[i] + infix + subinfix[i] + "_" + to_string(j) + ".root";// For Bdecay
                readTree(fileName);
            }
        }
    }
};

void rephrase() {
    Process process;
    process.readMatrix();
    process.loopOn();
    TFile file(outFile.c_str(), "RECREATE");
    TTree *outTree = new TTree("data", "data");
    // Define output tree variables
    // vector<Double_t> mu_pt;
    Double_t Jpsi_mass1, Jpsi_ctau1, Jpsi_pt1, Jpsi_y1;//, Jpsi_sigLxy, Jpsi_eta, Jpsi_phi;
    Double_t Jpsi_mass2, Jpsi_ctau2, Jpsi_pt2, Jpsi_y2;//, psi2S_sigLxy, psi2S_eta, psi2S_phi;
    Double_t evt_weight, evt_vtxProb;//, evt_d;
    Double_t evt_mass, evt_y, evt_pt, delta_y, delta_phi, delta_phi_GEN, evt_mass2;
    // Set output tree SetBranchAddress address
    // outTree->Branch("mu_pt", &mu_pt);
    outTree->Branch("Jpsi_mass1", &Jpsi_mass1, "Jpsi_mass1/D");
    outTree->Branch("Jpsi_ctau1", &Jpsi_ctau1, "Jpsi_ctau1/D");
    // outTree->Branch("Jpsi_sigLxy", &Jpsi_sigLxy, "Jpsi_sigLxy/D");
    outTree->Branch("Jpsi_pt1", &Jpsi_pt1, "Jpsi_pt1/D");
    outTree->Branch("Jpsi_y1", &Jpsi_y1, "Jpsi_y1/D");
    // outTree->Branch("Jpsi_eta", &Jpsi_eta, "Jpsi_eta/D");
    // outTree->Branch("Jpsi_phi", &Jpsi_phi, "Jpsi_phi/D");
    outTree->Branch("Jpsi_mass2", &Jpsi_mass2, "Jpsi_mass2/D");
    outTree->Branch("Jpsi_ctau2", &Jpsi_ctau2, "Jpsi_ctau2/D");
    // outTree->Branch("psi2S_sigLxy", &psi2S_sigLxy, "psi2S_sigLxy/D");
    outTree->Branch("Jpsi_pt2", &Jpsi_pt2, "Jpsi_pt2/D");
    outTree->Branch("Jpsi_y2", &Jpsi_y2, "Jpsi_y2/D");
    // outTree->Branch("psi2S_eta", &psi2S_eta, "psi2S_eta/D");
    // outTree->Branch("psi2S_phi", &psi2S_phi, "psi2S_phi/D");
    // outTree->Branch("evt_d", &evt_d, "evt_d/D");
    outTree->Branch("evt_vtxProb", &evt_vtxProb, "evt_vtxProb/D");
    outTree->Branch("evt_weight", &evt_weight, "evt_weight/D");
    outTree->Branch("evt_mass", &evt_mass, "evt_mass/D");
    outTree->Branch("evt_mass2", &evt_mass2, "evt_mass2/D");
    outTree->Branch("evt_y", &evt_y, "evt_y/D");
    outTree->Branch("evt_pt", &evt_pt, "evt_pt/D");
    outTree->Branch("delta_y", &delta_y, "delta_y/D");
    outTree->Branch("delta_phi", &delta_phi, "delta_phi/D");
    outTree->Branch("delta_phi_GEN", &delta_phi_GEN, "delta_phi_GEN/D");
    Double_t maxEvt_weight = 0;
    for(int i = 0; i < process.totEvent; i++) {
        // Change here to accommendate differential bins
        // mu_pt = process.mu_pt[i];
        Jpsi_mass1 = process.Jpsi_mass1[i];
        Jpsi_ctau1 = process.Jpsi_ctau1[i];
        // Jpsi_sigLxy = process.Jpsi_sigLxy[i];
        Jpsi_pt1 = process.Jpsi_pt1[i];
        Jpsi_y1 = process.Jpsi_y1[i];
        // Jpsi_eta = process.Jpsi_eta[i];
        // Jpsi_phi = process.Jpsi_phi[i];
        Jpsi_mass2 = process.Jpsi_mass2[i];
        Jpsi_ctau2 = process.Jpsi_ctau2[i];
        // psi2S_sigLxy = process.psi2S_sigLxy[i];
        Jpsi_pt2 = process.Jpsi_pt2[i];
        Jpsi_y2 = process.Jpsi_y2[i];
        // psi2S_eta = process.psi2S_eta[i];
        // psi2S_phi = process.psi2S_phi[i];
        // evt_d = process.evt_d[i];
        evt_vtxProb = process.evt_vtxProb[i];
        evt_weight = process.evt_weight[i];
        evt_mass = process.evt_mass[i];
        evt_mass2 = process.evt_mass2[i];
        evt_y = process.evt_y[i];
        evt_pt = process.evt_pt[i];
        delta_y = process.delta_y[i];
        delta_phi = process.delta_phi[i];
        delta_phi_GEN = process.delta_phi_GEN[i];
        if(process.evt_weight[i] > maxEvt_weight) maxEvt_weight = process.evt_weight[i];
        outTree->Fill();
    }
    outTree->Write();
    // file.Close();
    cout<<"Done! Total event number: "<<process.totEvent<<endl;
    cout<<"HLT event number: "<<process.hltEvent<<endl;
    cout<<"Valid fourmu vertex event number: "<<process.vtxEvent<<endl;
    cout<<"Total entry number: "<<process.totEntry<<endl;
    cout<<"maxEvt_weight = "<<maxEvt_weight<<endl;
    return;
}
