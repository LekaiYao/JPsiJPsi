#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
using namespace std;
#define PI 3.14159265359

void loadFile(vector<string>& filenames) {
    string prefix = "/eos/home-l/leyao/26JJ/MC_Maker/HelacOnia2016/CMSSW_10_6_20/src/NTUPLE/NLO_gpt0p8/Ntuple_2016_SPSstar_";
    for(int i = 1; i <= 96; i++) filenames.push_back(prefix + to_string(i) + ".root");
}

vector<Double_t> eff_pt, eff_y;
Double_t **nBin_Jpsi, **nVtx_Jpsi, **nVtx_evt, **nTrg_evt;
Double_t calWeight(Double_t Jpsi_pt1, Double_t Jpsi_y1, Double_t Jpsi_pt2, Double_t Jpsi_y2) {
    int i = upper_bound(eff_y.begin(), eff_y.end(), Jpsi_y1) - eff_y.begin() - 1;
    int j = upper_bound(eff_pt.begin(), eff_pt.end(), Jpsi_pt1) - eff_pt.begin() - 1;
    int k = upper_bound(eff_y.begin(), eff_y.end(), Jpsi_y2) - eff_y.begin() - 1;
    int l = upper_bound(eff_pt.begin(), eff_pt.end(), Jpsi_pt2) - eff_pt.begin() - 1;
    Double_t w = nBin_Jpsi[i][j] / nVtx_Jpsi[i][j] * nBin_Jpsi[k][l] / nVtx_Jpsi[k][l] * nVtx_evt[l][j] / nTrg_evt[l][j];
    return w;
}
void loadEff() {
    string line;
    ifstream effFile("txt/efficiency_0_0.6.txt");
    if(!effFile.is_open()) return;
    int eff_ptBin = 0, eff_yBin = 0, lineCnt = 0;
    while(getline(effFile, line)) {
        istringstream iss(line);
        if(lineCnt == 0) {
            int a, b, c;
            iss>>eff_ptBin>>eff_yBin>>a>>b>>c;
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

void count() {
    // Handle input files and construct TChain
    vector<string> filenames;
    loadFile(filenames);
    TChain *ch = new TChain("rootuple/oniaTree");
    for(int i = 0; i < (int)filenames.size(); i++) {
        cout<<"Processing "<<(i + 1)<<"th file: "<<filenames[i];
        ch->Add(filenames[i].c_str());
        cout<<". Done!"<<'\n';
    }
    loadEff();
    // Define sub-regions
    int nVars = 5, nBins[] = {6, 8, 5, 9, 7};
    double *vars = new double[nVars];
    string varNames[] = {"delta_y", "delta_phi", "evt_y", "evt_pt", "evt_mass"};
    vector<double> varBins[] = {
        {0, 0.5, 1, 1.5, 2, 2.5, 4},
        {0, 0.3927, 0.7854, 1.1781, 1.5708, 1.9635, 2.3562, 2.7489, 3.1416},
        {0, 0.4, 0.8, 1.2, 1.6, 2},
        {0, 5, 10, 15, 20, 25, 30, 35, 40, 80},
        {7.5, 17.5, 27.5, 37.5, 47.5, 57.5, 67.5, 107.5}
    };
    double **nCount = new double*[nVars], **nWeight = new double*[nVars];
    for(int i = 0; i < nVars; i++) {
        nCount[i] = new double[nBins[i]];
        nWeight[i] = new double[nBins[i]];
        for(int j = 0; j < nBins[i]; j++) {
            nCount[i][j] = 0;
            nWeight[i][j] = 0;
        }
    }
    // Loop on tree entries
    UChar_t GEevt_valid = 0, GEevt_passAcc = 0;
    Double_t GEJpsi1_pt = 0, GEJpsi1_eta = 0, GEJpsi1_phi = 0, GEJpsi1_mass = 0, GEJpsi1_y = 0;
    Double_t GEJpsi2_pt = 0, GEJpsi2_eta = 0, GEJpsi2_phi = 0, GEJpsi2_mass = 0, GEJpsi2_y = 0;
    Double_t GEevt_fourMuMass = 0;
    vector<Double_t> *REJpsi_pt = 0, *REJpsi_y = 0, *REmu_pt = 0;
    vector<int> *REevt_JpsiId1 = 0, *REevt_JpsiId2 = 0;//, *REJpsi_muId1 = 0, *REJpsi_muId2 = 0, *REmu_pvAsc = 0;
    vector<bool> *REevt_passHLT = 0, *REevt_matchTrg = 0, *REevt_samePV = 0;
    ch->SetBranchAddress("GEevt_valid", &GEevt_valid);
    ch->SetBranchAddress("GEevt_passAcc", &GEevt_passAcc);
    ch->SetBranchAddress("GEJpsi1_pt", &GEJpsi1_pt);
    ch->SetBranchAddress("GEJpsi1_eta", &GEJpsi1_eta);
    ch->SetBranchAddress("GEJpsi1_phi", &GEJpsi1_phi);
    ch->SetBranchAddress("GEJpsi1_mass", &GEJpsi1_mass);
    ch->SetBranchAddress("GEJpsi1_y", &GEJpsi1_y);
    ch->SetBranchAddress("GEJpsi2_pt", &GEJpsi2_pt);
    ch->SetBranchAddress("GEJpsi2_eta", &GEJpsi2_eta);
    ch->SetBranchAddress("GEJpsi2_phi", &GEJpsi2_phi);
    ch->SetBranchAddress("GEJpsi2_mass", &GEJpsi2_mass);
    ch->SetBranchAddress("GEJpsi2_y", &GEJpsi2_y);
    ch->SetBranchAddress("GEevt_fourMuMass", &GEevt_fourMuMass);
    ch->SetBranchAddress("REJpsi_pt", &REJpsi_pt);
    ch->SetBranchAddress("REJpsi_y", &REJpsi_y);
    ch->SetBranchAddress("REmu_pt", &REmu_pt);
    ch->SetBranchAddress("REevt_JpsiId1", &REevt_JpsiId1);
    ch->SetBranchAddress("REevt_JpsiId2", &REevt_JpsiId2);
    ch->SetBranchAddress("REevt_passHLT", &REevt_passHLT);
    ch->SetBranchAddress("REevt_matchTrg", &REevt_matchTrg);
    ch->SetBranchAddress("REevt_samePV", &REevt_samePV);
    Int_t nEntry = ch->GetEntries(), nPassAcc = 0, nMatchTrg = 0;
    Double_t totWeight = 0;
    for(int i = 0; i < nEntry; i++) {
        ch->GetEntry(i);
        if(!(i & 511)) cout<<"Processing: "<<(i+1)<<"th entry\r";
        if(!GEevt_valid || !GEevt_passAcc) continue;
        if(GEJpsi1_pt < 10 || GEJpsi1_pt > 40 || GEJpsi2_pt < 10 || GEJpsi2_pt > 40) continue;
        if(fabs(GEJpsi1_y) > 2 || fabs(GEJpsi2_y) > 2) continue;
        if(GEevt_fourMuMass < 7.5) continue;
        // if(REmu_pt.size() < 4) continue;
        nPassAcc++;
        double w = 0;
        for(int j = (int)REevt_matchTrg->size() - 1; j >= 0; j--) {
            if(!REevt_passHLT->at(j)) continue;
            if(!REevt_matchTrg->at(j)) continue;
            if(!REevt_samePV->at(j)) continue;
            int JpsiId1 = REevt_JpsiId1->at(j), JpsiId2 = REevt_JpsiId2->at(j);
            if(REJpsi_pt->at(JpsiId1) > 40 || REJpsi_pt->at(JpsiId1) < 10) continue;
            if(REJpsi_pt->at(JpsiId2) > 40 || REJpsi_pt->at(JpsiId2) < 10) continue;
            w = calWeight(REJpsi_pt->at(JpsiId1), REJpsi_y->at(JpsiId1), REJpsi_pt->at(JpsiId2), REJpsi_y->at(JpsiId2));
            nMatchTrg++;
            totWeight += w;
            break;
        }
        // delta_y
        vars[0] = fabs(GEJpsi1_y - GEJpsi2_y);
        // delta_phi
        vars[1] = PI - fabs(fabs(GEJpsi1_phi - GEJpsi2_phi) - PI);
        TLorentzVector JpsiLV1, JpsiLV2;
        JpsiLV1.SetPtEtaPhiM(GEJpsi1_pt, GEJpsi1_eta, GEJpsi1_phi, GEJpsi1_mass);
        JpsiLV2.SetPtEtaPhiM(GEJpsi2_pt, GEJpsi2_eta, GEJpsi2_phi, GEJpsi2_mass);
        // evt_y
        vars[2] = fabs((JpsiLV1 + JpsiLV2).Rapidity());
        // evt_pt
        vars[3] = (JpsiLV1 + JpsiLV2).Pt();
        // evt_mass
        vars[4] = GEevt_fourMuMass;
        for(int j = 0; j < nVars; j++) {
            if(vars[j] < varBins[j][0] || vars[j] >= varBins[j][nBins[j]]) continue;
            for(int k = 1; k <= nBins[j]; k++) {
                if(vars[j] >= varBins[j][k]) continue;
                nCount[j][k-1]++;
                nWeight[j][k-1] += w;
                break;
            }
        }
    }
    cout<<"\nnPassAcc="<<nPassAcc<<"\nnMatchTrg="<<nMatchTrg<<"\ntotWeight="<<totWeight<<endl;
    for(int i = 0; i < nVars; i++) {
        cout<<varNames[i]<<": {";
        for(int j = 0; j < nBins[i]; j++) cout<<(nCount[i][j] - nWeight[i][j]) / nWeight[i][j]<<", ";
        cout<<'}'<<endl;
    }
    return;
}
