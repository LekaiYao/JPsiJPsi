#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
using namespace std;

void loadFile(vector<string>& filenames) {
    string prefix = "/eos/home-c/chensh/JPsiJPsi/SKIM_tightfilter/DPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc_eff/";
    for(int i = 1; i <= 5; i++) {
        for(int j = 1; j <= 10; j++) filenames.push_back(prefix + to_string(i) + "/Ntuple_2016_DPS_" + to_string(j) + ".root");
    }
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
    ifstream effFile("/eos/home-c/chensh/JPsiJPsi/SKIM_tightfilter/SPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc_eff/txt/efficiency_0_0.6.txt");
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
    // loop on tree entries
    UChar_t GEevt_valid = 0, GEevt_passAcc = 0;
    Double_t GEJpsi1_pt = 0, GEJpsi1_y = 0, GEJpsi2_pt = 0, GEJpsi2_y = 0, GEevt_fourMuMass = 0;
    vector<Double_t> *REJpsi_pt = 0, *REJpsi_y = 0;
    vector<int> *REevt_JpsiId1 = 0, *REevt_JpsiId2 = 0;//, *REJpsi_muId1 = 0, *REJpsi_muId2 = 0, *REmu_pvAsc = 0;
    vector<bool> *REevt_passHLT = 0, *REevt_matchTrg = 0, *REevt_samePV = 0;
    ch->SetBranchAddress("GEevt_valid", &GEevt_valid);
    ch->SetBranchAddress("GEevt_passAcc", &GEevt_passAcc);
    ch->SetBranchAddress("GEJpsi1_pt", &GEJpsi1_pt);
    ch->SetBranchAddress("GEJpsi1_y", &GEJpsi1_y);
    ch->SetBranchAddress("GEJpsi2_pt", &GEJpsi2_pt);
    ch->SetBranchAddress("GEJpsi2_y", &GEJpsi2_y);
    ch->SetBranchAddress("GEevt_fourMuMass", &GEevt_fourMuMass);
    ch->SetBranchAddress("REJpsi_pt", &REJpsi_pt);
    ch->SetBranchAddress("REJpsi_y", &REJpsi_y);
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
        nPassAcc++;
        for(int j = (int)REevt_matchTrg->size() - 1; j >= 0; j--) {
            if(!REevt_passHLT->at(j)) continue;
            if(!REevt_matchTrg->at(j)) continue;
            if(!REevt_samePV->at(j)) continue;
            int JpsiId1 = REevt_JpsiId1->at(j), JpsiId2 = REevt_JpsiId2->at(j);
            if(REJpsi_pt->at(JpsiId1) > 40 || REJpsi_pt->at(JpsiId1) < 10) continue;
            if(REJpsi_pt->at(JpsiId2) > 40 || REJpsi_pt->at(JpsiId2) < 10) continue;
            nMatchTrg++;
            totWeight += calWeight(REJpsi_pt->at(JpsiId1), REJpsi_y->at(JpsiId1), REJpsi_pt->at(JpsiId2), REJpsi_y->at(JpsiId2));
            break;
        }
    }
    cout<<"\nnPassAcc="<<nPassAcc<<"\nnMatchTrg="<<nMatchTrg<<"\ntotWeight="<<totWeight<<endl;
    return;
}