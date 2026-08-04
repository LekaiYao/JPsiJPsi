#include <algorithm>
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
using namespace std;

void loadFile(vector<string>& filenames) {
    string prefix = "/eos/home-l/leyao/26JJ/MC_Maker/HelacOnia2016/CMSSW_10_6_20/src/NTUPLE/NLO_gpt0p8/Ntuple_2016_SPSstar_";
    for(int i = 1; i <= 96; i++) filenames.push_back(prefix + to_string(i) + ".root");
}

void efficiencyPlot() {
    // Handle input files and construct TChain
    vector<string> filenames;
    loadFile(filenames);
    TChain *ch = new TChain("rootuple/oniaTree");
    for(int i = 0; i < (int)filenames.size(); i++) {
        cout<<"Processing "<<(i + 1)<<"th file: "<<filenames[i];
        ch->Add(filenames[i].c_str());
        cout<<". Done!"<<'\n';
    }
    // Define efficiency maps
    vector<double> ptBin = {5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 28, 30, 35, 40};
    vector<double> yBin = {-2.25, -2, -1.75, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 1.75, 2, 2.25};
    int nPtBin = ptBin.size() - 1, nYBin = yBin.size() - 1;
    int **mBin_Jpsi = new int*[nYBin], **mRec_Jpsi = new int*[nYBin], **mIdt_Jpsi = new int*[nYBin], **mVtx_Jpsi = new int*[nYBin];
    int **mVtx_evt = new int*[nPtBin], **mHlt_evt = new int*[nPtBin], **mTrg_evt = new int*[nPtBin];
    for(int i = 0; i < nYBin; i++) {
        mBin_Jpsi[i] = new int[nPtBin];
        mRec_Jpsi[i] = new int[nPtBin];
        mIdt_Jpsi[i] = new int[nPtBin];
        mVtx_Jpsi[i] = new int[nPtBin];
        for(int j = 0; j < nPtBin; j++) {
            mBin_Jpsi[i][j] = 0;
            mRec_Jpsi[i][j] = 0;
            mIdt_Jpsi[i][j] = 0;
            mVtx_Jpsi[i][j] = 0;
        }
    }
    for(int i = 0; i < nPtBin; i++) {
        mVtx_evt[i] = new int[nPtBin];
        mHlt_evt[i] = new int[nPtBin];
        mTrg_evt[i] = new int[nPtBin];
        for(int j = 0; j < nPtBin; j++) {
            mVtx_evt[i][j] = 0;
            mHlt_evt[i][j] = 0;
            mTrg_evt[i][j] = 0;
        }
    }
    // Define branches
    UChar_t GEevt_valid = 0, GEevt_passAcc = 0, GEevt_passHLT = 0, GEevt_matchTrg = 0;
    UChar_t GEJpsi1_matchGEN = 0, GEJpsi1_passID = 0, GEJpsi1_passVtx = 0;
    UChar_t GEJpsi2_matchGEN = 0, GEJpsi2_passID = 0, GEJpsi2_passVtx = 0;
    Double_t GEJpsi1_pt = 0, GEJpsi1_y = 0;//, GEJpsi1_eta = 0, GEJpsi1_phi = 0, GEJpsi1_mass = 0;
    Double_t GEJpsi2_pt = 0, GEJpsi2_y = 0;//, GEJpsi2_eta = 0, GEJpsi2_phi = 0, GEJpsi2_mass = 0;
    Double_t GEevt_fourMuMass = 0;
    vector<Double_t> *GEmu_eta = 0, *GEmu_pt = 0;
    vector<Double_t> *REmu_pt = 0;// Temporary solution for flag initial value bug
    // vector<Double_t> *REJpsi_pt = 0, *REJpsi_y = 0;
    // vector<int> *REevt_JpsiId1 = 0, *REevt_JpsiId2 = 0;
    // vector<bool> *REevt_passHLT = 0, *REevt_matchTrg = 0, *REevt_samePV = 0;
    ch->SetBranchAddress("GEevt_valid", &GEevt_valid);
    ch->SetBranchAddress("GEevt_passAcc", &GEevt_passAcc);
    ch->SetBranchAddress("GEevt_passHLT", &GEevt_passHLT);
    ch->SetBranchAddress("GEevt_matchTrg", &GEevt_matchTrg);
    ch->SetBranchAddress("GEJpsi1_pt", &GEJpsi1_pt);
    // ch->SetBranchAddress("GEJpsi1_eta", &GEJpsi1_eta);
    // ch->SetBranchAddress("GEJpsi1_phi", &GEJpsi1_phi);
    // ch->SetBranchAddress("GEJpsi1_mass", &GEJpsi1_mass);
    ch->SetBranchAddress("GEJpsi1_y", &GEJpsi1_y);
    ch->SetBranchAddress("GEJpsi1_matchGEN", &GEJpsi1_matchGEN);
    ch->SetBranchAddress("GEJpsi1_passID", &GEJpsi1_passID);
    ch->SetBranchAddress("GEJpsi1_passVtx", &GEJpsi1_passVtx);
    ch->SetBranchAddress("GEJpsi2_pt", &GEJpsi2_pt);
    // ch->SetBranchAddress("GEJpsi2_eta", &GEJpsi2_eta);
    // ch->SetBranchAddress("GEJpsi2_phi", &GEJpsi2_phi);
    // ch->SetBranchAddress("GEJpsi2_mass", &GEJpsi2_mass);
    ch->SetBranchAddress("GEJpsi2_y", &GEJpsi2_y);
    ch->SetBranchAddress("GEJpsi2_matchGEN", &GEJpsi2_matchGEN);
    ch->SetBranchAddress("GEJpsi2_passID", &GEJpsi2_passID);
    ch->SetBranchAddress("GEJpsi2_passVtx", &GEJpsi2_passVtx);
    ch->SetBranchAddress("GEevt_fourMuMass", &GEevt_fourMuMass);
    ch->SetBranchAddress("GEmu_eta", &GEmu_eta);
    ch->SetBranchAddress("GEmu_pt", &GEmu_pt);
    ch->SetBranchAddress("REmu_pt", &REmu_pt);
    // ch->SetBranchAddress("REJpsi_pt", &REJpsi_pt);
    // ch->SetBranchAddress("REJpsi_y", &REJpsi_y);
    // ch->SetBranchAddress("REevt_JpsiId1", &REevt_JpsiId1);
    // ch->SetBranchAddress("REevt_JpsiId2", &REevt_JpsiId2);
    // ch->SetBranchAddress("REevt_passHLT", &REevt_passHLT);
    // ch->SetBranchAddress("REevt_matchTrg", &REevt_matchTrg);
    // ch->SetBranchAddress("REevt_samePV", &REevt_samePV);
    Int_t nEntry = ch->GetEntries(), nAcc = 0, nJpsi = 0;
    // Loop on tree entries
    for(int i = 0; i < nEntry; i++) {
        ch->GetEntry(i);
        if(!(i & 511)) cout<<"Processing: "<<(i+1)<<"th entry\r";
        if(!GEevt_valid || !GEevt_passAcc) continue;
        // if(!GEevt_valid) continue;
        // if(fabs(GEmu_eta->at(0)) > 2.4 || fabs(GEmu_eta->at(1)) > 2.4 || fabs(GEmu_eta->at(2)) > 2.4 || fabs(GEmu_eta->at(3)) > 2.4) continue;
        // if(GEmu_pt->at(0) < 2 || GEmu_pt->at(1) < 2 || GEmu_pt->at(2) < 2 || GEmu_pt->at(3) < 2) continue;
        if(GEJpsi1_pt < ptBin[0] || GEJpsi1_pt >= ptBin[nPtBin] || GEJpsi2_pt < ptBin[0] || GEJpsi2_pt >= ptBin[nPtBin]) continue;
        if(fabs(GEJpsi1_y) > yBin[nYBin] || fabs(GEJpsi2_y) > yBin[nYBin]) continue;
        if(GEevt_fourMuMass < 7.5) continue;
        // Locate the event on the map
        int iY1 = upper_bound(yBin.begin(), yBin.end(), GEJpsi1_y) - yBin.begin() - 1;
        int iPt1 = upper_bound(ptBin.begin(), ptBin.end(), GEJpsi1_pt) - ptBin.begin() - 1;
        int iY2 = upper_bound(yBin.begin(), yBin.end(), GEJpsi2_y) - yBin.begin() - 1;
        int iPt2 = upper_bound(ptBin.begin(), ptBin.end(), GEJpsi2_pt) - ptBin.begin() - 1;
        if(REmu_pt->size() < 4) continue;
        nAcc++;
        nJpsi += 2;
        mBin_Jpsi[iY1][iPt1]++;
        mBin_Jpsi[iY2][iPt2]++;
        mRec_Jpsi[iY1][iPt1] += GEJpsi1_matchGEN;
        mRec_Jpsi[iY2][iPt2] += GEJpsi2_matchGEN;
        mIdt_Jpsi[iY1][iPt1] += GEJpsi1_passID;
        mIdt_Jpsi[iY2][iPt2] += GEJpsi2_passID;
        mVtx_Jpsi[iY1][iPt1] += GEJpsi1_passVtx;
        mVtx_Jpsi[iY2][iPt2] += GEJpsi2_passVtx;
        if(!GEJpsi1_passVtx || !GEJpsi2_passVtx) continue;
        mVtx_evt[iPt2][iPt1]++;
        mHlt_evt[iPt2][iPt1] += GEevt_passHLT;
        mTrg_evt[iPt2][iPt1] += GEevt_matchTrg;
    }
    ofstream outFile("plot/raw_efficiency_NLO.txt");
    if(!outFile.is_open()) return;
    outFile<<nPtBin<<' '<<nYBin<<' '<<nAcc<<' '<<nJpsi<<' '<<1<<endl;
    for(int i = 0; i <= nPtBin; i++) outFile<<ptBin[i]<<' ';
    outFile<<endl;
    for(int i = 0; i <= nYBin; i++) outFile<<yBin[i]<<' ';
    outFile<<endl;
    for(int i = 0; i < nYBin; i++) {
        for(int j = 0; j < nPtBin; j++) {
            outFile<<mVtx_Jpsi[i][j]<<' '<<mIdt_Jpsi[i][j]<<' '<<mRec_Jpsi[i][j]<<' '<<mBin_Jpsi[i][j]<<' ';
        }
        outFile<<endl;
    }
    for(int i = 0; i < nPtBin; i++) {
        for(int j = 0; j < nPtBin; j++) {
            outFile<<mTrg_evt[i][j]<<' '<<mHlt_evt[i][j]<<' '<<mVtx_evt[i][j]<<' ';
        }
        outFile<<endl;
    }
    return;
}
