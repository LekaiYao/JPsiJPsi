#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
using namespace std;
#define PI 3.14159265359

void loadFile(vector<string>& filenames) {
    // GEN-only Pythia8 DPS samples (chensh open path; local direct/ is empty, read upstream directly)
    string prefix = "/eos/user/c/chensh/JPsiJPsi/GEN_nofilter/DPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc/direct/DPS_2016_JJ_";
    // string prefix = "direct/DPS_2016_JJ_";
    for(int i = 1; i <= 10; i++) filenames.push_back(prefix + to_string(i) + ".root");
}

vector<Double_t> acc_pt, acc_y;
Double_t **nGen_Jpsi, **nAcc_Jpsi;
Double_t calWeight(Double_t Jpsi_pt1, Double_t Jpsi_y1, Double_t Jpsi_pt2, Double_t Jpsi_y2) {
    int i = upper_bound(acc_y.begin(), acc_y.end(), Jpsi_y1) - acc_y.begin() - 1, j = upper_bound(acc_pt.begin(), acc_pt.end(), Jpsi_pt1) - acc_pt.begin() - 1;// i-J/psi y, j-J/psi pT
    int k = upper_bound(acc_y.begin(), acc_y.end(), Jpsi_y2) - acc_y.begin() - 1, l = upper_bound(acc_pt.begin(), acc_pt.end(), Jpsi_pt2) - acc_pt.begin() - 1;// k-psi(2S) y, l-psi(2S) pT
    Double_t w = nGen_Jpsi[i][j] / nAcc_Jpsi[i][j] * nGen_Jpsi[k][l] / nAcc_Jpsi[k][l];
    return w;
}
Double_t calJpsiWeight(Double_t Jpsi_pt1, Double_t Jpsi_y1) {
    int i = upper_bound(acc_y.begin(), acc_y.end(), Jpsi_y1) - acc_y.begin() - 1, j = upper_bound(acc_pt.begin(), acc_pt.end(), Jpsi_pt1) - acc_pt.begin() - 1;// i-J/psi y, j-J/psi pT
    Double_t w = nGen_Jpsi[i][j] / nAcc_Jpsi[i][j];
    return w;
}
void loadAcc(bool statis) {
    string line;
    // Save acc in arrays
    // ifstream accFile("acceptance_SPS+2DPS.txt");
    // DPS closure uses the SPS acceptance map (same map as efficiency closure, tests model independence)
    ifstream accFile("/eos/home-l/leyao/26JJ/JPsiJPsi/GEN_nofilter/SPS/CMSSW_10_2_5/src/4mu_acc/plot/acceptance.txt");
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
    return;
}

void count() {
    // Handle input files and construct TChain
    vector<string> filenames;
    loadFile(filenames);
    TChain *ch = new TChain("GenAnalyzer/gen_tree");
    for(int i = 0; i < (int)filenames.size(); i++) {
        cout<<"Processing "<<(i + 1)<<"th file: "<<filenames[i];
        ch->Add(filenames[i].c_str());
        cout<<". Done!"<<'\n';
    }
    loadAcc(false);
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
    vector<Double_t> *Jpsi_pt = 0, *Jpsi_eta = 0, *Jpsi_phi = 0, *Jpsi_mass = 0, *Jpsi_y = 0, *evt_mass = 0;
    vector<vector<Double_t>> *Jpsi_mu_pt = 0, *Jpsi_mu_eta = 0;
    ch->SetBranchAddress("GENjpsi_pt", &Jpsi_pt);
    ch->SetBranchAddress("GENjpsi_eta", &Jpsi_eta);
    ch->SetBranchAddress("GENjpsi_phi", &Jpsi_phi);
    ch->SetBranchAddress("GENjpsi_mass", &Jpsi_mass);
    ch->SetBranchAddress("GENjpsi_y", &Jpsi_y);
    ch->SetBranchAddress("GENjpsi_mu_pt", &Jpsi_mu_pt);
    ch->SetBranchAddress("GENjpsi_mu_eta", &Jpsi_mu_eta);
    ch->SetBranchAddress("GENevt_mass", &evt_mass);
    // Unified version
    Int_t nEntry = ch->GetEntries();
    vector<Double_t> J_pt, J_eta, J_phi, J_mass, J_y, J_w;
    vector<bool> J_passAcc;
    for(int i = 0; i < nEntry; i++) {
        ch->GetEntry(i);
        cout<<"Processing: "<<(i+1)<<"th entry\r";
        // if(Jpsi_pt->size() < 2) continue;
        // if(Jpsi_pt->at(0) < 10 || Jpsi_pt->at(0) > 40 || Jpsi_pt->at(1) < 10 || Jpsi_pt->at(1) > 40) continue;
        // if(fabs(Jpsi_y->at(0)) > 2 || fabs(Jpsi_y->at(1)) > 2) continue;
        for(int j = 0; j < Jpsi_pt->size() && j < 2; j++) {
            if(Jpsi_pt->at(j) < 10 || Jpsi_pt->at(j) > 40 || fabs(Jpsi_y->at(j)) > 2) continue;
            J_pt.push_back(Jpsi_pt->at(j));
            J_eta.push_back(Jpsi_eta->at(j));
            J_phi.push_back(Jpsi_phi->at(j));
            J_mass.push_back(Jpsi_mass->at(j));
            J_y.push_back(Jpsi_y->at(j));
            J_w.push_back(calJpsiWeight(Jpsi_pt->at(j), Jpsi_y->at(j)));
            bool pAc = Jpsi_mu_pt->at(j)[0] > 3.5 && Jpsi_mu_pt->at(j)[1] > 3.5 && fabs(Jpsi_mu_eta->at(j)[0]) < 2.4 && fabs(Jpsi_mu_eta->at(j)[1]) < 2.4;
            J_passAcc.push_back(pAc);
        }
    }
    cout<<endl;
    Int_t nPool = J_pt.size(), nDraw = 10000, nEvent = 0, nPassAcc = 0;
    Double_t totWeight = 0;
    srand(time(0));
    delete gRandom;
    gRandom = new TRandom3(rand());
    for(int i = 0; i < nDraw; i++) {
        cout<<"Generating: "<<(i+1)<<"th mixing event\r";
        int i1 = (int)(gRandom->Rndm() * nPool), i2 = (int)(gRandom->Rndm() * nPool);
        if(i1 == i2) continue;
        TLorentzVector JpsiLV1, JpsiLV2;
        JpsiLV1.SetPtEtaPhiM(J_pt[i1], J_eta[i1], J_phi[i1], J_mass[i1]);
        JpsiLV2.SetPtEtaPhiM(J_pt[i2], J_eta[i2], J_phi[i2], J_mass[i2]);
        if((JpsiLV1 + JpsiLV2).M() < 7.5) continue;
        nEvent++;
        double w = 0;
        if(J_passAcc[i1] && J_passAcc[i2]) w = J_w[i1] * J_w[i2];
        nPassAcc += (int)(w != 0);
        totWeight += w;
        // delta_y
        vars[0] = fabs(J_y[i1] - J_y[i2]);
        // delta_phi
        vars[1] = PI - fabs(fabs(J_phi[i1] - J_phi[i2]) - PI);
        // evt_y
        vars[2] = fabs((JpsiLV1 + JpsiLV2).Rapidity());
        // evt_pt
        vars[3] = (JpsiLV1 + JpsiLV2).Pt();
        // evt_mass
        vars[4] = (JpsiLV1 + JpsiLV2).M();
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
    cout<<"\nnEvent="<<nEvent<<"\nnPassAcc="<<nPassAcc<<"\ntotWeight="<<totWeight<<endl;
    // Separated version
    // Int_t nEntry = ch->GetEntries(), nEvent = 0, nPassAcc = 0;
    // Double_t totWeight = 0;
    // for(int i = 0; i < nEntry; i++) {
    //     ch->GetEntry(i);
    //     cout<<"Processing: "<<(i+1)<<"th entry\r";
    //     if(Jpsi_pt->empty()) continue;
    //     if(Jpsi_pt->at(0) < 10 || Jpsi_pt->at(0) > 40 || Jpsi_pt->at(1) < 10 || Jpsi_pt->at(1) > 40) continue;
    //     if(fabs(Jpsi_y->at(0)) > 2 || fabs(Jpsi_y->at(1)) > 2) continue;
    //     if(evt_mass->at(0) < 7.5) continue;
    //     nEvent++;
    //     if(Jpsi_mu_pt->at(0)[0] > 3.5 && Jpsi_mu_pt->at(0)[1] > 3.5 && fabs(Jpsi_mu_eta->at(0)[0]) < 2.4 && fabs(Jpsi_mu_eta->at(0)[1]) < 2.4) {
    //         nPassAcc++;
    //         totWeight += calWeight_Jpsi(Jpsi_pt->at(0), Jpsi_y->at(0));
    //     }
    // }
    // cout<<"\n[J/psi]\nnEvent="<<nEvent[0]<<"\nnPassAcc="<<nPassAcc[0]<<"\ntotWeight="<<totWeight[0];
    for(int i = 0; i < nVars; i++) {
        cout<<varNames[i]<<": {";
        for(int j = 0; j < nBins[i]; j++) {
            if(nWeight[i][j] > 0) cout<<(nCount[i][j] - nWeight[i][j]) / nWeight[i][j]<<", ";
            else cout<<'['<<nCount[i][j]<<"], ";
        }
        cout<<'}'<<endl;
    }
    return;
}
