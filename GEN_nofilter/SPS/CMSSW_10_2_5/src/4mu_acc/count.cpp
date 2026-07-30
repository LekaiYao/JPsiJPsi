#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
using namespace std;
#define PI 3.14159265359

void loadFile(vector<string>& filenames) {
    // GEN-only Pythia8 SPS samples (chensh open path; local direct/ is empty, read upstream directly)
    string prefix = "/eos/user/c/chensh/JPsiJPsi/GEN_nofilter/SPS/CMSSW_10_2_5/src/4mu_acc/direct/SPS_2016_JJ_";
    // string prefix = "direct/SPS_2016_JJ_";
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
void loadAcc(bool statis) {
    string line;
    // Save acc in arrays
    // ifstream accFile("acceptance_SPS+2DPS.txt");
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
    Int_t nEntry = ch->GetEntries(), nEvent = 0, nPassAcc = 0;
    Double_t totWeight = 0;
    for(int i = 0; i < nEntry; i++) {
        ch->GetEntry(i);
        cout<<"Processing: "<<(i+1)<<"th entry\r";
        if(Jpsi_pt->at(0) < 10 || Jpsi_pt->at(0) > 40 || Jpsi_pt->at(1) < 10 || Jpsi_pt->at(1) > 40) continue;
        if(fabs(Jpsi_y->at(0)) > 2 || fabs(Jpsi_y->at(1)) > 2) continue;
        if(evt_mass->at(0) < 7.5) continue;
        nEvent++;
        double w = 0;
        if(Jpsi_mu_pt->at(0)[0] < 3.5 || Jpsi_mu_pt->at(0)[1] < 3.5 || Jpsi_mu_pt->at(1)[0] < 3.5 || Jpsi_mu_pt->at(1)[1] < 3.5) w = 0;
        else if(fabs(Jpsi_mu_eta->at(0)[0]) > 2.4 || fabs(Jpsi_mu_eta->at(0)[1]) > 2.4 || fabs(Jpsi_mu_eta->at(1)[0]) > 2.4 || fabs(Jpsi_mu_eta->at(1)[1]) > 2.4) w = 0;
        else w = calWeight(Jpsi_pt->at(0), Jpsi_y->at(0), Jpsi_pt->at(1), Jpsi_y->at(1));
        nPassAcc += (int)(w != 0);
        totWeight += w;
        // delta_y
        vars[0] = fabs(Jpsi_y->at(0) - Jpsi_y->at(1));
        // delta_phi
        vars[1] = PI - fabs(fabs(Jpsi_phi->at(0) - Jpsi_phi->at(1)) - PI);
        TLorentzVector JpsiLV1, JpsiLV2;
        JpsiLV1.SetPtEtaPhiM(Jpsi_pt->at(0), Jpsi_eta->at(0), Jpsi_phi->at(0), Jpsi_mass->at(0));
        JpsiLV2.SetPtEtaPhiM(Jpsi_pt->at(1), Jpsi_eta->at(1), Jpsi_phi->at(1), Jpsi_mass->at(1));
        // evt_y
        vars[2] = fabs((JpsiLV1 + JpsiLV2).Rapidity());
        // evt_pt
        vars[3] = (JpsiLV1 + JpsiLV2).Pt();
        // evt_mass
        vars[4] = evt_mass->at(0);
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
    // Int_t nEntry = ch->GetEntries(), nEvent[] = {0, 0}, nPassAcc[] = {0, 0};
    // Double_t totWeight[] = {0, 0};
    // for(int i = 0; i < nEntry && i < 60000; i++) {
    //     ch->GetEntry(i);
    //     cout<<"Processing: "<<(i+1)<<"th entry\r";
    //     if(Jpsi_pt->empty() || psi2S_pt->empty()) continue;
    //     if(Jpsi_pt->at(0) < 10 || Jpsi_pt->at(0) > 40 || psi2S_pt->at(0) < 10 || psi2S_pt->at(0) > 40) continue;
    //     if(fabs(Jpsi_y->at(0)) > 2 || fabs(psi2S_y->at(0)) > 2) continue;
    //     nEvent[0]++;
    //     if(Jpsi_mu_pt->at(0)[0] > 3.5 && Jpsi_mu_pt->at(0)[1] > 3.5 && fabs(Jpsi_mu_eta->at(0)[0]) < 2.4 && fabs(Jpsi_mu_eta->at(0)[1]) < 2.4) {
    //         nPassAcc[0]++;
    //         totWeight[0] += calWeight_Jpsi(Jpsi_pt->at(0), Jpsi_y->at(0));
    //     }
    //     nEvent[1]++;
    //     if(psi2S_mu_pt->at(0)[0] > 3.5 && psi2S_mu_pt->at(0)[1] > 3.5 && fabs(psi2S_mu_eta->at(0)[0]) < 2.4 && fabs(psi2S_mu_eta->at(0)[1]) < 2.4) {
    //         nPassAcc[1]++;
    //         totWeight[1] += calWeight_psi2S(psi2S_pt->at(0), psi2S_y->at(0));
    //     }
    // }
    // cout<<"\n[J/psi]\nnEvent="<<nEvent[0]<<"\nnPassAcc="<<nPassAcc[0]<<"\ntotWeight="<<totWeight[0];
    // cout<<"\n[psi(2S)]\nnEvent="<<nEvent[1]<<"\nnPassAcc="<<nPassAcc[1]<<"\ntotWeight="<<totWeight[1]<<endl;
    for(int i = 0; i < nVars; i++) {
        cout<<varNames[i]<<": {";
        for(int j = 0; j < nBins[i]; j++) {
            if(nWeight[i][j]) cout<<(nCount[i][j] - nWeight[i][j]) / nWeight[i][j]<<", ";
            else cout<<'['<<nCount[i][j]<<"], ";
        }
        cout<<'}'<<endl;
    }
    return;
}
