#include <iostream>
#include <fstream>
#include <utility>
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TLorentzVector.h"
using namespace std;
#define PI 3.14159265359

struct MixedJpsi {
    Double_t pt, eta, phi, mass, y, weight;
    Long64_t sourceEvent;
    bool passAcc;
};

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
    vector<Double_t> *Jpsi_pt = 0, *Jpsi_eta = 0, *Jpsi_phi = 0, *Jpsi_mass = 0, *Jpsi_y = 0;
    vector<vector<Double_t>> *Jpsi_mu_pt = 0, *Jpsi_mu_eta = 0;
    vector<pair<int, int>> *pairId = 0;
    ch->SetBranchAddress("GENjpsi_pt", &Jpsi_pt);
    ch->SetBranchAddress("GENjpsi_eta", &Jpsi_eta);
    ch->SetBranchAddress("GENjpsi_phi", &Jpsi_phi);
    ch->SetBranchAddress("GENjpsi_mass", &Jpsi_mass);
    ch->SetBranchAddress("GENjpsi_y", &Jpsi_y);
    ch->SetBranchAddress("GENjpsi_mu_pt", &Jpsi_mu_pt);
    ch->SetBranchAddress("GENjpsi_mu_eta", &Jpsi_mu_eta);
    ch->SetBranchAddress("GEN_pair_id", &pairId);

    // Keep the two slots selected by GEN_pair_id separate.  Their empirical
    // spectra are not exchangeable, so DPS mixing must preserve one J/psi
    // from each slot rather than draw twice from a merged pool.
    Long64_t nEntry = ch->GetEntries();
    vector<MixedJpsi> pools[2];
    for(Long64_t i = 0; i < nEntry; i++) {
        ch->GetEntry(i);
        if(!pairId || pairId->empty()) continue;
        int selected[] = {pairId->at(0).first, pairId->at(0).second};
        for(int slot = 0; slot < 2; slot++) {
            int j = selected[slot];
            if(j < 0 || j >= (int)Jpsi_pt->size()) continue;
            if(Jpsi_pt->at(j) < 10 || Jpsi_pt->at(j) > 40 || fabs(Jpsi_y->at(j)) > 2) continue;
            MixedJpsi candidate;
            candidate.pt = Jpsi_pt->at(j);
            candidate.eta = Jpsi_eta->at(j);
            candidate.phi = Jpsi_phi->at(j);
            candidate.mass = Jpsi_mass->at(j);
            candidate.y = Jpsi_y->at(j);
            candidate.weight = calJpsiWeight(candidate.pt, candidate.y);
            candidate.sourceEvent = i;
            candidate.passAcc = Jpsi_mu_pt->at(j)[0] > 3.5 && Jpsi_mu_pt->at(j)[1] > 3.5
                && fabs(Jpsi_mu_eta->at(j)[0]) < 2.4 && fabs(Jpsi_mu_eta->at(j)[1]) < 2.4;
            pools[slot].push_back(candidate);
        }
    }

    Long64_t nEvent = 0, nPassAcc = 0, nTried = 0;
    Double_t totWeight = 0;
    vector<Double_t> removedCount(nEntry, 0), removedWeight(nEntry, 0);
    for(size_t i1 = 0; i1 < pools[0].size(); i1++) {
        const MixedJpsi &jpsi1 = pools[0][i1];
        for(size_t i2 = 0; i2 < pools[1].size(); i2++) {
            const MixedJpsi &jpsi2 = pools[1][i2];
            if(jpsi1.sourceEvent == jpsi2.sourceEvent) continue;
            nTried++;
            TLorentzVector JpsiLV1, JpsiLV2;
            JpsiLV1.SetPtEtaPhiM(jpsi1.pt, jpsi1.eta, jpsi1.phi, jpsi1.mass);
            JpsiLV2.SetPtEtaPhiM(jpsi2.pt, jpsi2.eta, jpsi2.phi, jpsi2.mass);
            if((JpsiLV1 + JpsiLV2).M() < 7.5) continue;
            nEvent++;
            double w = 0;
            if(jpsi1.passAcc && jpsi2.passAcc) w = jpsi1.weight * jpsi2.weight;
            nPassAcc += (int)(w != 0);
            totWeight += w;
            removedCount[jpsi1.sourceEvent] += 1;
            removedCount[jpsi2.sourceEvent] += 1;
            removedWeight[jpsi1.sourceEvent] += w;
            removedWeight[jpsi2.sourceEvent] += w;
            // delta_y
            vars[0] = fabs(jpsi1.y - jpsi2.y);
            // delta_phi
            vars[1] = PI - fabs(fabs(jpsi1.phi - jpsi2.phi) - PI);
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
    }
    Double_t closure = (nEvent - totWeight) / totWeight;

    // Delete-one-source-event jackknife.  Mixed pairs are correlated because
    // each input J/psi is reused, so pair-level Poisson errors are invalid.
    Double_t jackknifeMean = 0;
    for(Long64_t i = 0; i < nEntry; i++) {
        Double_t count_i = nEvent - removedCount[i];
        Double_t weight_i = totWeight - removedWeight[i];
        jackknifeMean += (count_i - weight_i) / weight_i;
    }
    jackknifeMean /= nEntry;
    Double_t jackknifeSum = 0;
    for(Long64_t i = 0; i < nEntry; i++) {
        Double_t count_i = nEvent - removedCount[i];
        Double_t weight_i = totWeight - removedWeight[i];
        Double_t closure_i = (count_i - weight_i) / weight_i;
        jackknifeSum += (closure_i - jackknifeMean) * (closure_i - jackknifeMean);
    }
    Double_t closureStat = sqrt((nEntry - 1.) / nEntry * jackknifeSum);

    cout<<"pool1="<<pools[0].size()<<"\npool2="<<pools[1].size()
        <<"\nnTried="<<nTried<<"\nnEvent="<<nEvent<<"\nnPassAcc="<<nPassAcc
        <<"\ntotWeight="<<totWeight<<"\nclosure="<<closure
        <<"\nclosure_jackknife_stat="<<closureStat<<endl;
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
