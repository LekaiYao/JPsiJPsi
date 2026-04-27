#include <iostream>
#include <fstream>
using namespace std;

vector<Double_t> acc_pt, acc_eta;
Double_t **nGen_Jpsi_1, **nAcc_Jpsi_1;
Double_t **nGen_Jpsi_2, **nAcc_Jpsi_2;

void merge() {
    string line;
    // Read SPS acceptance map
    ifstream accFile("plot/acceptance.txt");
    if(!accFile.is_open()) return;
    int acc_ptBin = 0, acc_etaBin = 0, lineCnt = 0;
    while(getline(accFile, line)) {
        istringstream iss(line);
        if(lineCnt == 0) {
            iss>>acc_ptBin>>acc_etaBin;
            nGen_Jpsi_1 = new Double_t*[acc_etaBin];
            nAcc_Jpsi_1 = new Double_t*[acc_etaBin];
        }else if(lineCnt == 1) {
            Double_t pt;
            for(int i = 0; i <= acc_ptBin; i++) {
                iss>>pt;
                acc_pt.push_back(pt);
            }
        }else if(lineCnt == 2) {
            Double_t eta;
            for(int i = 0; i <= acc_etaBin; i++) {
                iss>>eta;
                acc_eta.push_back(eta);
            }
        }else if(lineCnt <= 2 + acc_etaBin) {
            int i = lineCnt - 3;
            nAcc_Jpsi_1[i] = new Double_t[acc_ptBin];
            nGen_Jpsi_1[i] = new Double_t[acc_ptBin];
            for(int j = 0; j < acc_ptBin; j++) iss>>nAcc_Jpsi_1[i][j]>>nGen_Jpsi_1[i][j];
        }else break;
        lineCnt++;
    }
    accFile.close();
    // Add DPS acceptance map with a factor
    accFile = ifstream("/eos/home-c/chensh/JPsiPsi2s/GEN_nofilter/DPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc/plot/acceptance.txt");
    if(!accFile.is_open()) return;
    lineCnt = 0;
    while(getline(accFile, line)) {
        istringstream iss(line);
        if(lineCnt == 0) {
            iss>>acc_ptBin>>acc_etaBin;
            nGen_Jpsi_2 = new Double_t*[acc_etaBin];
            nAcc_Jpsi_2 = new Double_t*[acc_etaBin];
        }else if(lineCnt == 1) {
            Double_t pt;
            for(int i = 0; i <= acc_ptBin; i++) {
                iss>>pt;
                acc_pt.push_back(pt);
            }
        }else if(lineCnt == 2) {
            Double_t eta;
            for(int i = 0; i <= acc_etaBin; i++) {
                iss>>eta;
                acc_eta.push_back(eta);
            }
        }else if(lineCnt <= 2 + acc_etaBin) {
            int i = lineCnt - 3;
            nAcc_Jpsi_2[i] = new Double_t[acc_ptBin];
            nGen_Jpsi_2[i] = new Double_t[acc_ptBin];
            for(int j = 0; j < acc_ptBin; j++) iss>>nAcc_Jpsi_2[i][j]>>nGen_Jpsi_2[i][j];
        }else break;
        lineCnt++;
    }
    accFile.close();
    Double_t Jpsi_SPS = 0, Jpsi_DPS = 0;
    for(int i = 0; i < acc_etaBin; i++) {
        for(int j = 0; j < acc_ptBin; j++) {
            Jpsi_SPS += nGen_Jpsi_1[i][j];
            Jpsi_DPS += nGen_Jpsi_2[i][j];
        }
    }
    for(int i = 0; i < acc_etaBin; i++) {
        for(int j = 0; j < acc_ptBin; j++) {
            nGen_Jpsi_1[i][j] += 2 * Jpsi_SPS / Jpsi_DPS * nGen_Jpsi_2[i][j];
            nAcc_Jpsi_1[i][j] += 2 * Jpsi_SPS / Jpsi_DPS * nAcc_Jpsi_2[i][j];
            nGen_psi2S_1[i][j] += 2 * psi2S_SPS / psi2S_DPS * nGen_psi2S_2[i][j];
            nAcc_psi2S_1[i][j] += 2 * psi2S_SPS / psi2S_DPS * nAcc_psi2S_2[i][j];
        }
    }
    // Output merged efficiency map
    ofstream outFile("acceptance_2SPS+DPS.txt");
    outFile<<acc_ptBin<<' '<<acc_etaBin<<endl;
    for(int i = 0; i <= acc_ptBin; i++) outFile<<acc_pt[i]<<' ';
    outFile<<'\n';
    for(int i = 0; i <= acc_etaBin; i++) outFile<<acc_eta[i]<<' ';
    outFile<<'\n';
    for(int i = 0; i < acc_etaBin; i++) {
        for(int j = 0; j < acc_ptBin; j++) outFile<<nAcc_Jpsi_1[i][j]<<' '<<nGen_Jpsi_1[i][j]<<' ';
        outFile<<'\n';
    }
    for(int i = 0; i < acc_etaBin; i++) {
        for(int j = 0; j < acc_ptBin; j++) outFile<<nAcc_psi2S_1[i][j]<<' '<<nGen_psi2S_1[i][j]<<' ';
        outFile<<'\n';
    }
    outFile.close();
    return;
}
