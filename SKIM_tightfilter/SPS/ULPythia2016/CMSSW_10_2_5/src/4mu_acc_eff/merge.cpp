#include <iostream>
#include <fstream>
using namespace std;

void merge() {
    vector<Double_t> eff_pt, eff_eta;
    Double_t **nBin_Jpsi, **nVtx_Jpsi, **nVtx_evt, **nTrg_evt;
    // vector<vector<Double_t>> nBin_Jpsi, nVtx_Jpsi, nVtx_evt, nTrg_evt;
    string line;
    // Read SPS efficiency map
    ifstream effFile("plot/efficiency.txt");
    if(!effFile.is_open()) return;
    int eff_ptBin = 0, eff_etaBin = 0, lineCnt = 0, n_SPS_J, n_SPS;
    while(getline(effFile, line)) {
        istringstream iss(line);
        if(lineCnt == 0) {
            iss>>eff_ptBin>>eff_etaBin>>n_SPS_J>>n_SPS;
            nBin_Jpsi = new Double_t*[eff_etaBin];
            nVtx_Jpsi = new Double_t*[eff_etaBin];
            nVtx_evt = new Double_t*[eff_ptBin];
            nTrg_evt = new Double_t*[eff_ptBin];
            // for(int i = 0; i < eff_etaBin) {
            //     nBin_Jpsi.push_back(vector<Double_t>());
            //     nVtx_Jpsi.push_back(vector<Double_t>());
            // }
            // for(int i = 0; i < eff_ptBin) {
            //     nVtx_evt.push_back(vector<Double_t>());
            //     nTrg_evt.push_back(vector<Double_t>());
            // }
        }else if(lineCnt == 1) {
            Double_t pt;
            for(int i = 0; i <= eff_ptBin; i++) {
                iss>>pt;
                eff_pt.push_back(pt);
                // cout<<pt<<' ';
            }
        }else if(lineCnt == 2) {
            Double_t eta;
            for(int i = 0; i <= eff_etaBin; i++) {
                iss>>eta;
                eff_eta.push_back(eta);
                // cout<<eta<<' ';
            }
        }else if(lineCnt <= 2 + eff_etaBin) {
            int i = lineCnt - 3;
            // Double_t x, y;
            nVtx_Jpsi[i] = new Double_t[eff_ptBin];
            nBin_Jpsi[i] = new Double_t[eff_ptBin];
            for(int j = 0; j < eff_ptBin; j++) {
                iss>>nVtx_Jpsi[i][j]>>nBin_Jpsi[i][j];
                // iss>>x>>y;
                // nVtx_Jpsi[i].push_back(x);
                // nBin_Jpsi[i].push_back(y);
            }
        }else if(lineCnt <= 2 + eff_etaBin + eff_ptBin) {
            int i = lineCnt - 3 - eff_etaBin;
            // Double_t x, y;
            nVtx_evt[i] = new Double_t[eff_ptBin];
            nTrg_evt[i] = new Double_t[eff_ptBin];
            for(int j = 0; j < eff_ptBin; j++) {
                iss>>nTrg_evt[i][j]>>nVtx_evt[i][j];
                // iss>>x>>y;
                // nVtx_Jpsi[i].push_back(x);
                // nBin_Jpsi[i].push_back(y);
            }
        }else break;
        lineCnt++;
    }
    effFile.close();
    // Add DPS efficiency map with a factor
    effFile = ifstream("/eos/home-c/chensh/JPsiJPsi/SKIM_tightfilter/DPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc_eff/plot/efficiency.txt");
    if(!effFile.is_open()) return;
    int n_DPS_J, n_DPS;
    eff_ptBin = 0;
    eff_etaBin = 0;
    lineCnt = 0;
    Double_t fac_DPS = 0.5;
    while(getline(effFile, line)) {
        istringstream iss(line);
        if(lineCnt == 0) {
            iss>>eff_ptBin>>eff_etaBin>>n_DPS_J>>n_DPS;
            fac_DPS *= (double)n_SPS / n_DPS * 0.18933 / 0.33072;
        }else if(lineCnt == 1) {
            Double_t pt;
            for(int i = 0; i <= eff_ptBin; i++) iss>>pt;
        }else if(lineCnt == 2) {
            Double_t eta;
            for(int i = 0; i <= eff_etaBin; i++) iss>>eta;
        }else if(lineCnt <= 2 + eff_etaBin) {
            int i = lineCnt - 3;
            Double_t nVtx_Jpsi_ij, nBin_Jpsi_ij;
            for(int j = 0; j < eff_ptBin; j++) {
                iss>>nVtx_Jpsi_ij>>nBin_Jpsi_ij;
                nVtx_Jpsi[i][j] += fac_DPS * nVtx_Jpsi_ij;
                nBin_Jpsi[i][j] += fac_DPS * nBin_Jpsi_ij;
            }
        }else if(lineCnt <= 2 + eff_etaBin + eff_ptBin) {
            int i = lineCnt - 3 - eff_etaBin;
            Double_t nVtx_evt_ij, nTrg_evt_ij;
            for(int j = 0; j < eff_ptBin; j++) {
                iss>>nTrg_evt_ij>>nVtx_evt_ij;
                nVtx_evt[i][j] += fac_DPS * nVtx_evt_ij;
                nTrg_evt[i][j] += fac_DPS * nTrg_evt_ij;
            }
        }else break;
        lineCnt++;
    }
    effFile.close();
    // Output merged efficiency map
    ofstream outFile("efficiency_2SPS+DPS.txt");
    outFile<<eff_ptBin<<' '<<eff_etaBin<<endl;
    for(int i = 0; i <= eff_ptBin; i++) outFile<<eff_pt[i]<<' ';
    outFile<<'\n';
    for(int i = 0; i <= eff_etaBin; i++) outFile<<eff_eta[i]<<' ';
    outFile<<'\n';
    for(int i = 0; i < eff_etaBin; i++) {
        for(int j = 0; j < eff_ptBin; j++) outFile<<nVtx_Jpsi[i][j]<<' '<<nBin_Jpsi[i][j]<<' ';
        outFile<<'\n';
    }
    for(int i = 0; i < eff_ptBin; i++) {
        for(int j = 0; j < eff_ptBin; j++) outFile<<nTrg_evt[i][j]<<' '<<nVtx_evt[i][j]<<' ';
        outFile<<'\n';
    }
    outFile.close();
    return;
}
