// Description: Randomly mix LO and NLO ntuple events while preserving the input tree format.
// Run with: root -l -b -q mix.cpp
#include <algorithm>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TSystem.h"
#include "TTree.h"

using namespace std;

// I/O settings area. The file-name rule follows rephrase.cpp:
// prefix[i] + infix + prefix[i] + "_" + fileIndex + ".root"
#define N_LO_DIR 1
string loPrefix[N_LO_DIR] = {""};
string loInfix = "/eos/home-l/leyao/26JJ/MC_Maker/HelacOnia2016/CMSSW_10_6_20/src/NTUPLE/LO/Ntuple_2016_SPSLO";
int loSuffix[N_LO_DIR] = {1};

#define N_NLO_DIR 1
string nloPrefix[N_NLO_DIR] = {""};
string nloInfix = "/eos/home-l/leyao/26JJ/MC_Maker/HelacOnia2016/CMSSW_10_6_20/src/NTUPLE/NLO/Ntuple_2016_SPSNLO";
int nloSuffix[N_NLO_DIR] = {1};

Long64_t n_NLO = 100;
Double_t r = 1.0;  // n_LO = n_NLO * r
ULong64_t randomSeed = 20260505;
string outFile = "Ntuple_2016_SPS_LO_NLO_mix.root";

const char *treeName = "rootuple/oniaTree";

void addInputFiles(TChain &chain, string prefix[], const string &infix, int suffix[], int nDir, bool verbose = true) {
    for(int i = 0; i < nDir; i++) {
        for(int j = 1; j <= suffix[i]; j++) {
            string fileName = prefix[i] + infix + prefix[i] + "_" + to_string(j) + ".root";
            if(gSystem->AccessPathName(fileName.c_str())) {
                cout<<"Skip missing file: "<<fileName<<endl;
                continue;
            }
            chain.Add(fileName.c_str());
            if(verbose) cout<<"Add file: "<<fileName<<endl;
        }
    }
}

bool sameBranchLayout(TChain &aChain, TChain &bChain) {
    aChain.LoadTree(0);
    bChain.LoadTree(0);
    TTree *a = aChain.GetTree();
    TTree *b = bChain.GetTree();
    if(!a || !b) return false;
    TObjArray *aBranches = a->GetListOfBranches();
    TObjArray *bBranches = b->GetListOfBranches();
    if(aBranches->GetEntries() != bBranches->GetEntries()) {
        cout<<"Different branch counts: "<<aBranches->GetEntries()<<" vs "<<bBranches->GetEntries()<<endl;
        return false;
    }
    for(int i = 0; i < aBranches->GetEntries(); i++) {
        string name = aBranches->At(i)->GetName();
        if(!bBranches->FindObject(name.c_str())) {
            cout<<"Branch only in one input: "<<name<<endl;
            return false;
        }
    }
    return true;
}

vector<Long64_t> sampleEntries(Long64_t nTotal, Long64_t nPick, mt19937_64 &rng) {
    vector<Long64_t> entries(nTotal);
    iota(entries.begin(), entries.end(), 0);
    shuffle(entries.begin(), entries.end(), rng);
    entries.resize(nPick);
    sort(entries.begin(), entries.end());
    return entries;
}

void mix() {
    Long64_t n_LO = (Long64_t)(n_NLO * r);
    if(n_NLO < 0 || n_LO < 0) {
        cout<<"n_NLO and n_LO must be non-negative."<<endl;
        return;
    }

    TChain loChain(treeName);
    TChain nloChain(treeName);
    TChain mixedChain(treeName);
    addInputFiles(loChain, loPrefix, loInfix, loSuffix, N_LO_DIR);
    addInputFiles(nloChain, nloPrefix, nloInfix, nloSuffix, N_NLO_DIR);
    addInputFiles(mixedChain, nloPrefix, nloInfix, nloSuffix, N_NLO_DIR, false);
    addInputFiles(mixedChain, loPrefix, loInfix, loSuffix, N_LO_DIR, false);

    Long64_t loEntries = loChain.GetEntries();
    Long64_t nloEntries = nloChain.GetEntries();
    cout<<"Available LO events: "<<loEntries<<endl;
    cout<<"Available NLO events: "<<nloEntries<<endl;

    if(loEntries == 0 || nloEntries == 0) {
        cout<<"Input LO or NLO chain is empty."<<endl;
        return;
    }
    if(n_NLO > nloEntries) {
        cout<<"Requested NLO events exceed available entries: "<<n_NLO<<" > "<<nloEntries<<endl;
        return;
    }
    if(n_LO > loEntries) {
        cout<<"Requested LO events exceed available entries: "<<n_LO<<" > "<<loEntries<<endl;
        return;
    }
    if(!sameBranchLayout(loChain, nloChain)) {
        cout<<"LO and NLO branch layouts are not identical. Stop before mixing."<<endl;
        return;
    }

    mt19937_64 rng(randomSeed);
    vector<Long64_t> nloSelected = sampleEntries(nloEntries, n_NLO, rng);
    vector<Long64_t> loSelected = sampleEntries(loEntries, n_LO, rng);
    vector<Long64_t> selected;
    selected.reserve(nloSelected.size() + loSelected.size());
    for(Long64_t entry : nloSelected) selected.push_back(entry);
    for(Long64_t entry : loSelected) selected.push_back(nloEntries + entry);
    shuffle(selected.begin(), selected.end(), rng);

    TFile out(outFile.c_str(), "RECREATE");
    TDirectory *rootuple = out.mkdir("rootuple");
    rootuple->cd();

    mixedChain.LoadTree(0);
    TTree *outTree = mixedChain.CloneTree(0);
    if(!outTree) {
        cout<<"Failed to clone output tree."<<endl;
        return;
    }
    outTree->SetName("oniaTree");
    outTree->SetTitle("Tree of FourMuon");
    for(Long64_t entry : selected) {
        mixedChain.GetEntry(entry);
        outTree->Fill();
    }
    Long64_t mixedEntries = outTree->GetEntries();
    outTree->Write();
    out.Close();

    cout<<"Output file: "<<outFile<<endl;
    cout<<"Selected NLO events: "<<n_NLO<<endl;
    cout<<"Selected LO events: "<<n_LO<<endl;
    cout<<"Mixed events: "<<mixedEntries<<endl;
    cout<<"Random seed: "<<randomSeed<<endl;
}

void mix(Long64_t inputNNLO, Double_t inputRatio, const char *inputOutFile, ULong64_t inputSeed = 20260505) {
    n_NLO = inputNNLO;
    r = inputRatio;
    outFile = inputOutFile;
    randomSeed = inputSeed;
    mix();
}
