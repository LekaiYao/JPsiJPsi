#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "TChain.h"
#include "TLorentzVector.h"

using namespace std;

namespace {

const vector<double> kPtEdges = {
    10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14,
    15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 28, 30, 35, 40};
const vector<double> kYEdges = {
    -2, -1.75, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 1.75, 2};

struct Candidate {
    double pt;
    double eta;
    double phi;
    double mass;
    Long64_t sourceEvent;
    bool passAcceptance;
    int bin;
};

vector<string> readInputList(const string &path) {
    ifstream input(path);
    if (!input.is_open()) throw runtime_error("Cannot open DPS input list: " + path);
    vector<string> files;
    string line;
    while (getline(input, line)) {
        if (!line.empty() && line[0] != '#') files.push_back(line);
    }
    if (files.empty()) throw runtime_error("DPS input list is empty: " + path);
    return files;
}

int findBin(double pt, double y) {
    if (pt <= kPtEdges.front() || pt >= kPtEdges.back()
        || y <= kYEdges.front() || y >= kYEdges.back()) return -1;
    const int ptBin = static_cast<int>(upper_bound(kPtEdges.begin(), kPtEdges.end(), pt)
                                       - kPtEdges.begin()) - 1;
    const int yBin = static_cast<int>(upper_bound(kYEdges.begin(), kYEdges.end(), y)
                                      - kYEdges.begin()) - 1;
    return yBin * static_cast<int>(kPtEdges.size() - 1) + ptBin;
}

bool passAcceptance(const vector<vector<Double_t>> &muPt,
                    const vector<vector<Double_t>> &muEta, int index) {
    return index >= 0 && index < static_cast<int>(muPt.size())
        && index < static_cast<int>(muEta.size()) && muPt.at(index).size() >= 2
        && muEta.at(index).size() >= 2
        && muPt.at(index)[0] > 3.5 && muPt.at(index)[1] > 3.5
        && fabs(muEta.at(index)[0]) < 2.4 && fabs(muEta.at(index)[1]) < 2.4;
}

void writeMap(const string &path, const vector<unsigned long long> &accepted,
              const vector<unsigned long long> &fiducial,
              unsigned long long nEvent, unsigned long long nAcceptedEvent) {
    ofstream output(path);
    if (!output.is_open()) throw runtime_error("Cannot create DPS acceptance table: " + path);
    output << kPtEdges.size() - 1 << ' ' << kYEdges.size() - 1 << ' '
           << nEvent << ' ' << nAcceptedEvent << '\n' << setprecision(12);
    for (double edge : kPtEdges) output << edge << ' ';
    output << '\n';
    for (double edge : kYEdges) output << edge << ' ';
    output << '\n';
    const size_t nPt = kPtEdges.size() - 1;
    for (size_t iy = 0; iy + 1 < kYEdges.size(); ++iy) {
        for (size_t ipt = 0; ipt < nPt; ++ipt) {
            const size_t bin = iy * nPt + ipt;
            output << accepted[bin] << ' ' << fiducial[bin] << ' ';
        }
        output << '\n';
    }
}

}  // namespace

void build_dps_allcross_acceptance_map(const char *outputPath,
                                       const char *inputListPath) {
    TChain chain("GenAnalyzer/gen_tree");
    const vector<string> files = readInputList(inputListPath);
    for (const string &path : files) {
        if (chain.Add(path.c_str()) == 0) throw runtime_error("Cannot add DPS input: " + path);
    }

    vector<Double_t> *jpsiPt = nullptr;
    vector<Double_t> *jpsiEta = nullptr;
    vector<Double_t> *jpsiPhi = nullptr;
    vector<Double_t> *jpsiMass = nullptr;
    vector<Double_t> *jpsiY = nullptr;
    vector<vector<Double_t>> *muPt = nullptr;
    vector<vector<Double_t>> *muEta = nullptr;
    vector<pair<int, int>> *pairId = nullptr;
    chain.SetBranchAddress("GENjpsi_pt", &jpsiPt);
    chain.SetBranchAddress("GENjpsi_eta", &jpsiEta);
    chain.SetBranchAddress("GENjpsi_phi", &jpsiPhi);
    chain.SetBranchAddress("GENjpsi_mass", &jpsiMass);
    chain.SetBranchAddress("GENjpsi_y", &jpsiY);
    chain.SetBranchAddress("GENjpsi_mu_pt", &muPt);
    chain.SetBranchAddress("GENjpsi_mu_eta", &muEta);
    chain.SetBranchAddress("GEN_pair_id", &pairId);

    vector<Candidate> pools[2];
    Long64_t malformed = 0;
    const Long64_t entries = chain.GetEntries();
    for (Long64_t entry = 0; entry < entries; ++entry) {
        chain.GetEntry(entry);
        if (!pairId || pairId->empty() || !jpsiPt || !jpsiEta || !jpsiPhi
            || !jpsiMass || !jpsiY || !muPt || !muEta) {
            ++malformed;
            continue;
        }
        const int selected[2] = {pairId->at(0).first, pairId->at(0).second};
        for (int slot = 0; slot < 2; ++slot) {
            const int index = selected[slot];
            if (index < 0 || index >= static_cast<int>(jpsiPt->size())
                || index >= static_cast<int>(jpsiEta->size())
                || index >= static_cast<int>(jpsiPhi->size())
                || index >= static_cast<int>(jpsiMass->size())
                || index >= static_cast<int>(jpsiY->size())
                || index >= static_cast<int>(muPt->size())
                || index >= static_cast<int>(muEta->size())) {
                ++malformed;
                continue;
            }
            const int bin = findBin(jpsiPt->at(index), jpsiY->at(index));
            if (bin < 0 || muPt->at(index).size() < 2 || muEta->at(index).size() < 2) continue;
            pools[slot].push_back({
                jpsiPt->at(index), jpsiEta->at(index), jpsiPhi->at(index),
                jpsiMass->at(index), entry,
                passAcceptance(*muPt, *muEta, index), bin});
        }
    }

    const size_t nBins = (kPtEdges.size() - 1) * (kYEdges.size() - 1);
    vector<unsigned long long> fiducial(nBins, 0);
    vector<unsigned long long> accepted(nBins, 0);
    unsigned long long nTried = 0;
    unsigned long long nEvent = 0;
    unsigned long long nAcceptedEvent = 0;
    for (const Candidate &first : pools[0]) {
        TLorentzVector firstP4;
        firstP4.SetPtEtaPhiM(first.pt, first.eta, first.phi, first.mass);
        for (const Candidate &second : pools[1]) {
            if (first.sourceEvent == second.sourceEvent) continue;
            ++nTried;
            TLorentzVector secondP4;
            secondP4.SetPtEtaPhiM(second.pt, second.eta, second.phi, second.mass);
            if ((firstP4 + secondP4).M() < 7.5) continue;
            ++nEvent;
            ++fiducial[first.bin];
            ++fiducial[second.bin];
            if (first.passAcceptance) ++accepted[first.bin];
            if (second.passAcceptance) ++accepted[second.bin];
            if (first.passAcceptance && second.passAcceptance) ++nAcceptedEvent;
        }
    }

    writeMap(outputPath, accepted, fiducial, nEvent, nAcceptedEvent);
    unsigned int emptyDenominatorBins = 0;
    unsigned int zeroNumeratorBins = 0;
    for (size_t bin = 0; bin < nBins; ++bin) {
        if (fiducial[bin] == 0) ++emptyDenominatorBins;
        if (accepted[bin] == 0) ++zeroNumeratorBins;
    }
    cout << setprecision(17)
         << "entries=" << entries << '\n'
         << "pool1=" << pools[0].size() << '\n'
         << "pool2=" << pools[1].size() << '\n'
         << "nTried=" << nTried << '\n'
         << "nEvent=" << nEvent << '\n'
         << "nAcceptedEvent=" << nAcceptedEvent << '\n'
         << "singleJpsiDenominator=" << 2 * nEvent << '\n'
         << "pairAcceptance=" << static_cast<double>(nAcceptedEvent) / nEvent << '\n'
         << "emptyComponentDenominatorBins=" << emptyDenominatorBins << '\n'
         << "zeroComponentNumeratorBins=" << zeroNumeratorBins << '\n'
         << "malformed=" << malformed << endl;
}
