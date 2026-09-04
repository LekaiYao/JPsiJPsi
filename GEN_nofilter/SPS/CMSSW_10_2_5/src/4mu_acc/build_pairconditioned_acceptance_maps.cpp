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

using namespace std;

namespace {

const vector<double> kCoarsePtEdges = {
    10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
    20, 21, 22, 23, 24, 26, 28, 30, 35, 40};
const vector<double> kRefinedPtEdges = {
    10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14,
    15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 28, 30, 35, 40};
const vector<double> kQuarterPtEdges = {
    10, 10.25, 10.5, 10.75, 11, 11.25, 11.5, 11.75,
    12, 12.25, 12.5, 12.75, 13, 13.25, 13.5, 13.75, 14,
    15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 28, 30, 35, 40};
const vector<double> kPiecewisePtEdges = {
    10, 10.25, 10.5, 10.75, 11, 11.25, 11.5, 11.75,
    12, 12.25, 12.5, 12.75, 13, 13.25, 13.5, 13.75, 14,
    14.5, 15, 15.5, 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5,
    20, 20.5, 21, 21.5, 22, 22.5, 23, 23.5, 24, 24.5, 25,
    25.5, 26, 26.5, 27, 27.5, 28, 28.5, 29, 29.5, 30,
    31, 32, 33, 34, 35, 36, 37, 38, 39, 40};
const vector<double> kYEdges = {
    -2, -1.75, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 1.75, 2};

struct AcceptanceMap {
    vector<double> ptEdges;
    vector<Long64_t> fiducial;
    vector<Long64_t> accepted;

    explicit AcceptanceMap(const vector<double> &edges)
        : ptEdges(edges),
          fiducial((edges.size() - 1) * (kYEdges.size() - 1), 0),
          accepted((edges.size() - 1) * (kYEdges.size() - 1), 0) {}
};

vector<string> readInputList(const string &path) {
    ifstream input(path);
    if (!input.is_open()) throw runtime_error("Cannot open SPS input list: " + path);
    vector<string> files;
    string line;
    while (getline(input, line)) {
        if (!line.empty() && line[0] != '#') files.push_back(line);
    }
    if (files.empty()) throw runtime_error("Empty SPS input list: " + path);
    return files;
}

int findBin(double value, const vector<double> &edges) {
    if (value <= edges.front() || value >= edges.back()) return -1;
    return static_cast<int>(upper_bound(edges.begin(), edges.end(), value) - edges.begin()) - 1;
}

bool passFiducial(double pt, double y) {
    return pt > 10 && pt < 40 && y > -2 && y < 2;
}

bool passAcceptance(const vector<vector<Double_t>> &muPt,
                    const vector<vector<Double_t>> &muEta, int index) {
    if (index < 0 || index >= static_cast<int>(muPt.size())
        || index >= static_cast<int>(muEta.size()) || muPt.at(index).size() < 2
        || muEta.at(index).size() < 2) return false;
    return muPt.at(index)[0] > 3.5 && muPt.at(index)[1] > 3.5
        && fabs(muEta.at(index)[0]) < 2.4 && fabs(muEta.at(index)[1]) < 2.4;
}

void fillMap(AcceptanceMap &map, double pt, double y, bool accepted) {
    const int ptBin = findBin(pt, map.ptEdges);
    const int yBin = findBin(y, kYEdges);
    if (ptBin < 0 || yBin < 0) throw runtime_error("Fiducial J/psi is outside map binning");
    const size_t bin = static_cast<size_t>(yBin) * (map.ptEdges.size() - 1) + ptBin;
    ++map.fiducial[bin];
    if (accepted) ++map.accepted[bin];
}

void writeMap(const string &path, const AcceptanceMap &map,
              Long64_t nEvent, Long64_t nAcceptedEvent) {
    ofstream output(path);
    if (!output.is_open()) throw runtime_error("Cannot create acceptance table: " + path);
    output << (map.ptEdges.size() - 1) << ' ' << (kYEdges.size() - 1) << ' '
           << nEvent << ' ' << nAcceptedEvent << '\n' << setprecision(12);
    for (double edge : map.ptEdges) output << edge << ' ';
    output << '\n';
    for (double edge : kYEdges) output << edge << ' ';
    output << '\n';
    for (size_t iy = 0; iy + 1 < kYEdges.size(); ++iy) {
        for (size_t ipt = 0; ipt + 1 < map.ptEdges.size(); ++ipt) {
            const size_t bin = iy * (map.ptEdges.size() - 1) + ipt;
            output << map.accepted[bin] << ' ' << map.fiducial[bin] << ' ';
        }
        output << '\n';
    }
}

}  // namespace

void build_pairconditioned_acceptance_maps(
    const char *coarseOutput,
    const char *refinedOutput,
    const char *spsInputList,
    const char *quarterOutput = "",
    const char *piecewiseOutput = "") {
    TChain chain("GenAnalyzer/gen_tree");
    for (const string &path : readInputList(spsInputList)) {
        if (chain.Add(path.c_str()) == 0) throw runtime_error("Cannot add SPS input: " + path);
    }

    vector<Double_t> *jpsiPt = nullptr;
    vector<Double_t> *jpsiY = nullptr;
    vector<vector<Double_t>> *muPt = nullptr;
    vector<vector<Double_t>> *muEta = nullptr;
    vector<Double_t> *eventMass = nullptr;
    vector<pair<int, int>> *pairId = nullptr;
    chain.SetBranchAddress("GENjpsi_pt", &jpsiPt);
    chain.SetBranchAddress("GENjpsi_y", &jpsiY);
    chain.SetBranchAddress("GENjpsi_mu_pt", &muPt);
    chain.SetBranchAddress("GENjpsi_mu_eta", &muEta);
    chain.SetBranchAddress("GENevt_mass", &eventMass);
    chain.SetBranchAddress("GEN_pair_id", &pairId);

    AcceptanceMap coarse(kCoarsePtEdges);
    AcceptanceMap refined(kRefinedPtEdges);
    AcceptanceMap quarter(kQuarterPtEdges);
    AcceptanceMap piecewise(kPiecewisePtEdges);
    Long64_t nEvent = 0;
    Long64_t nAcceptedEvent = 0;
    Long64_t malformed = 0;
    const Long64_t entries = chain.GetEntries();
    for (Long64_t entry = 0; entry < entries; ++entry) {
        chain.GetEntry(entry);
        if (!pairId || pairId->empty() || !eventMass || eventMass->empty()
            || !jpsiPt || !jpsiY || !muPt || !muEta) {
            ++malformed;
            continue;
        }
        const int index[2] = {pairId->at(0).first, pairId->at(0).second};
        if (index[0] < 0 || index[1] < 0 || index[0] >= static_cast<int>(jpsiPt->size())
            || index[1] >= static_cast<int>(jpsiPt->size())
            || index[0] >= static_cast<int>(jpsiY->size())
            || index[1] >= static_cast<int>(jpsiY->size())) {
            ++malformed;
            continue;
        }
        if (eventMass->at(0) <= 7.5) continue;
        if (!passFiducial(jpsiPt->at(index[0]), jpsiY->at(index[0]))
            || !passFiducial(jpsiPt->at(index[1]), jpsiY->at(index[1]))) continue;
        ++nEvent;
        const bool accepted[2] = {
            passAcceptance(*muPt, *muEta, index[0]),
            passAcceptance(*muPt, *muEta, index[1])};
        if (accepted[0] && accepted[1]) ++nAcceptedEvent;
        for (int slot = 0; slot < 2; ++slot) {
            fillMap(coarse, jpsiPt->at(index[slot]), jpsiY->at(index[slot]), accepted[slot]);
            fillMap(refined, jpsiPt->at(index[slot]), jpsiY->at(index[slot]), accepted[slot]);
            fillMap(quarter, jpsiPt->at(index[slot]), jpsiY->at(index[slot]), accepted[slot]);
            fillMap(piecewise, jpsiPt->at(index[slot]), jpsiY->at(index[slot]), accepted[slot]);
        }
    }

    writeMap(coarseOutput, coarse, nEvent, nAcceptedEvent);
    writeMap(refinedOutput, refined, nEvent, nAcceptedEvent);
    if (quarterOutput && quarterOutput[0]) writeMap(quarterOutput, quarter, nEvent, nAcceptedEvent);
    if (piecewiseOutput && piecewiseOutput[0]) writeMap(piecewiseOutput, piecewise, nEvent, nAcceptedEvent);
    cout << setprecision(17)
         << "entries=" << entries << '\n'
         << "nEvent=" << nEvent << '\n'
         << "nAcceptedEvent=" << nAcceptedEvent << '\n'
         << "singleJpsiDenominator=" << 2 * nEvent << '\n'
         << "malformed=" << malformed << endl;
}
