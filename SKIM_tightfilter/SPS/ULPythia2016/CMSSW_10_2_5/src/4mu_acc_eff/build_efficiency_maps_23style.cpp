#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
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
const vector<double> kYEdges = {
    -2, -1.75, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 1.75, 2};

struct Counts {
    vector<double> ptEdges;
    Long64_t entries = 0;
    Long64_t nAccEvt = 0;
    Long64_t nFidEvt = 0;
    vector<Long64_t> jBin, jReco, jId, jVtx;
    vector<Long64_t> evtVtx, evtHlt, evtTrg;

    explicit Counts(const vector<double> &edges) : ptEdges(edges) {
        const size_t nPt = ptEdges.size() - 1;
        const size_t nY = kYEdges.size() - 1;
        jBin.assign(nPt * nY, 0);
        jReco.assign(nPt * nY, 0);
        jId.assign(nPt * nY, 0);
        jVtx.assign(nPt * nY, 0);
        evtVtx.assign(nPt * nPt, 0);
        evtHlt.assign(nPt * nPt, 0);
        evtTrg.assign(nPt * nPt, 0);
    }
};

vector<string> readInputs(const string &path) {
    ifstream input(path);
    if (!input.is_open()) throw runtime_error("Cannot open input list: " + path);
    vector<string> files;
    string line;
    while (getline(input, line)) if (!line.empty() && line[0] != '#') files.push_back(line);
    if (files.empty()) throw runtime_error("Empty input list: " + path);
    return files;
}

int strictBin(double value, const vector<double> &edges) {
    if (value <= edges.front() || value >= edges.back()) return -1;
    const vector<double>::const_iterator found = lower_bound(edges.begin(), edges.end(), value);
    if (found != edges.end() && *found == value) return -1;
    return static_cast<int>(upper_bound(edges.begin(), edges.end(), value) - edges.begin()) - 1;
}

void fillCounts(Counts &counts, double pt1Value, double y1Value,
                double pt2Value, double y2Value, bool match1, bool id1,
                bool vtx1, bool match2, bool id2, bool vtx2,
                bool hasFourthMuon, bool passHlt, bool matchTrigger) {
    const int pt1 = strictBin(pt1Value, counts.ptEdges);
    const int pt2 = strictBin(pt2Value, counts.ptEdges);
    const int y1 = strictBin(y1Value, kYEdges);
    const int y2 = strictBin(y2Value, kYEdges);
    if (pt1 < 0 || pt2 < 0 || y1 < 0 || y2 < 0)
        throw runtime_error("Selected event falls outside or exactly on an internal map edge");
    const size_t nPt = counts.ptEdges.size() - 1;
    const size_t bins[2] = {
        static_cast<size_t>(y1) * nPt + pt1,
        static_cast<size_t>(y2) * nPt + pt2};
    const bool matched[2] = {match1, match2};
    const bool identified[2] = {id1, id2};
    const bool vertexed[2] = {vtx1, vtx2};
    for (int slot = 0; slot < 2; ++slot) {
        ++counts.jBin[bins[slot]];
        if (matched[slot] && hasFourthMuon) ++counts.jReco[bins[slot]];
        if (identified[slot] && hasFourthMuon) ++counts.jId[bins[slot]];
        if (vertexed[slot]) ++counts.jVtx[bins[slot]];
    }
    const size_t eventBin = static_cast<size_t>(pt2) * nPt + pt1;
    if (vtx1 && vtx2) ++counts.evtVtx[eventBin];
    if (passHlt) ++counts.evtHlt[eventBin];
    if (matchTrigger) ++counts.evtTrg[eventBin];
}

Long64_t sumCounts(const vector<Long64_t> &values) {
    Long64_t total = 0;
    for (Long64_t value : values) total += value;
    return total;
}

void writeRaw(const string &path, const Counts &counts, bool includeFid) {
    ofstream output(path);
    if (!output.is_open()) throw runtime_error("Cannot create raw map: " + path);
    const size_t nPt = counts.ptEdges.size() - 1;
    const size_t nY = kYEdges.size() - 1;
    output << nPt << ' ' << nY << ' ' << counts.nAccEvt << ' '
           << sumCounts(counts.jBin);
    if (includeFid) output << ' ' << counts.nFidEvt;
    output << '\n';
    for (double edge : counts.ptEdges) output << edge << ' ';
    output << '\n';
    for (double edge : kYEdges) output << edge << ' ';
    output << '\n';
    for (size_t iy = 0; iy < nY; ++iy) {
        for (size_t ip = 0; ip < nPt; ++ip) {
            const size_t bin = iy * nPt + ip;
            output << counts.jVtx[bin] << ' ' << counts.jId[bin] << ' '
                   << counts.jReco[bin] << ' ' << counts.jBin[bin] << ' ';
        }
        output << '\n';
    }
    for (size_t ip2 = 0; ip2 < nPt; ++ip2) {
        for (size_t ip1 = 0; ip1 < nPt; ++ip1) {
            const size_t bin = ip2 * nPt + ip1;
            output << counts.evtTrg[bin] << ' ' << counts.evtHlt[bin] << ' '
                   << counts.evtVtx[bin] << ' ';
        }
        output << '\n';
    }
}

}  // namespace

void build_efficiency_maps_23style(const char *sampleArg, const char *inputListArg,
                                   const char *coarseOutputArg,
                                   const char *refinedOutputArg) {
    const string sample(sampleArg ? sampleArg : "");
    if (sample != "sps" && sample != "dps")
        throw runtime_error("sample must be sps or dps");
    TChain chain("rootuple/oniaTree");
    for (const string &path : readInputs(inputListArg)) {
        if (chain.Add(path.c_str()) != 1) throw runtime_error("Cannot add input: " + path);
    }

    UChar_t eventValid = 0, eventPassAcc = 0;
    UChar_t jpsi1Match = 0, jpsi1Id = 0, jpsi1Vtx = 0;
    UChar_t jpsi2Match = 0, jpsi2Id = 0, jpsi2Vtx = 0;
    UChar_t eventHlt = 0, eventTrigger = 0;
    Double_t eventMass = 0, pt1 = 0, y1 = 0, pt2 = 0, y2 = 0;
    vector<double> *recoMuonPt = nullptr;
    chain.SetBranchStatus("*", 0);
#define ENABLE(branch, target) chain.SetBranchStatus(branch, 1); chain.SetBranchAddress(branch, target)
    ENABLE("GEevt_valid", &eventValid);
    ENABLE("GEevt_passAcc", &eventPassAcc);
    ENABLE("GEevt_fourMuMass", &eventMass);
    ENABLE("GEJpsi1_pt", &pt1);
    ENABLE("GEJpsi1_y", &y1);
    ENABLE("GEJpsi2_pt", &pt2);
    ENABLE("GEJpsi2_y", &y2);
    ENABLE("GEJpsi1_matchGEN", &jpsi1Match);
    ENABLE("GEJpsi1_passID", &jpsi1Id);
    ENABLE("GEJpsi1_passVtx", &jpsi1Vtx);
    ENABLE("GEJpsi2_matchGEN", &jpsi2Match);
    ENABLE("GEJpsi2_passID", &jpsi2Id);
    ENABLE("GEJpsi2_passVtx", &jpsi2Vtx);
    ENABLE("GEevt_passHLT", &eventHlt);
    ENABLE("GEevt_matchTrg", &eventTrigger);
    ENABLE("REmu_pt", &recoMuonPt);
#undef ENABLE

    Counts coarse(kCoarsePtEdges);
    Counts refined(kRefinedPtEdges);
    coarse.entries = refined.entries = chain.GetEntries();
    for (Long64_t entry = 0; entry < coarse.entries; ++entry) {
        chain.GetEntry(entry);
        const bool fid1 = pt1 > 10 && pt1 < 40 && y1 > -2 && y1 < 2;
        const bool fid2 = pt2 > 10 && pt2 < 40 && y2 > -2 && y2 < 2;
        if (eventValid && eventMass > 7.5 && fid1 && fid2) {
            ++coarse.nFidEvt;
            ++refined.nFidEvt;
        }
        if (!(eventValid && eventPassAcc && eventMass > 7.5 && fid1 && fid2)) continue;
        ++coarse.nAccEvt;
        ++refined.nAccEvt;
        const bool hasFourthMuon = recoMuonPt && recoMuonPt->size() > 3
            && recoMuonPt->at(3) != 0;
        fillCounts(coarse, pt1, y1, pt2, y2, jpsi1Match, jpsi1Id, jpsi1Vtx,
                   jpsi2Match, jpsi2Id, jpsi2Vtx, hasFourthMuon,
                   eventHlt, eventTrigger);
        fillCounts(refined, pt1, y1, pt2, y2, jpsi1Match, jpsi1Id, jpsi1Vtx,
                   jpsi2Match, jpsi2Id, jpsi2Vtx, hasFourthMuon,
                   eventHlt, eventTrigger);
    }

    writeRaw(coarseOutputArg, coarse, sample == "sps");
    writeRaw(refinedOutputArg, refined, sample == "sps");
    cout << "sample=" << sample << '\n'
         << "entries=" << coarse.entries << '\n'
         << "nAccEvt=" << coarse.nAccEvt << '\n'
         << "nFidEvt=" << coarse.nFidEvt << '\n'
         << "nJpsi=" << sumCounts(coarse.jBin) << endl;
}
