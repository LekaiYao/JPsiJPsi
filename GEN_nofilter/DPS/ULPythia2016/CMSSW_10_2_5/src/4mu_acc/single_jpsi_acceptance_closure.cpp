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

const vector<double> kPtEdges = {
    10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
    20, 21, 22, 23, 24, 26, 28, 30, 35, 40};
const vector<double> kYEdges = {
    -2, -1.75, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 1.75, 2};

struct SampleCounts {
    string label;
    Long64_t entries = 0;
    Long64_t emptyJpsiEvent = 0;
    Long64_t inconsistentJpsiVectors = 0;
    Long64_t jpsiObjects = 0;
    Long64_t malformedMuon = 0;
    Long64_t fiducial = 0;
    Long64_t accepted = 0;
    vector<Long64_t> binFiducial;
    vector<Long64_t> binAccepted;

    explicit SampleCounts(const string &name)
        : label(name),
          binFiducial((kPtEdges.size() - 1) * (kYEdges.size() - 1), 0),
          binAccepted((kPtEdges.size() - 1) * (kYEdges.size() - 1), 0) {}

    Long64_t totalFiducial() const { return fiducial; }
    Long64_t totalAccepted() const { return accepted; }
};

vector<string> readInputList(const string &path) {
    ifstream input(path);
    if (!input.is_open()) throw runtime_error("Cannot open input list: " + path);
    vector<string> files;
    string line;
    while (getline(input, line)) {
        if (!line.empty() && line[0] != '#') files.push_back(line);
    }
    if (files.empty()) throw runtime_error("Empty input list: " + path);
    return files;
}

int findBin(double value, const vector<double> &edges) {
    if (value < edges.front() || value >= edges.back()) return -1;
    return static_cast<int>(upper_bound(edges.begin(), edges.end(), value) - edges.begin()) - 1;
}

size_t flatIndex(int ptBin, int yBin) {
    return static_cast<size_t>(yBin) * (kPtEdges.size() - 1) + ptBin;
}

SampleCounts countSample(const string &label, const string &inputListPath) {
    const vector<string> files = readInputList(inputListPath);
    TChain chain("GenAnalyzer/gen_tree");
    for (const string &path : files) {
        if (chain.Add(path.c_str()) == 0) throw runtime_error("Cannot add input ROOT file: " + path);
    }

    vector<Double_t> *jpsiPt = nullptr;
    vector<Double_t> *jpsiY = nullptr;
    vector<vector<Double_t>> *muPt = nullptr;
    vector<vector<Double_t>> *muEta = nullptr;
    chain.SetBranchAddress("GENjpsi_pt", &jpsiPt);
    chain.SetBranchAddress("GENjpsi_y", &jpsiY);
    chain.SetBranchAddress("GENjpsi_mu_pt", &muPt);
    chain.SetBranchAddress("GENjpsi_mu_eta", &muEta);

    SampleCounts counts(label);
    counts.entries = chain.GetEntries();
    for (Long64_t entry = 0; entry < counts.entries; ++entry) {
        chain.GetEntry(entry);
        if (!jpsiPt || jpsiPt->empty()) {
            ++counts.emptyJpsiEvent;
            continue;
        }
        if (!jpsiY || jpsiY->size() != jpsiPt->size()) {
            ++counts.inconsistentJpsiVectors;
            continue;
        }
        counts.jpsiObjects += jpsiPt->size();
        for (size_t index = 0; index < jpsiPt->size(); ++index) {

            // This is deliberately a single-J/psi selection.  There is no
            // GEN_pair_id, partner-J/psi, or pair-mass requirement.
            const int ptBin = findBin(jpsiPt->at(index), kPtEdges);
            const int yBin = findBin(jpsiY->at(index), kYEdges);
            if (ptBin < 0 || yBin < 0) continue;
            const size_t bin = flatIndex(ptBin, yBin);
            ++counts.fiducial;
            ++counts.binFiducial[bin];

            if (!muPt || !muEta || index >= muPt->size()
                || index >= muEta->size() || muPt->at(index).size() < 2
                || muEta->at(index).size() < 2) {
                ++counts.malformedMuon;
                continue;
            }
            const bool pass = muPt->at(index)[0] > 3.5 && muPt->at(index)[1] > 3.5
                && fabs(muEta->at(index)[0]) < 2.4 && fabs(muEta->at(index)[1]) < 2.4;
            if (!pass) continue;
            ++counts.accepted;
            ++counts.binAccepted[bin];
        }
    }
    return counts;
}

void writeAcceptanceTable(const string &path, const SampleCounts &counts) {
    ofstream output(path);
    if (!output.is_open()) throw runtime_error("Cannot create acceptance table: " + path);
    output << (kPtEdges.size() - 1) << ' ' << (kYEdges.size() - 1) << ' '
           << counts.totalFiducial() << ' ' << counts.totalAccepted() << '\n';
    output << setprecision(12);
    for (double edge : kPtEdges) output << edge << ' ';
    output << '\n';
    for (double edge : kYEdges) output << edge << ' ';
    output << '\n';
    for (size_t iy = 0; iy + 1 < kYEdges.size(); ++iy) {
        for (size_t ipt = 0; ipt + 1 < kPtEdges.size(); ++ipt) {
            const size_t bin = flatIndex(static_cast<int>(ipt), static_cast<int>(iy));
            output << counts.binAccepted[bin] << ' ' << counts.binFiducial[bin] << ' ';
        }
        output << '\n';
    }
}

struct ClosureResult {
    double corrected = 0;
    double residual = 0;
    double statDps = 0;
    double statSpsMap = 0;
    double statTotal = 0;
    int unsupportedBins = 0;
};

ClosureResult calculateSpsSelfClosure(const SampleCounts &sps) {
    ClosureResult result;
    for (size_t bin = 0; bin < sps.binFiducial.size(); ++bin) {
        if (sps.binFiducial[bin] == 0) continue;
        if (sps.binAccepted[bin] == 0) {
            ++result.unsupportedBins;
            continue;
        }
        const double aSps = static_cast<double>(sps.binAccepted[bin]) / sps.binFiducial[bin];
        result.corrected += sps.binAccepted[bin] / aSps;
    }
    if (result.unsupportedBins != 0) throw runtime_error("SPS population found in unsupported SPS map bins");
    if (result.corrected <= 0) throw runtime_error("Non-positive corrected SPS single-J/psi count");
    result.residual = (sps.totalFiducial() - result.corrected) / result.corrected;
    return result;
}

ClosureResult calculateDpsClosure(const SampleCounts &sps, const SampleCounts &dps) {
    ClosureResult result;
    double varianceDps = 0;
    double varianceSpsMap = 0;
    for (size_t bin = 0; bin < sps.binFiducial.size(); ++bin) {
        if (dps.binFiducial[bin] == 0) continue;
        if (sps.binFiducial[bin] == 0 || sps.binAccepted[bin] == 0) {
            ++result.unsupportedBins;
            continue;
        }
        const double aSps = static_cast<double>(sps.binAccepted[bin]) / sps.binFiducial[bin];
        const double aDps = static_cast<double>(dps.binAccepted[bin]) / dps.binFiducial[bin];
        result.corrected += dps.binAccepted[bin] / aSps;
        varianceDps += dps.binFiducial[bin] * aDps * (1.0 - aDps) / (aSps * aSps);
        const double varianceAcceptance = aSps * (1.0 - aSps) / sps.binFiducial[bin];
        const double derivative = dps.binAccepted[bin] / (aSps * aSps);
        varianceSpsMap += derivative * derivative * varianceAcceptance;
    }
    if (result.unsupportedBins != 0) throw runtime_error("DPS population found in unsupported SPS map bins");
    if (result.corrected <= 0) throw runtime_error("Non-positive corrected DPS single-J/psi count");
    const double nDps = dps.totalFiducial();
    result.residual = (nDps - result.corrected) / result.corrected;
    const double closureDerivative = nDps / (result.corrected * result.corrected);
    result.statDps = closureDerivative * sqrt(varianceDps);
    result.statSpsMap = closureDerivative * sqrt(varianceSpsMap);
    result.statTotal = hypot(result.statDps, result.statSpsMap);
    return result;
}

void writeBinComparison(const string &path, const SampleCounts &sps, const SampleCounts &dps) {
    ofstream output(path);
    if (!output.is_open()) throw runtime_error("Cannot create bin comparison: " + path);
    output << "pt_low,pt_high,y_low,y_high,sps_fiducial,sps_accepted,dps_fiducial,dps_accepted,"
              "acceptance_sps,acceptance_dps,dps_over_sps_minus_one,ratio_stat\n";
    output << setprecision(17);
    for (size_t iy = 0; iy + 1 < kYEdges.size(); ++iy) {
        for (size_t ipt = 0; ipt + 1 < kPtEdges.size(); ++ipt) {
            const size_t bin = flatIndex(static_cast<int>(ipt), static_cast<int>(iy));
            const double aSps = sps.binFiducial[bin]
                ? static_cast<double>(sps.binAccepted[bin]) / sps.binFiducial[bin] : 0;
            const double aDps = dps.binFiducial[bin]
                ? static_cast<double>(dps.binAccepted[bin]) / dps.binFiducial[bin] : 0;
            output << kPtEdges[ipt] << ',' << kPtEdges[ipt + 1] << ','
                   << kYEdges[iy] << ',' << kYEdges[iy + 1] << ','
                   << sps.binFiducial[bin] << ',' << sps.binAccepted[bin] << ','
                   << dps.binFiducial[bin] << ',' << dps.binAccepted[bin] << ',';
            if (aSps > 0) output << aSps;
            output << ',';
            if (aDps > 0) output << aDps;
            output << ',';
            if (aSps > 0 && aDps > 0) {
                const double ratio = aDps / aSps;
                const double relativeVariance = (1.0 - aSps) / sps.binAccepted[bin]
                    + (1.0 - aDps) / dps.binAccepted[bin];
                output << (ratio - 1.0) << ',' << ratio * sqrt(relativeVariance);
            } else {
                output << ',';
            }
            output << '\n';
        }
    }
}

void writeSummary(const string &path, const SampleCounts &sps, const SampleCounts &dps,
                  const ClosureResult &spsClosure, const ClosureResult &dpsClosure) {
    ofstream output(path);
    if (!output.is_open()) throw runtime_error("Cannot create summary: " + path);
    output << "sps_entries,sps_jpsi_objects,sps_single_fiducial,sps_single_accepted,"
              "dps_entries,dps_jpsi_objects,dps_single_fiducial,dps_single_accepted,"
              "sps_raw_acceptance,dps_raw_acceptance,"
              "sps_corrected_single,sps_closure,dps_corrected_single,dps_closure,dps_closure_stat_dps,"
              "dps_closure_stat_sps_map,dps_closure_stat_total,unsupported_bins,"
              "sps_empty_jpsi_event,dps_empty_jpsi_event,sps_inconsistent_jpsi_vectors,"
              "dps_inconsistent_jpsi_vectors,sps_malformed_muon,dps_malformed_muon\n";
    output << setprecision(17)
           << sps.entries << ',' << sps.jpsiObjects << ',' << sps.totalFiducial() << ','
           << sps.totalAccepted() << ',' << dps.entries << ',' << dps.jpsiObjects << ','
           << dps.totalFiducial() << ',' << dps.totalAccepted() << ','
           << static_cast<double>(sps.totalAccepted()) / sps.totalFiducial() << ','
           << static_cast<double>(dps.totalAccepted()) / dps.totalFiducial() << ','
           << spsClosure.corrected << ',' << spsClosure.residual << ','
           << dpsClosure.corrected << ',' << dpsClosure.residual << ','
           << dpsClosure.statDps << ',' << dpsClosure.statSpsMap << ',' << dpsClosure.statTotal << ','
           << dpsClosure.unsupportedBins << ',' << sps.emptyJpsiEvent << ',' << dps.emptyJpsiEvent << ','
           << sps.inconsistentJpsiVectors << ',' << dps.inconsistentJpsiVectors << ','
           << sps.malformedMuon << ',' << dps.malformedMuon << '\n';
}

}  // namespace

void single_jpsi_acceptance_closure(
    const char *outputDirectory,
    const char *spsInputList,
    const char *dpsInputList) {
    cout << setprecision(17);
    const SampleCounts sps = countSample("SPS", spsInputList);
    const SampleCounts dps = countSample("DPS", dpsInputList);
    const ClosureResult spsClosure = calculateSpsSelfClosure(sps);
    const ClosureResult dpsClosure = calculateDpsClosure(sps, dps);
    const string output = outputDirectory;
    writeAcceptanceTable(output + "/acceptance_sps_single_jpsi.txt", sps);
    writeAcceptanceTable(output + "/acceptance_dps_single_jpsi.txt", dps);
    writeBinComparison(output + "/bin_comparison.csv", sps, dps);
    writeSummary(output + "/summary.csv", sps, dps, spsClosure, dpsClosure);

    cout << "sps_single_fiducial=" << sps.totalFiducial()
         << "\nsps_single_accepted=" << sps.totalAccepted()
         << "\ndps_single_fiducial=" << dps.totalFiducial()
         << "\ndps_single_accepted=" << dps.totalAccepted()
         << "\nsps_raw_acceptance="
         << static_cast<double>(sps.totalAccepted()) / sps.totalFiducial()
         << "\ndps_raw_acceptance="
         << static_cast<double>(dps.totalAccepted()) / dps.totalFiducial()
         << "\nsps_corrected_single=" << spsClosure.corrected
         << "\nsps_closure=" << spsClosure.residual
         << "\ndps_corrected_single=" << dpsClosure.corrected
         << "\ndps_closure=" << dpsClosure.residual
         << "\ndps_closure_stat_total=" << dpsClosure.statTotal << endl;
}
