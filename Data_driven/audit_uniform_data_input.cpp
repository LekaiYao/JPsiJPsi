#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"

void audit_uniform_data_input(
    const char *newPath =
        "Data_driven/results/unified_input/WeightData_mJJ7p5.root",
    const char *outPath =
        "Data_driven/results/unified_input/selection_audit.txt",
    const char *oldPath =
        "Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/WeightData.root") {
    TFile newFile(newPath, "READ");
    TTree *newTree = static_cast<TTree *>(newFile.Get("data"));
    if (!newTree) {
        std::cerr << "Missing new input tree" << std::endl;
        return;
    }
    TFile *oldFile = nullptr;
    TTree *oldTree = nullptr;
    if (oldPath && std::string(oldPath).size()) {
        oldFile = TFile::Open(oldPath, "READ");
        if (oldFile && !oldFile->IsZombie())
            oldTree = static_cast<TTree *>(oldFile->Get("data"));
    }

    const std::vector<std::string> common = {
        "Jpsi_mass1", "Jpsi_ctau1", "Jpsi_pt1", "Jpsi_y1",
        "Jpsi_mass2", "Jpsi_ctau2", "Jpsi_pt2", "Jpsi_y2",
        "evt_vtxProb", "evt_weight", "evt_mass", "evt_mass2",
        "evt_y", "evt_pt", "delta_y", "delta_phi"};
    std::vector<double> oldValue(common.size(), 0.0);
    std::vector<double> newValue(common.size(), 0.0);
    std::vector<Long64_t> changed(common.size(), 0);
    for (std::size_t i = 0; i < common.size(); ++i) {
        if (!newTree->GetBranch(common[i].c_str())) {
            std::cerr << "Missing new branch " << common[i] << std::endl;
            if (oldFile) delete oldFile;
            return;
        }
        newTree->SetBranchAddress(common[i].c_str(), &newValue[i]);
        if (oldTree && oldTree->GetBranch(common[i].c_str()))
            oldTree->SetBranchAddress(common[i].c_str(), &oldValue[i]);
        else
            oldTree = nullptr;
    }

    Long64_t changedAny = 0;
    Long64_t changed4D = 0;
    Long64_t oldBadMass = 0;
    Long64_t badMass = 0;
    Long64_t badPt = 0;
    Long64_t badWeight = 0;
    Long64_t massDifferenceAboveMeV = 0;
    Long64_t weightDifferenceAbovePermille = 0;
    double maximumMassDifference = 0.0;
    double maximumRelativeWeightDifference = 0.0;
    const Long64_t entries = newTree->GetEntries();
    const Long64_t comparedEntries = oldTree
        ? std::min(oldTree->GetEntries(), newTree->GetEntries()) : 0;
    for (Long64_t entry = 0; entry < entries; ++entry) {
        newTree->GetEntry(entry);
        if (entry < comparedEntries) {
            oldTree->GetEntry(entry);
            bool any = false;
            bool fourD = false;
            for (std::size_t i = 0; i < common.size(); ++i) {
                const double scale = std::max(
                    1.0, std::max(std::abs(oldValue[i]), std::abs(newValue[i])));
                if (std::abs(oldValue[i] - newValue[i]) > 1e-10 * scale) {
                    ++changed[i];
                    any = true;
                    if (i == 0 || i == 1 || i == 4 || i == 5)
                        fourD = true;
                }
            }
            changedAny += any;
            changed4D += fourD;
            oldBadMass += oldValue[10] < 7.5;
            const double massDifference = std::abs(oldValue[10] - newValue[10]);
            massDifferenceAboveMeV += massDifference > 0.001;
            maximumMassDifference = std::max(maximumMassDifference, massDifference);
            const double denominator = std::max(std::abs(oldValue[9]), 1e-300);
            const double relativeWeightDifference =
                std::abs(oldValue[9] - newValue[9]) / denominator;
            weightDifferenceAbovePermille += relativeWeightDifference > 0.001;
            maximumRelativeWeightDifference =
                std::max(maximumRelativeWeightDifference, relativeWeightDifference);
        }
        badMass += newValue[10] < 7.5;
        badPt += newValue[2] < 10.0 || newValue[2] > 40.0 ||
                 newValue[6] < 10.0 || newValue[6] > 40.0;
        badWeight += !std::isfinite(newValue[9]) || newValue[9] <= 0.0;
    }

    std::ofstream out(outPath);
    out << std::setprecision(12)
        << "source_files=119\n"
        << "source_entries=6132\n"
        << "selection=passHLT && matchTrg && 10<=pt(Jpsi1,Jpsi2)<=40"
           " && samePV && m(JJ)>=7.5\n"
        << "reference_comparison=" << (oldTree ? "available" : "not_available") << "\n"
        << "reference_entries=" << (oldTree ? oldTree->GetEntries() : 0) << "\n"
        << "compared_entries=" << comparedEntries << "\n"
        << "new_entries=" << newTree->GetEntries() << "\n"
        << "old_bad_mass=" << oldBadMass << "\n"
        << "new_bad_mass=" << badMass << "\n"
        << "new_bad_pt=" << badPt << "\n"
        << "new_bad_weight=" << badWeight << "\n"
        << "rows_changed_any_common_branch=" << changedAny << "\n"
        << "rows_changed_4d_observable=" << changed4D << "\n";
    out << "rows_mass_difference_above_1MeV=" << massDifferenceAboveMeV << "\n"
        << "maximum_mass_difference_GeV=" << maximumMassDifference << "\n"
        << "rows_weight_difference_above_0p1percent="
        << weightDifferenceAbovePermille << "\n"
        << "maximum_relative_weight_difference="
        << maximumRelativeWeightDifference << "\n";
    for (std::size_t i = 0; i < common.size(); ++i)
        out << "rows_changed_" << common[i] << '=' << changed[i] << "\n";
    out.close();
    if (oldFile) delete oldFile;

    std::cout << "Saved " << outPath << std::endl;
}
