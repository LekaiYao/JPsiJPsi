#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "TFile.h"
#include "TLorentzVector.h"
#include "TNamed.h"
#include "TSystem.h"
#include "TTree.h"

#include "build_mixed_dps.cpp"

namespace {

struct SelectedPair {
  Long64_t selectedEntry = -1;
  Long64_t sourceEntry = -1;
  Jpsi first;
  Jpsi second;
};

bool closeEnough(double a, double b) {
  return std::abs(a - b) <= 1e-8 * std::max({1.0, std::abs(a), std::abs(b)});
}

std::vector<Long64_t> readSelectedEntries(const char *path) {
  std::ifstream input(path);
  if (!input)
    throw std::runtime_error("Cannot open selected-entry CSV " +
                             std::string(path));
  std::string line;
  std::getline(input, line);
  if (line != "input_entry,random_priority")
    throw std::runtime_error("Unexpected selected-entry CSV header");
  std::vector<Long64_t> entries;
  while (std::getline(input, line)) {
    if (line.empty())
      continue;
    const std::size_t comma = line.find(',');
    if (comma == std::string::npos)
      throw std::runtime_error("Malformed selected-entry CSV row");
    entries.push_back(std::stoll(line.substr(0, comma)));
  }
  if (entries.empty())
    throw std::runtime_error("Selected-entry CSV is empty");
  for (std::size_t i = 0; i < entries.size(); ++i) {
    if (entries[i] < 0 || (i && entries[i] <= entries[i - 1]))
      throw std::runtime_error("Selected entries are not strictly increasing");
  }
  return entries;
}

std::vector<SelectedPair> loadSelectedPairs(const char *candidatePath,
                                            const std::vector<Long64_t> &entries,
                                            Long64_t &candidateEntries) {
  TFile input(candidatePath, "READ");
  if (input.IsZombie())
    throw std::runtime_error("Cannot open candidate file " +
                             std::string(candidatePath));
  auto *first = dynamic_cast<TTree *>(input.Get("jpsi1"));
  auto *second = dynamic_cast<TTree *>(input.Get("jpsi2"));
  if (!first || !second || first->GetEntries() != second->GetEntries())
    throw std::runtime_error("Invalid jpsi1/jpsi2 candidate trees");
  candidateEntries = first->GetEntries();
  if (entries.back() >= candidateEntries)
    throw std::runtime_error("Selected entry exceeds candidate tree");

  Long64_t source1 = -1, source2 = -1;
  double pt1 = 0, eta1 = 0, y1 = 0, phi1 = 0, mass1 = 0;
  double pt2 = 0, eta2 = 0, y2 = 0, phi2 = 0, mass2 = 0;
  first->SetBranchAddress("source_entry", &source1);
  first->SetBranchAddress("Jpsi_pt1", &pt1);
  first->SetBranchAddress("Jpsi_eta1", &eta1);
  first->SetBranchAddress("Jpsi_y1", &y1);
  first->SetBranchAddress("Jpsi_phi1", &phi1);
  first->SetBranchAddress("Jpsi_mass1", &mass1);
  second->SetBranchAddress("source_entry", &source2);
  second->SetBranchAddress("Jpsi_pt2", &pt2);
  second->SetBranchAddress("Jpsi_eta2", &eta2);
  second->SetBranchAddress("Jpsi_y2", &y2);
  second->SetBranchAddress("Jpsi_phi2", &phi2);
  second->SetBranchAddress("Jpsi_mass2", &mass2);

  std::vector<SelectedPair> selected;
  selected.reserve(entries.size());
  for (Long64_t entry : entries) {
    first->GetEntry(entry);
    second->GetEntry(entry);
    if (source1 != source2)
      throw std::runtime_error("Slot source-entry mismatch");
    selected.push_back(
        {entry, source1, {pt1, eta1, y1, phi1, mass1},
         {pt2, eta2, y2, phi2, mass2}});
  }
  input.Close();
  return selected;
}

Long64_t validateAgainstDirectHalf(const char *directPath,
                                   const std::vector<SelectedPair> &selected) {
  TFile input(directPath, "READ");
  if (input.IsZombie())
    throw std::runtime_error("Cannot open direct half-MC file " +
                             std::string(directPath));
  auto *tree = dynamic_cast<TTree *>(input.Get("mix"));
  if (!tree || tree->GetEntries() != static_cast<Long64_t>(selected.size()))
    throw std::runtime_error("Direct half-MC entry count mismatch");
  double pt1 = 0, y1 = 0, pt2 = 0, y2 = 0, mass = 0, pairY = 0;
  double dy = 0, dphi = 0;
  tree->SetBranchAddress("Jpsi_pt1", &pt1);
  tree->SetBranchAddress("Jpsi_y1", &y1);
  tree->SetBranchAddress("Jpsi_pt2", &pt2);
  tree->SetBranchAddress("Jpsi_y2", &y2);
  tree->SetBranchAddress("evt_mass", &mass);
  tree->SetBranchAddress("evt_y", &pairY);
  tree->SetBranchAddress("delta_y", &dy);
  tree->SetBranchAddress("delta_phi", &dphi);

  Long64_t mismatches = 0;
  for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    const auto &pair = selected[i];
    TLorentzVector a, b;
    a.SetPtEtaPhiM(pair.first.pt, pair.first.eta, pair.first.phi,
                   pair.first.mass);
    b.SetPtEtaPhiM(pair.second.pt, pair.second.eta, pair.second.phi,
                   pair.second.mass);
    const TLorentzVector combined = a + b;
    const double expectedDy = std::abs(pair.first.y - pair.second.y);
    const double expectedDphi =
        kPi - std::abs(std::abs(pair.first.phi - pair.second.phi) - kPi);
    if (!closeEnough(pt1, pair.first.pt) || !closeEnough(y1, pair.first.y) ||
        !closeEnough(pt2, pair.second.pt) || !closeEnough(y2, pair.second.y) ||
        !closeEnough(mass, combined.M()) ||
        !closeEnough(pairY, std::abs(combined.Rapidity())) ||
        !closeEnough(dy, expectedDy) || !closeEnough(dphi, expectedDphi))
      ++mismatches;
  }
  input.Close();
  return mismatches;
}

} // namespace

void mix_dps_mc_half_crossslot(const char *candidatePath,
                               const char *selectedEntriesPath,
                               const char *directHalfPath,
                               const char *outputPath,
                               const char *summaryPath,
                               const char *acceptancePath = kAccFile,
                               const char *efficiencyPath = kEffFile) {
  try {
    const std::vector<Long64_t> entries =
        readSelectedEntries(selectedEntriesPath);
    Long64_t candidateEntries = 0;
    const std::vector<SelectedPair> selected =
        loadSelectedPairs(candidatePath, entries, candidateEntries);
    const Long64_t validationMismatches =
        validateAgainstDirectHalf(directHalfPath, selected);
    if (validationMismatches != 0)
      throw std::runtime_error("Selected candidates do not match direct half MC");

    Correction correction;
    if (!correction.load(acceptancePath, efficiencyPath))
      throw std::runtime_error("Cannot load frozen correction tables");

    TString output(outputPath);
    gSystem->mkdir(gSystem->DirName(output), true);
    TFile target(outputPath, "RECREATE");
    if (target.IsZombie())
      throw std::runtime_error("Cannot create output " +
                               std::string(outputPath));
    TTree tree("mix", "deterministic fixed-half DPS-MC cross-slot mixing");
    ULong64_t run1 = 0, lumi1 = 0, event1 = 0;
    ULong64_t run2 = 0, lumi2 = 0, event2 = 0;
    Long64_t sourceEntry1 = -1, sourceEntry2 = -1;
    double pt1 = 0, y1 = 0, pt2 = 0, y2 = 0;
    double pairMass = 0, pairY = 0, dy = 0, dphi = 0, weight = 0;
    tree.Branch("run1", &run1);
    tree.Branch("lumi1", &lumi1);
    tree.Branch("event1", &event1);
    tree.Branch("run2", &run2);
    tree.Branch("lumi2", &lumi2);
    tree.Branch("event2", &event2);
    tree.Branch("source_entry1", &sourceEntry1);
    tree.Branch("source_entry2", &sourceEntry2);
    tree.Branch("Jpsi_pt1", &pt1);
    tree.Branch("Jpsi_y1", &y1);
    tree.Branch("Jpsi_pt2", &pt2);
    tree.Branch("Jpsi_y2", &y2);
    tree.Branch("evt_mass", &pairMass);
    tree.Branch("evt_y", &pairY);
    tree.Branch("delta_y", &dy);
    tree.Branch("delta_phi", &dphi);
    tree.Branch("evt_weight", &weight);

    Long64_t attempts = 0, sameSourceRejected = 0, massRejected = 0;
    Long64_t correctionRejected = 0, accepted = 0;
    double sumWeights = 0, sumWeights2 = 0;
    for (const auto &first : selected) {
      for (const auto &second : selected) {
        ++attempts;
        if (first.sourceEntry == second.sourceEntry) {
          ++sameSourceRejected;
          continue;
        }
        TLorentzVector a, b;
        a.SetPtEtaPhiM(first.first.pt, first.first.eta, first.first.phi,
                       first.first.mass);
        b.SetPtEtaPhiM(second.second.pt, second.second.eta,
                       second.second.phi, second.second.mass);
        const TLorentzVector combined = a + b;
        if (combined.M() < 7.5) {
          ++massRejected;
          continue;
        }
        weight = correction.weight(first.first, second.second);
        if (!(weight > 0) || !std::isfinite(weight)) {
          ++correctionRejected;
          continue;
        }
        sourceEntry1 = first.sourceEntry;
        sourceEntry2 = second.sourceEntry;
        event1 = static_cast<ULong64_t>(sourceEntry1);
        event2 = static_cast<ULong64_t>(sourceEntry2);
        pt1 = first.first.pt;
        y1 = first.first.y;
        pt2 = second.second.pt;
        y2 = second.second.y;
        pairMass = combined.M();
        pairY = std::abs(combined.Rapidity());
        dy = std::abs(y1 - y2);
        dphi = kPi -
               std::abs(std::abs(first.first.phi - second.second.phi) - kPi);
        tree.Fill();
        ++accepted;
        sumWeights += weight;
        sumWeights2 += weight * weight;
      }
    }

    std::ostringstream metadataText;
    metadataText << std::setprecision(12)
                 << "candidate_input=" << candidatePath << "\n"
                 << "selected_entries_input=" << selectedEntriesPath << "\n"
                 << "direct_half_input=" << directHalfPath << "\n"
                 << "mixing=deterministic_all_cross_Jpsi1_A_x_Jpsi2_B\n"
                 << "same_source_constraint=A_not_equal_B\n"
                 << "pair_weight=recomputed_frozen_acceptance_efficiency_"
                    "correction\n"
                 << "full_candidate_events=" << candidateEntries << "\n"
                 << "selected_pair_events=" << selected.size() << "\n"
                 << "slot1_candidates=" << selected.size() << "\n"
                 << "slot2_candidates=" << selected.size() << "\n"
                 << "direct_half_validation_mismatches="
                 << validationMismatches << "\n"
                 << "cartesian_attempts=" << attempts << "\n"
                 << "same_source_rejected=" << sameSourceRejected << "\n"
                 << "mass_rejected=" << massRejected << "\n"
                 << "correction_rejected=" << correctionRejected << "\n"
                 << "accepted_pairs=" << accepted << "\n"
                 << "sum_weights=" << sumWeights << "\n"
                 << "sum_weights2=" << sumWeights2 << "\n";
    TNamed metadata("run_metadata", metadataText.str().c_str());
    tree.Write();
    metadata.Write();
    target.Close();
    std::ofstream summary(summaryPath);
    if (!summary)
      throw std::runtime_error("Cannot create summary " +
                               std::string(summaryPath));
    summary << "status=complete\n" << metadataText.str()
            << "output=" << outputPath << "\n";
    summary.close();
    std::cout << "DPS_MC_HALF_CROSSSLOT_MIXING_COMPLETE\n"
              << metadataText.str() << "output=" << outputPath << std::endl;
  } catch (const std::exception &error) {
    std::cerr << "DPS_MC_HALF_CROSSSLOT_MIXING_FAILED " << error.what()
              << std::endl;
    gSystem->Exit(1);
  }
}
