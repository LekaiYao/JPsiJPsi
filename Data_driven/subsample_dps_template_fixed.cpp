#include "TFile.h"
#include "TNamed.h"
#include "TSystem.h"
#include "TTree.h"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Priority {
  std::uint64_t key = 0;
  Long64_t entry = -1;
};

void requireBranch(TTree *tree, const char *name) {
  if (!tree->GetBranch(name))
    throw std::runtime_error("Missing branch " + std::string(name));
}

}  // namespace

void subsample_dps_template_fixed(
    const char *input,
    const char *output,
    const char *summaryPath,
    const char *selectedEntriesPath,
    ULong64_t seed = 20260902) {
  try {
    TFile source(input, "READ");
    if (source.IsZombie())
      throw std::runtime_error("Cannot open input " + std::string(input));
    auto *tree = dynamic_cast<TTree *>(source.Get("mix"));
    if (!tree)
      throw std::runtime_error("Missing mix tree in " + std::string(input));
    requireBranch(tree, "evt_weight");
    requireBranch(tree, "delta_y");
    requireBranch(tree, "delta_phi");

    const Long64_t inputEntries = tree->GetEntries();
    if (inputEntries < 2)
      throw std::runtime_error("Need at least two DPS MC entries");
    const Long64_t selectedEntries = inputEntries / 2;
    const std::string sampling =
        inputEntries % 2 == 0
            ? "exact_half_without_replacement_by_smallest_mt19937_64_priorities"
            : "floor_half_without_replacement_by_smallest_mt19937_64_priorities";

    std::mt19937_64 generator(seed);
    std::vector<Priority> priorities;
    priorities.reserve(inputEntries);
    for (Long64_t entry = 0; entry < inputEntries; ++entry)
      priorities.push_back({generator(), entry});
    std::sort(priorities.begin(), priorities.end(),
              [](const Priority &a, const Priority &b) {
                return a.key < b.key ||
                       (a.key == b.key && a.entry < b.entry);
              });

    std::vector<unsigned char> selected(inputEntries, 0);
    std::vector<std::uint64_t> keys(inputEntries, 0);
    Long64_t keyCollisions = 0;
    for (Long64_t i = 1; i < inputEntries; ++i)
      if (priorities[i - 1].key == priorities[i].key)
        ++keyCollisions;
    for (Long64_t i = 0; i < selectedEntries; ++i) {
      selected[priorities[i].entry] = 1;
      keys[priorities[i].entry] = priorities[i].key;
    }

    double weight = 0, deltaY = 0, deltaPhi = 0;
    tree->SetBranchAddress("evt_weight", &weight);
    tree->SetBranchAddress("delta_y", &deltaY);
    tree->SetBranchAddress("delta_phi", &deltaPhi);

    TString outputName(output);
    gSystem->mkdir(gSystem->DirName(outputName), true);
    TFile target(output, "RECREATE");
    if (target.IsZombie())
      throw std::runtime_error("Cannot create output " + std::string(output));
    auto *half = tree->CloneTree(0);
    half->SetName("mix");
    half->SetTitle("fixed-seed half-sample DPS MC template");
    tree->CopyAddresses(half);

    TTree selection("selected_entries",
                    "source entries selected by fixed random priority");
    Long64_t selectedEntry = -1;
    ULong64_t randomPriority = 0;
    selection.Branch("input_entry", &selectedEntry);
    selection.Branch("random_priority", &randomPriority);

    std::ofstream selectedCsv(selectedEntriesPath);
    if (!selectedCsv)
      throw std::runtime_error("Cannot create selected-entry CSV");
    selectedCsv << "input_entry,random_priority\n";

    double inputSumW = 0, selectedSumW = 0;
    double inputCrSumW = 0, selectedCrSumW = 0;
    Long64_t selectedCount = 0;
    for (Long64_t entry = 0; entry < inputEntries; ++entry) {
      tree->GetEntry(entry);
      inputSumW += weight;
      if (deltaY >= 1.8 && deltaPhi <= 1.5707963267948966)
        inputCrSumW += weight;
      if (!selected[entry])
        continue;
      half->Fill();
      selectedEntry = entry;
      randomPriority = keys[entry];
      selection.Fill();
      selectedCsv << entry << "," << randomPriority << "\n";
      selectedSumW += weight;
      if (deltaY >= 1.8 && deltaPhi <= 1.5707963267948966)
        selectedCrSumW += weight;
      ++selectedCount;
    }
    selectedCsv.close();

    if (selectedCount != selectedEntries ||
        half->GetEntries() != selectedEntries ||
        selection.GetEntries() != selectedEntries)
      throw std::runtime_error("Selected-entry count mismatch");

    std::ostringstream metadataText;
    metadataText << std::setprecision(12)
                 << "source=" << input << "\n"
                 << "sampling=" << sampling << "\n"
                 << "half_count_policy=floor_integer_division\n"
                 << "engine=std_mt19937_64\n"
                 << "seed=" << seed << "\n"
                 << "input_entries=" << inputEntries << "\n"
                 << "selected_entries=" << selectedEntries << "\n"
                 << "selected_fraction="
                 << static_cast<double>(selectedEntries) / inputEntries << "\n"
                 << "priority_collisions=" << keyCollisions << "\n"
                 << "output_order=source_entry_order\n";
    TNamed metadata("subsample_metadata", metadataText.str().c_str());
    target.cd();
    half->Write();
    selection.Write();
    metadata.Write();
    target.Close();
    source.Close();

    std::ofstream summary(summaryPath);
    if (!summary)
      throw std::runtime_error("Cannot create sampling summary");
    summary << std::setprecision(12)
            << "status=complete\n"
            << "source=" << input << "\n"
            << "output=" << output << "\n"
            << "selected_entries_csv=" << selectedEntriesPath << "\n"
            << "sampling=" << sampling << "\n"
            << "half_count_policy=floor_integer_division\n"
            << "engine=std_mt19937_64\n"
            << "seed=" << seed << "\n"
            << "input_entries=" << inputEntries << "\n"
            << "selected_entries=" << selectedEntries << "\n"
            << "selected_fraction="
            << static_cast<double>(selectedEntries) / inputEntries << "\n"
            << "priority_collisions=" << keyCollisions << "\n"
            << "input_sum_weights=" << inputSumW << "\n"
            << "selected_sum_weights=" << selectedSumW << "\n"
            << "selected_weight_fraction=" << selectedSumW / inputSumW << "\n"
            << "input_nominal_cr_sum_weights=" << inputCrSumW << "\n"
            << "selected_nominal_cr_sum_weights=" << selectedCrSumW << "\n"
            << "input_nominal_cr_fraction=" << inputCrSumW / inputSumW << "\n"
            << "selected_nominal_cr_fraction="
            << selectedCrSumW / selectedSumW << "\n";
    summary.close();
    std::cout << "DPS_FIXED_HALF_SAMPLE_COMPLETE input_entries="
              << inputEntries << " selected_entries=" << selectedEntries
              << " seed=" << seed << std::endl;
  } catch (const std::exception &error) {
    std::cerr << "DPS_FIXED_HALF_SAMPLE_FAILED " << error.what() << std::endl;
    gSystem->Exit(1);
  }
}
