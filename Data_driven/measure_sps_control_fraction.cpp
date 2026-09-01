#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"

void measure_sps_control_fraction(
    const char *output,
    double controlDyMin = 1.8,
    double controlPhiMax = 1.5707963267948966,
    const char *input =
        "Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/WeightSPSstar.root") {
  TFile file(input, "READ");
  TTree *tree = dynamic_cast<TTree *>(file.Get("data"));
  if (!tree) {
    std::cerr << "Missing data tree in " << input << std::endl;
    return;
  }
  double dy = 0, phi = 0, weight = 0;
  tree->SetBranchAddress("delta_y", &dy);
  tree->SetBranchAddress("delta_phi", &phi);
  tree->SetBranchAddress("evt_weight", &weight);
  Long64_t accepted = 0, controlEntries = 0;
  double total = 0, total2 = 0, control = 0, control2 = 0;
  for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    if (!std::isfinite(weight) || !(weight > 0) || dy < 0 || dy > 4 ||
        phi < 0 || phi > 3.14159265358979323846)
      continue;
    ++accepted;
    total += weight;
    total2 += weight * weight;
    if (dy >= controlDyMin && phi <= controlPhiMax) {
      ++controlEntries;
      control += weight;
      control2 += weight * weight;
    }
  }
  TString path(output);
  gSystem->mkdir(gSystem->DirName(path), true);
  std::ofstream out(output);
  out << std::setprecision(12)
      << "input=" << input << "\n"
      << "tree_entries=" << tree->GetEntries() << "\n"
      << "accepted_entries=" << accepted << "\n"
      << "control_entries=" << controlEntries << "\n"
      << "control_dy_min=" << controlDyMin << "\n"
      << "control_phi_max=" << controlPhiMax << "\n"
      << "total_sum_weights=" << total << "\n"
      << "total_sum_weights2=" << total2 << "\n"
      << "control_sum_weights=" << control << "\n"
      << "control_sum_weights2=" << control2 << "\n"
      << "sps_control_fraction=" << control / total << "\n";
  std::cout << "sps_control_fraction=" << control / total << std::endl;
}
