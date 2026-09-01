#include "TFile.h"
#include "TNamed.h"
#include "TSystem.h"
#include "TTree.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

void prepare_dps_template_input(
    const char *input =
        "Data_driven/results/dps_reference/"
        "WeightDPS_current.root",
    const char *output =
        "Data_driven/results/dps_reference/"
        "nominal_dps_template_input.root") {
  TFile source(input, "READ");
  auto *data = dynamic_cast<TTree *>(source.Get("data"));
  if (!data) {
    std::cerr << "Missing data tree in " << input << std::endl;
    return;
  }

  const char *required[] = {"Jpsi_pt1", "Jpsi_y1", "Jpsi_pt2", "Jpsi_y2",
                            "evt_mass", "evt_y", "delta_y", "delta_phi",
                            "evt_weight"};
  for (const char *name : required) {
    if (!data->GetBranch(name)) {
      std::cerr << "Missing branch " << name << " in " << input << std::endl;
      return;
    }
  }

  double pt1 = 0, y1 = 0, pt2 = 0, y2 = 0, mass = 0, pairY = 0, dy = 0,
         dphi = 0, weight = 0;
  data->SetBranchAddress("Jpsi_pt1", &pt1);
  data->SetBranchAddress("Jpsi_y1", &y1);
  data->SetBranchAddress("Jpsi_pt2", &pt2);
  data->SetBranchAddress("Jpsi_y2", &y2);
  data->SetBranchAddress("evt_mass", &mass);
  data->SetBranchAddress("evt_y", &pairY);
  data->SetBranchAddress("delta_y", &dy);
  data->SetBranchAddress("delta_phi", &dphi);
  data->SetBranchAddress("evt_weight", &weight);

  TString outputPath(output);
  gSystem->mkdir(gSystem->DirName(outputPath), true);
  TFile target(output, "RECREATE");
  TTree tree("mix", "nominal DPS MC adapted to the template input interface");
  ULong64_t run1 = 0, lumi1 = 0, event1 = 0, run2 = 0, lumi2 = 0, event2 = 0;
  tree.Branch("run1", &run1);
  tree.Branch("lumi1", &lumi1);
  tree.Branch("event1", &event1);
  tree.Branch("run2", &run2);
  tree.Branch("lumi2", &lumi2);
  tree.Branch("event2", &event2);
  tree.Branch("Jpsi_pt1", &pt1);
  tree.Branch("Jpsi_y1", &y1);
  tree.Branch("Jpsi_pt2", &pt2);
  tree.Branch("Jpsi_y2", &y2);
  tree.Branch("evt_mass", &mass);
  tree.Branch("evt_y", &pairY);
  tree.Branch("delta_y", &dy);
  tree.Branch("delta_phi", &dphi);
  tree.Branch("evt_weight", &weight);

  Long64_t accepted = 0, rejected = 0;
  double sumWeights = 0, sumWeights2 = 0;
  for (Long64_t entry = 0; entry < data->GetEntries(); ++entry) {
    data->GetEntry(entry);
    if (!(weight > 0) || !std::isfinite(weight) || mass < 7.5 ||
        dy < 0 || dy > 4 || dphi < 0 || dphi > 3.14159265358979323846) {
      ++rejected;
      continue;
    }
    event1 = static_cast<ULong64_t>(2 * entry);
    event2 = static_cast<ULong64_t>(2 * entry + 1);
    tree.Fill();
    ++accepted;
    sumWeights += weight;
    sumWeights2 += weight * weight;
  }

  std::ostringstream text;
  text << "source=" << input << "\n"
       << "source_tree=data\n"
       << "output_tree=mix\n"
       << "adapter_only=true\n"
       << "synthetic_source_ids=true\n"
       << "source_entries=" << data->GetEntries() << "\n"
       << "accepted_entries=" << accepted << "\n"
       << "rejected_entries=" << rejected << "\n"
       << "sum_weights=" << sumWeights << "\n"
       << "sum_weights2=" << sumWeights2 << "\n";
  TNamed metadata("run_metadata", text.str().c_str());
  tree.Write();
  metadata.Write();
  target.Close();
  source.Close();
  std::cout << text.str() << "output=" << output << std::endl;
}
