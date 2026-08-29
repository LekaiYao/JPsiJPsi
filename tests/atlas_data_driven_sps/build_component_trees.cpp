#include "TFile.h"
#include "TNamed.h"
#include "TSystem.h"
#include "TTree.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
namespace {
constexpr double pi = 3.14159265358979323846;
struct Cell {
  int id, iy, ip;
  double yl, yh, pl, ph, mix, dps, wd, ws;
  bool has(double y, double p) const {
    return y >= yl && (y < yh || (yh == 4 && y <= yh)) && p >= pl &&
           (p < ph || (ph >= pi - 1e-12 && p <= ph));
  }
};
std::vector<Cell> load(TFile &f) {
  auto *t = (TTree *)f.Get("cells");
  int iy, ip;
  double yl, yh, pl, ph, mix, dps, wd, ws;
  t->SetBranchAddress("iy", &iy);
  t->SetBranchAddress("iphi", &ip);
  t->SetBranchAddress("dy_low", &yl);
  t->SetBranchAddress("dy_high", &yh);
  t->SetBranchAddress("dphi_low", &pl);
  t->SetBranchAddress("dphi_high", &ph);
  t->SetBranchAddress("mix", &mix);
  t->SetBranchAddress("dps", &dps);
  t->SetBranchAddress("weight_dps", &wd);
  t->SetBranchAddress("weight_sps", &ws);
  std::vector<Cell> v;
  for (Long64_t i = 0; i < t->GetEntries(); i++) {
    t->GetEntry(i);
    v.push_back({int(i), iy, ip, yl, yh, pl, ph, mix, dps, wd, ws});
  }
  t->ResetBranchAddresses();
  return v;
}
int find(const std::vector<Cell> &v, double y, double p) {
  for (auto &c : v)
    if (c.has(y, p))
      return c.id;
  return -1;
}
} // namespace
void build_component_trees(
    const char *outdir =
        "tests/atlas_data_driven_sps/results/component_trees_adaptive14",
    const char *templateInput =
        "tests/atlas_data_driven_sps/results/templates_2d_adaptive14/"
        "templates_2d.root",
    const char *dataInput =
        "Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/WeightData.root",
    const char *mixedInput =
        "tests/atlas_data_driven_sps/results/mixed_dps_5M_seed20260728.root") {
  gSystem->mkdir(outdir, true);
  std::string base(outdir);
  TFile fw(templateInput);
  auto cells = load(fw);
  if (cells.empty()) {
    std::cerr << "empty cell table\n";
    return;
  }
  double sumMixCells = 0, expectedDps = 0;
  for (const auto &cell : cells) {
    sumMixCells += cell.mix;
    expectedDps += cell.dps;
  }
  const double alpha = sumMixCells ? expectedDps / sumMixCells : 0;
  TFile fi(dataInput);
  auto *ti = (TTree *)fi.Get("data");
  const Long64_t dataInputEntries = ti->GetEntries();
  double dy, dp, bw;
  ti->SetBranchAddress("delta_y", &dy);
  ti->SetBranchAddress("delta_phi", &dp);
  ti->SetBranchAddress("evt_weight", &bw);
  TFile fd((base + "/component_data.root").c_str(), "RECREATE");
  auto *td = ti->CloneTree(0);
  td->SetName("data_components");
  int cell_id;
  double cell_ws, cell_wd, weight_sps, weight_dps;
  td->Branch("component_cell", &cell_id);
  td->Branch("component_fraction_sps", &cell_ws);
  td->Branch("component_fraction_dps", &cell_wd);
  td->Branch("evt_weight_sps", &weight_sps);
  td->Branch("evt_weight_dps", &weight_dps);
  Long64_t dataMiss = 0;
  double dataBase = 0, dataSps = 0, dataDps = 0, maxDataClosure = 0;
  for (Long64_t i = 0; i < dataInputEntries; i++) {
    ti->GetEntry(i);
    cell_id = find(cells, dy, dp);
    if (cell_id < 0) {
      dataMiss++;
      cell_ws = cell_wd = weight_sps = weight_dps = 0;
    } else {
      cell_ws = cells[cell_id].ws;
      cell_wd = cells[cell_id].wd;
      weight_sps = bw * cell_ws;
      weight_dps = bw * cell_wd;
      maxDataClosure =
          std::max(maxDataClosure, std::abs(weight_sps + weight_dps - bw));
    }
    dataBase += bw;
    dataSps += weight_sps;
    dataDps += weight_dps;
    td->Fill();
  }
  fw.cd();
  auto *cellTree = (TTree *)fw.Get("cells");
  fd.cd();
  cellTree->CloneTree()->Write("component_cells");
  std::string dataMetaText =
      std::string("fractions_from=") + templateInput +
      "; data_source=" + dataInput +
      "; do_not_interpret_weighted_raw_entries_as_PP_yields";
  TNamed dataMeta("component_metadata", dataMetaText.c_str());
  dataMeta.Write();
  const Long64_t dataOutputEntries = td->GetEntries();
  td->Write();
  fd.Close();
  fi.Close();
  TFile fm(mixedInput);
  auto *tm = (TTree *)fm.Get("mix");
  const Long64_t mixInputEntries = tm->GetEntries();
  ULong64_t r1, l1, e1, r2, l2, e2;
  double pt1, y1, pt2, y2, mass, ey, mw;
  tm->SetBranchAddress("run1", &r1);
  tm->SetBranchAddress("lumi1", &l1);
  tm->SetBranchAddress("event1", &e1);
  tm->SetBranchAddress("run2", &r2);
  tm->SetBranchAddress("lumi2", &l2);
  tm->SetBranchAddress("event2", &e2);
  tm->SetBranchAddress("Jpsi_pt1", &pt1);
  tm->SetBranchAddress("Jpsi_y1", &y1);
  tm->SetBranchAddress("Jpsi_pt2", &pt2);
  tm->SetBranchAddress("Jpsi_y2", &y2);
  tm->SetBranchAddress("evt_mass", &mass);
  tm->SetBranchAddress("evt_y", &ey);
  tm->SetBranchAddress("delta_y", &dy);
  tm->SetBranchAddress("delta_phi", &dp);
  tm->SetBranchAddress("evt_weight", &mw);
  TFile fo((base + "/component_mix.root").c_str(), "RECREATE");
  TTree mo("mix_components", "normalized mixed-event DPS components");
  double nw;
  mo.Branch("run1", &r1);
  mo.Branch("lumi1", &l1);
  mo.Branch("event1", &e1);
  mo.Branch("run2", &r2);
  mo.Branch("lumi2", &l2);
  mo.Branch("event2", &e2);
  mo.Branch("Jpsi_pt1", &pt1);
  mo.Branch("Jpsi_y1", &y1);
  mo.Branch("Jpsi_pt2", &pt2);
  mo.Branch("Jpsi_y2", &y2);
  mo.Branch("evt_mass", &mass);
  mo.Branch("evt_y", &ey);
  mo.Branch("delta_y", &dy);
  mo.Branch("delta_phi", &dp);
  mo.Branch("evt_weight", &mw);
  mo.Branch("component_cell", &cell_id);
  mo.Branch("evt_weight_dps_normalized", &nw);
  Long64_t mixMiss = 0;
  double mixRaw = 0, mixNorm = 0;
  std::vector<double> cellNorm(cells.size());
  for (Long64_t i = 0; i < mixInputEntries; i++) {
    tm->GetEntry(i);
    cell_id = find(cells, dy, dp);
    nw = alpha * mw;
    if (cell_id < 0) {
      mixMiss++;
      nw = 0;
    } else
      cellNorm[cell_id] += nw;
    mixRaw += mw;
    mixNorm += nw;
    mo.Fill();
  }
  fw.cd();
  cellTree = (TTree *)fw.Get("cells");
  fo.cd();
  cellTree->CloneTree()->Write("component_cells");
  std::ostringstream mixMetaText;
  mixMetaText << std::setprecision(12) << "cells=" << cells.size()
              << "; alpha=" << alpha
              << "; source=" << mixedInput;
  TNamed mixMeta("component_metadata", mixMetaText.str().c_str());
  mixMeta.Write();
  const Long64_t mixOutputEntries = mo.GetEntries();
  mo.Write();
  fo.Close();
  fm.Close();
  double maxCell = 0;
  fw.cd();
  cellTree = (TTree *)fw.Get("cells");
  for (Long64_t i = 0; i < cellTree->GetEntries(); i++) {
    maxCell = std::max(maxCell, std::abs(cellNorm[i] - cells[i].dps));
  }
  std::ofstream s(base + "/summary.txt");
  s << std::setprecision(12)
    << "cells=" << cells.size() << "\ncell_source=" << templateInput
    << "\ndata_source=" << dataInput << "\nmixed_source=" << mixedInput
    << "\nalpha="
    << alpha << "\ndata_input_entries=" << dataInputEntries
    << "\ndata_output_entries=" << dataOutputEntries
    << "\ndata_unassigned_entries=" << dataMiss
    << "\ndata_sum_evt_weight=" << dataBase
    << "\ndata_sum_evt_weight_sps=" << dataSps
    << "\ndata_sum_evt_weight_dps=" << dataDps
    << "\ndata_max_abs_weight_closure=" << maxDataClosure
    << "\nmix_input_entries=" << mixInputEntries
    << "\nmix_output_entries=" << mixOutputEntries
    << "\nmix_unassigned_entries=" << mixMiss
    << "\nmix_sum_raw_weight=" << mixRaw
    << "\nmix_sum_normalized_weight=" << mixNorm
    << "\nexpected_dps_template_integral=" << expectedDps
    << "\nmix_integral_minus_expected=" << mixNorm - expectedDps
    << "\nmax_abs_mixed_cell_minus_expected=" << maxCell
    << "\nwarning=data_component_weighted_raw_sums_are_not_PP_yields_rerun_4D_"
       "fits_in_step4=true\n";
  std::cout << "data=" << dataOutputEntries << " missing=" << dataMiss
            << " closure=" << maxDataClosure << "\nmix=" << mixOutputEntries
            << " missing=" << mixMiss << " normalized=" << mixNorm
            << " expected=" << expectedDps << " max cell diff=" << maxCell
            << "\noutput=" << outdir << std::endl;
}
