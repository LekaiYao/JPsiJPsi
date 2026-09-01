#include "TCanvas.h"
#include "TFile.h"
#include "TH2Poly.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include <cstdio>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
namespace {
constexpr double pi = 3.14159265358979323846;
struct Cell {
  int iy, ip;
  double yl, yh, pl, ph, n, en, m = 0, vm = 0, d = 0, ed = 0, s = 0, es = 0,
                                wd = 0, ws = 0;
  bool cr(double dyMin, double phiMax) const {
    return yl >= dyMin - 1e-8 && ph <= phiMax + 1e-8;
  }
  bool has(double y, double p) const {
    return y >= yl && (y < yh || (yh == 4 && y <= yh)) && p >= pl &&
           (p < ph || (ph == pi && p <= ph));
  }
};
std::vector<std::string> split(const std::string &s) {
  std::vector<std::string> v;
  std::stringstream q(s);
  std::string x;
  while (std::getline(q, x, ","[0]))
    v.push_back(x);
  return v;
}
std::vector<Cell> read(const char *fitInput) {
  std::ifstream f(fitInput);
  std::string l;
  std::getline(f, l);
  std::vector<Cell> v;
  while (std::getline(f, l)) {
    auto x = split(l);
    if (x.size() < 19)
      continue;
    v.push_back({std::stoi(x[0]), std::stoi(x[1]), std::stod(x[2]),
                 std::stod(x[3]), std::stod(x[4]), std::stod(x[5]),
                 std::stod(x[11]), std::stod(x[12])});
  }
  return v;
}
int find(const std::vector<Cell> &v, double y, double p) {
  for (size_t i = 0; i < v.size(); i++)
    if (v[i].has(y, p))
      return i;
  return -1;
}
void draw(TH2Poly &h, const std::vector<Cell> &cells, const std::string &p,
          const char *l, const std::vector<double> &values,
          const std::vector<double> &errors, bool showErrors,
          const char *valueFormat) {
  TCanvas c("ca", "", 900, 700);
  c.SetRightMargin(.16);
  c.SetLeftMargin(.12);
  c.SetBottomMargin(.12);
  h.SetStats(0);
  h.SetTitle("");
  h.SetMarkerSize(1.0);
  h.GetXaxis()->SetTitle("|#Delta y|");
  h.GetYaxis()->SetTitle("|#Delta#phi|");
  h.Draw("COLZ");
  TLatex valueText;
  valueText.SetTextAlign(22);
  valueText.SetTextFont(42);
  valueText.SetTextSize(.023);
  char buffer[128], valueBuffer[48], errorBuffer[48];
  for (size_t i = 0; i < cells.size(); ++i) {
    if (showErrors) {
      std::snprintf(valueBuffer, sizeof(valueBuffer), valueFormat, values[i]);
      std::snprintf(errorBuffer, sizeof(errorBuffer), valueFormat, errors[i]);
      std::snprintf(buffer, sizeof(buffer), "%s #pm %s", valueBuffer,
                    errorBuffer);
    } else {
      std::snprintf(buffer, sizeof(buffer), valueFormat, values[i]);
    }
    valueText.DrawLatex(0.5 * (cells[i].yl + cells[i].yh),
                        0.5 * (cells[i].pl + cells[i].ph), buffer);
  }
  TLatex t;
  t.SetNDC();
  t.SetTextFont(42);
  t.SetTextSize(.034);
  t.DrawLatex(.12, .94, l);
  c.SaveAs(p.c_str());
}
void configure(TH2Poly &h, const char *n, const std::vector<Cell> &v) {
  h.SetName(n);
  for (const auto &c : v)
    h.AddBin(c.yl, c.pl, c.yh, c.ph);
}
double significance(const Cell &c) {
  return c.es > 0 ? c.s / c.es : 0.0;
}
} // namespace
void build_2d_templates_adaptive(
    const char *outdir =
        "Data_driven/results/templates_2d_adaptive14",
    const char *fitInput =
        "Data_driven/results/pp_data_2d_adaptive14/fit_results.csv",
    const char *mixedInput =
        "Data_driven/results/mixed_dps_5M_seed20260728.root",
    double controlDyMin = 1.8, double controlPhiMax = pi / 2,
    double spsLeakageFraction = 0.0, bool allowSignedWeights = false) {
  gSystem->mkdir(outdir, true);
  auto cells = read(fitInput);
  if (cells.empty()) {
    std::cerr << "no fit cells\n";
    return;
  }
  TFile fm(mixedInput);
  auto *tr = (TTree *)fm.Get("mix");
  double y, p, w;
  tr->SetBranchAddress("delta_y", &y);
  tr->SetBranchAddress("delta_phi", &p);
  tr->SetBranchAddress("evt_weight", &w);
  Long64_t used = 0;
  for (Long64_t i = 0; i < tr->GetEntries(); i++) {
    tr->GetEntry(i);
    int b = find(cells, y, p);
    if (b >= 0 && std::isfinite(w) &&
        (allowSignedWeights ? w != 0 : w > 0)) {
      cells[b].m += w;
      cells[b].vm += w * w;
      used++;
    }
  }
  double A = 0, vA = 0, B = 0, vB = 0;
  for (auto &c : cells)
    if (c.cr(controlDyMin, controlPhiMax)) {
      A += c.n;
      vA += c.en * c.en;
      B += c.m;
      vB += c.vm;
    }
  double ntInput = 0, mtInput = 0;
  for (const auto &c : cells) {
    ntInput += c.n;
    mtInput += c.m;
  }
  const double numerator = A - spsLeakageFraction * ntInput;
  const double denominator = B - spsLeakageFraction * mtInput;
  double a = numerator / denominator;
  double va = vA / (denominator * denominator) +
              numerator * numerator * vB / std::pow(denominator, 4);
  int neg = 0;
  double nt = 0, dt = 0, st = 0, maxc = 0, maxw = 0;
  for (size_t i = 0; i < cells.size(); i++) {
    auto &c = cells[i];
    c.d = a * c.m;
    c.s = c.n - c.d;
    double vd = 0, vs = 0;
    for (size_t j = 0; j < cells.size(); j++) {
      bool same = i == j;
      bool inControl = cells[j].cr(controlDyMin, controlPhiMax);
      double ddn = inControl ? c.m / B : 0,
             ddm = a * (same ? 1 : 0) - (inControl ? a * c.m / B : 0),
             dsn = (same ? 1 : 0) - ddn, dsm = -ddm;
      vd += ddn * ddn * cells[j].en * cells[j].en + ddm * ddm * cells[j].vm;
      vs += dsn * dsn * cells[j].en * cells[j].en + dsm * dsm * cells[j].vm;
    }
    c.ed = std::sqrt(vd);
    c.es = std::sqrt(vs);
    c.wd = c.d / c.n;
    c.ws = c.s / c.n;
    if (c.s < 0)
      neg++;
    nt += c.n;
    dt += c.d;
    st += c.s;
    maxc = std::max(maxc, std::abs(c.s + c.d - c.n));
    maxw = std::max(maxw, std::abs(c.wd + c.ws - 1));
  }
  TH2Poly hd, hm, hp, hs, hwd, hws, hz;
  configure(hd, "h_data_dy_dphi", cells);
  configure(hm, "h_mix_dy_dphi", cells);
  configure(hp, "h_template_dps_dy_dphi", cells);
  configure(hs, "h_template_sps_dy_dphi", cells);
  configure(hwd, "h_weight_dps_dy_dphi", cells);
  configure(hws, "h_weight_sps_dy_dphi", cells);
  configure(hz, "h_sps_significance_dy_dphi", cells);
  for (size_t i = 0; i < cells.size(); i++) {
    int b = i + 1;
    auto &c = cells[i];
    hd.SetBinContent(b, c.n);
    hd.SetBinError(b, c.en);
    hm.SetBinContent(b, c.m);
    hm.SetBinError(b, std::sqrt(c.vm));
    hp.SetBinContent(b, c.d);
    hp.SetBinError(b, c.ed);
    hs.SetBinContent(b, c.s);
    hs.SetBinError(b, c.es);
    hwd.SetBinContent(b, c.wd);
    hws.SetBinContent(b, c.ws);
    hz.SetBinContent(b, significance(c));
  }
  std::string base(outdir);
  TFile fo((base + "/templates_2d.root").c_str(), "RECREATE");
  hd.Write();
  hm.Write();
  hp.Write();
  hs.Write();
  hwd.Write();
  hws.Write();
  hz.Write();
  TTree ct("cells", "adaptive 2D component cells");
  Cell q;
  int control;
  ct.Branch("iy", &q.iy);
  ct.Branch("iphi", &q.ip);
  ct.Branch("dy_low", &q.yl);
  ct.Branch("dy_high", &q.yh);
  ct.Branch("dphi_low", &q.pl);
  ct.Branch("dphi_high", &q.ph);
  ct.Branch("data", &q.n);
  ct.Branch("data_error", &q.en);
  ct.Branch("mix", &q.m);
  ct.Branch("mix_pair_variance", &q.vm);
  ct.Branch("dps", &q.d);
  ct.Branch("dps_error", &q.ed);
  ct.Branch("sps", &q.s);
  ct.Branch("sps_error", &q.es);
  ct.Branch("weight_dps", &q.wd);
  ct.Branch("weight_sps", &q.ws);
  ct.Branch("in_control", &control);
  for (const auto &c : cells) {
    q = c;
    control = c.cr(controlDyMin, controlPhiMax);
    ct.Fill();
  }
  ct.Write();
  fo.Close();
  std::ofstream csv(base + "/cells.csv");
  csv << "cell,iy,iphi,dy_low,dy_high,dphi_low,dphi_high,in_control,data,data_"
         "error,mix,mix_pair_error,dps,dps_error,sps,sps_error,weight_dps,"
         "weight_sps,sps_significance\n"
      << std::setprecision(12);
  for (size_t i = 0; i < cells.size(); i++) {
    auto &c = cells[i];
    csv << i << "," << c.iy << "," << c.ip << "," << c.yl << "," << c.yh << ","
        << c.pl << "," << c.ph << ","
        << c.cr(controlDyMin, controlPhiMax) << "," << c.n << "," << c.en
        << "," << c.m << "," << std::sqrt(c.vm) << "," << c.d << "," << c.ed
        << "," << c.s << "," << c.es << "," << c.wd << "," << c.ws << ","
        << significance(c) << "\n";
  }
  std::ofstream o(base + "/summary.txt");
  o << std::setprecision(12) << "cells=" << cells.size()
    << "\nfit_input=" << fitInput << "\nmixed_input=" << mixedInput
    << "\nmixed_pairs="
    << used
    << "\ncontrol_dy_min=" << controlDyMin
    << "\ncontrol_phi_max=" << controlPhiMax
    << "\nsps_leakage_fraction=" << spsLeakageFraction
    << "\nallow_signed_weights=" << allowSignedWeights
    << "\ndata_control="
    << A << "\ndata_control_error=" << std::sqrt(vA) << "\nmix_control=" << B
    << "\nalpha=" << a << "\nalpha_pair_level_error=" << std::sqrt(va)
    << "\ndata_total=" << nt << "\ndps_total=" << dt << "\nsps_total=" << st
    << "\nf_dps=" << dt / nt << "\nf_sps=" << st / nt
    << "\nnegative_sps_bins=" << neg << "\nmax_abs_bin_closure=" << maxc
    << "\nmax_abs_weight_sum_minus_one=" << maxw << "\n";
  std::vector<double> dataValues, dataErrors, dpsValues, dpsErrors, spsValues,
      spsErrors, dpsWeights, spsWeights, significances, noErrors;
  for (const auto &c : cells) {
    dataValues.push_back(c.n);
    dataErrors.push_back(c.en);
    dpsValues.push_back(c.d);
    dpsErrors.push_back(c.ed);
    spsValues.push_back(c.s);
    spsErrors.push_back(c.es);
    dpsWeights.push_back(c.wd);
    spsWeights.push_back(c.ws);
    significances.push_back(significance(c));
    noErrors.push_back(0);
  }
  draw(hd, cells, base + "/pp_yield_map.pdf",
       "Data prompt-prompt fitted yield: adaptive 14 cells", dataValues,
       dataErrors, true, "%.1f");
  draw(hp, cells, base + "/dps_normalized_map.pdf",
       "Control-region normalized mixed-DPS yield", dpsValues, dpsErrors,
       true, "%.1f");
  draw(hs, cells, base + "/sps_data_driven_map.pdf",
       "Data-driven SPS yield: Data PP - normalized DPS", spsValues,
       spsErrors, true, "%.1f");
  draw(hwd, cells, base + "/dps_weight_map.pdf",
       "DPS component weight (central value only)", dpsWeights, noErrors,
       false, "%.3f");
  draw(hws, cells, base + "/sps_weight_map.pdf",
       "SPS component weight (central value only)", spsWeights, noErrors,
       false, "%.3f");
  draw(hz, cells, base + "/sps_significance_map.pdf",
       "SPS significance: adaptive 14 cells", significances, noErrors, false,
       "%.2f");
  std::ofstream meta(base + "/metadata.txt");
  meta << std::setprecision(12)
       << "artifact=adaptive_2d_SPS_DPS_template_decomposition\n"
       << "status=generated_analysis_artifact\n"
       << "pp_data_selection=defined_by_fit_input\n"
       << "dps_template_construction=defined_by_mixed_input\n"
       << "binning=adaptive_14_cells_abs_delta_y_vs_abs_delta_phi\n"
       << "control_dy_min=" << controlDyMin << "\n"
       << "control_phi_max=" << controlPhiMax << "\n"
       << "sps_leakage_fraction=" << spsLeakageFraction << "\n"
       << "pp_yield_map=center_cells_csv_data;error_cells_csv_data_error\n"
       << "dps_normalized_map=center_cells_csv_dps;error_cells_csv_dps_error\n"
       << "dps_error=includes_control_region_normalization_correlation_in_"
          "existing_cell_propagation\n"
       << "sps_data_driven_map=center_cells_csv_sps;error_cells_csv_sps_error\n"
       << "negative_sps_cell=retained_signed_and_compatible_with_zero\n"
       << "sps_weight_map=center_cells_csv_weight_sps;error_not_displayed\n"
       << "dps_weight_map=center_cells_csv_weight_dps;error_not_displayed\n"
       << "weight_uncertainty_reason=full_ratio_covariance_not_encoded_in_"
          "cells_csv;ROOT_default_sqrt_content_is_not_used\n"
       << "statistical_provenance=not_computed_by_this_generator;see_route_"
          "metadata_and_producer_package\n"
       << "generator=Data_driven/build_2d_templates_adaptive.cpp\n"
       << "fit_input=" << fitInput << "\n"
       << "mixed_input=" << mixedInput << "\n"
       << "allow_signed_weights=" << allowSignedWeights << "\n"
       << "git_commit=" << gSystem->GetFromPipe("git rev-parse HEAD").Data()
       << "\n";
  std::cout << "cells=" << cells.size() << " alpha=" << a << " data=" << nt
            << " dps=" << dt << " sps=" << st << " negative=" << neg
            << "\noutput=" << outdir << std::endl;
}
