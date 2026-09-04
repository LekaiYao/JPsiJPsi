#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TNamed.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVectorD.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

constexpr double kPi = 3.14159265358979323846;

struct EventKey {
  ULong64_t run = 0, lumi = 0, event = 0;
  bool operator==(const EventKey &other) const {
    return run == other.run && lumi == other.lumi && event == other.event;
  }
};

struct EventKeyHash {
  std::size_t operator()(const EventKey &key) const {
    auto mix = [](std::size_t seed, ULong64_t value) {
      return seed ^ (std::hash<ULong64_t>{}(value) + 0x9e3779b97f4a7c15ULL +
                     (seed << 6) + (seed >> 2));
    };
    std::size_t seed = 0;
    seed = mix(seed, key.run);
    seed = mix(seed, key.lumi);
    return mix(seed, key.event);
  }
};

struct Values {
  double pt1 = 0, y1 = 0, pt2 = 0, y2 = 0;
  double mass = 0, pairY = 0, dy = 0, dphi = 0;
};

struct Variable {
  std::string key;
  std::string label;
  std::vector<double> edges;
  std::function<double(const Values &)> value;
  bool logY = false;
};

struct RawSample {
  std::vector<std::vector<double>> sumw, sumw2;
  std::vector<Long64_t> accepted, outside;
  std::vector<double> outsideWeight;
  Long64_t inputEntries = 0, invalidWeights = 0;
  std::vector<std::vector<double>> incident;
  std::unordered_map<EventKey, std::size_t, EventKeyHash> clusterIndex;
};

struct Shape {
  std::vector<double> fraction;
  TMatrixDSym covariance;
  double totalWeight = 0, totalWeight2 = 0, effectiveEntries = 0;
  Shape(int bins) : fraction(bins, 0), covariance(bins) { covariance.Zero(); }
};

struct Metrics {
  double maxCdf = 0, maxCdfEdge = 0;
  double jsDistance = std::numeric_limits<double>::quiet_NaN();
  double chi2 = 0, pValue = 0, maxRatioDeviation = 0;
  int covarianceRank = 0, droppedBin = -1;
};

std::vector<Variable> variables() {
  return {
      {"jpsi1_pt", "p_{T}(J/#psi_{1}) [GeV]",
       {10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 35, 40},
       [](const Values &v) { return v.pt1; }},
      {"jpsi1_abs_y", "|y(J/#psi_{1})|",
       {0, .2, .4, .6, .8, 1., 1.2, 1.4, 1.6, 1.8, 2.},
       [](const Values &v) { return std::abs(v.y1); }},
      {"jpsi2_pt", "p_{T}(J/#psi_{2}) [GeV]",
       {10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 35, 40},
       [](const Values &v) { return v.pt2; }},
      {"jpsi2_abs_y", "|y(J/#psi_{2})|",
       {0, .2, .4, .6, .8, 1., 1.2, 1.4, 1.6, 1.8, 2.},
       [](const Values &v) { return std::abs(v.y2); }},
      {"delta_y", "|#Delta y|",
       {0, .25, .5, .75, 1., 1.5, 2., 2.5, 3., 4.},
       [](const Values &v) { return v.dy; }},
      {"delta_phi", "|#Delta#phi|",
       {0, kPi / 8, 2 * kPi / 8, 3 * kPi / 8, 4 * kPi / 8,
        5 * kPi / 8, 6 * kPi / 8, 7 * kPi / 8, kPi},
       [](const Values &v) { return v.dphi; }},
      {"evt_mass", "m(J/#psi J/#psi) [GeV]",
       {7.5, 10, 12.5, 15, 20, 25, 30, 40, 50, 70, 100, 150, 220, 320},
       [](const Values &v) { return v.mass; }, true},
      {"evt_abs_y", "|y(J/#psi J/#psi)|",
       {0, .2, .4, .6, .8, 1., 1.2, 1.4, 1.6, 1.8, 2.},
       [](const Values &v) { return std::abs(v.pairY); }},
      {"evt_pt", "p_{T}(J/#psi J/#psi) [GeV]",
       {0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 80},
       [](const Values &v) {
         const double pt2 = v.pt1 * v.pt1 + v.pt2 * v.pt2 +
                            2 * v.pt1 * v.pt2 * std::cos(v.dphi);
         return std::sqrt(std::max(0.0, pt2));
       }}};
}

int findBin(const std::vector<double> &edges, double value) {
  if (!std::isfinite(value) || value < edges.front() || value > edges.back())
    return -1;
  if (value == edges.back()) return static_cast<int>(edges.size()) - 2;
  return static_cast<int>(std::upper_bound(edges.begin(), edges.end(), value) -
                          edges.begin()) - 1;
}

void requireBranches(TTree *tree, bool sourceIds) {
  const char *common[] = {"Jpsi_pt1", "Jpsi_y1", "Jpsi_pt2", "Jpsi_y2",
                          "evt_mass", "evt_y", "delta_y", "delta_phi",
                          "evt_weight"};
  for (const char *name : common)
    if (!tree->GetBranch(name)) throw std::runtime_error("Missing branch " + std::string(name));
  if (sourceIds) {
    const char *ids[] = {"run1", "lumi1", "event1", "run2", "lumi2", "event2"};
    for (const char *name : ids)
      if (!tree->GetBranch(name)) throw std::runtime_error("Missing branch " + std::string(name));
  }
}

RawSample readSample(const char *path, const std::vector<Variable> &vars,
                     bool sourceIds) {
  TFile file(path, "READ");
  if (file.IsZombie()) throw std::runtime_error("Cannot open " + std::string(path));
  auto *tree = dynamic_cast<TTree *>(file.Get("mix"));
  if (!tree) throw std::runtime_error("Missing mix tree in " + std::string(path));
  requireBranches(tree, sourceIds);

  RawSample out;
  out.sumw.resize(vars.size()); out.sumw2.resize(vars.size());
  out.accepted.assign(vars.size(), 0); out.outside.assign(vars.size(), 0);
  out.outsideWeight.assign(vars.size(), 0);
  std::vector<int> offsets(vars.size(), 0);
  int totalBins = 0;
  for (std::size_t iv = 0; iv < vars.size(); ++iv) {
    offsets[iv] = totalBins;
    const int n = static_cast<int>(vars[iv].edges.size()) - 1;
    totalBins += n;
    out.sumw[iv].assign(n, 0); out.sumw2[iv].assign(n, 0);
  }

  Values v; double weight = 0;
  EventKey a, b;
  tree->SetBranchAddress("Jpsi_pt1", &v.pt1); tree->SetBranchAddress("Jpsi_y1", &v.y1);
  tree->SetBranchAddress("Jpsi_pt2", &v.pt2); tree->SetBranchAddress("Jpsi_y2", &v.y2);
  tree->SetBranchAddress("evt_mass", &v.mass); tree->SetBranchAddress("evt_y", &v.pairY);
  tree->SetBranchAddress("delta_y", &v.dy); tree->SetBranchAddress("delta_phi", &v.dphi);
  tree->SetBranchAddress("evt_weight", &weight);
  if (sourceIds) {
    tree->SetBranchAddress("run1", &a.run); tree->SetBranchAddress("lumi1", &a.lumi);
    tree->SetBranchAddress("event1", &a.event); tree->SetBranchAddress("run2", &b.run);
    tree->SetBranchAddress("lumi2", &b.lumi); tree->SetBranchAddress("event2", &b.event);
  }
  auto cluster = [&](const EventKey &key) {
    auto found = out.clusterIndex.find(key);
    if (found != out.clusterIndex.end()) return found->second;
    const std::size_t id = out.incident.size();
    out.clusterIndex.emplace(key, id);
    out.incident.emplace_back(totalBins, 0);
    return id;
  };

  out.inputEntries = tree->GetEntries();
  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
    tree->GetEntry(entry);
    if (!std::isfinite(weight)) { ++out.invalidWeights; continue; }
    std::size_t ia = 0, ib = 0;
    if (sourceIds) { ia = cluster(a); ib = cluster(b); }
    for (std::size_t iv = 0; iv < vars.size(); ++iv) {
      const int bin = findBin(vars[iv].edges, vars[iv].value(v));
      if (bin < 0) { ++out.outside[iv]; out.outsideWeight[iv] += weight; continue; }
      ++out.accepted[iv]; out.sumw[iv][bin] += weight;
      out.sumw2[iv][bin] += weight * weight;
      if (sourceIds) {
        out.incident[ia][offsets[iv] + bin] += weight;
        if (ib != ia) out.incident[ib][offsets[iv] + bin] += weight;
      }
    }
  }
  file.Close();
  return out;
}

Shape normalizedIndependent(const RawSample &raw, std::size_t iv) {
  const int n = raw.sumw[iv].size(); Shape out(n);
  out.totalWeight = 0; out.totalWeight2 = 0;
  for (int i = 0; i < n; ++i) { out.totalWeight += raw.sumw[iv][i]; out.totalWeight2 += raw.sumw2[iv][i]; }
  if (!(out.totalWeight > 0)) throw std::runtime_error("Non-positive MC normalization");
  out.effectiveEntries = out.totalWeight * out.totalWeight / out.totalWeight2;
  for (int i = 0; i < n; ++i) out.fraction[i] = raw.sumw[iv][i] / out.totalWeight;
  const double s = out.totalWeight, s2 = s * s, s3 = s2 * s, s4 = s2 * s2;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      const double diagonal = i == j ? raw.sumw2[iv][i] / s2 : 0;
      out.covariance(i, j) = diagonal -
          (raw.sumw2[iv][i] * raw.sumw[iv][j] + raw.sumw2[iv][j] * raw.sumw[iv][i]) / s3 +
          raw.sumw[iv][i] * raw.sumw[iv][j] * out.totalWeight2 / s4;
    }
  return out;
}

Shape normalizedClusterJackknife(const RawSample &raw,
                                 const std::vector<Variable> &vars,
                                 std::size_t iv) {
  const int n = raw.sumw[iv].size(); Shape out(n);
  for (int i = 0; i < n; ++i) { out.totalWeight += raw.sumw[iv][i]; out.totalWeight2 += raw.sumw2[iv][i]; }
  if (!(out.totalWeight > 0) || raw.incident.size() < 3)
    throw std::runtime_error("Invalid mixing normalization or cluster count");
  out.effectiveEntries = out.totalWeight * out.totalWeight / out.totalWeight2;
  for (int i = 0; i < n; ++i) out.fraction[i] = raw.sumw[iv][i] / out.totalWeight;
  int offset = 0;
  for (std::size_t k = 0; k < iv; ++k) offset += vars[k].edges.size() - 1;
  std::vector<std::vector<double>> leave;
  leave.reserve(raw.incident.size());
  std::vector<double> mean(n, 0);
  for (const auto &cluster : raw.incident) {
    double removed = 0;
    for (int i = 0; i < n; ++i) removed += cluster[offset + i];
    const double denominator = out.totalWeight - removed;
    if (!(denominator > 0) || !std::isfinite(denominator))
      throw std::runtime_error("Invalid leave-one-source-event normalization");
    leave.emplace_back(n, 0);
    for (int i = 0; i < n; ++i) {
      leave.back()[i] = (raw.sumw[iv][i] - cluster[offset + i]) / denominator;
      mean[i] += leave.back()[i];
    }
  }
  for (double &x : mean) x /= leave.size();
  const double scale = (leave.size() - 1.0) / leave.size();
  for (const auto &sample : leave)
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        out.covariance(i, j) += scale * (sample[i] - mean[i]) * (sample[j] - mean[j]);
  return out;
}

Metrics metrics(const Shape &mix, const Shape &mc, const Variable &var) {
  const int n = mix.fraction.size(); Metrics out;
  double cumulative = 0;
  for (int i = 0; i < n - 1; ++i) {
    cumulative += mix.fraction[i] - mc.fraction[i];
    if (std::abs(cumulative) > out.maxCdf) {
      out.maxCdf = std::abs(cumulative); out.maxCdfEdge = var.edges[i + 1];
    }
  }
  bool positive = true; double js = 0;
  for (int i = 0; i < n; ++i) {
    const double p = mix.fraction[i], q = mc.fraction[i];
    if (p < -1e-12 || q < -1e-12) { positive = false; break; }
    const double pp = std::max(0.0, p), qq = std::max(0.0, q), m = .5 * (pp + qq);
    if (pp > 0) js += .5 * pp * std::log(pp / m);
    if (qq > 0) js += .5 * qq * std::log(qq / m);
  }
  if (positive) out.jsDistance = std::sqrt(std::max(0.0, js));

  out.droppedBin = 0;
  for (int i = 1; i < n; ++i)
    if (mix.fraction[i] + mc.fraction[i] > mix.fraction[out.droppedBin] + mc.fraction[out.droppedBin])
      out.droppedBin = i;
  std::vector<int> kept;
  for (int i = 0; i < n; ++i) if (i != out.droppedBin) kept.push_back(i);
  TMatrixDSym covariance(n - 1); covariance.Zero(); TVectorD difference(n - 1);
  for (int i = 0; i < n - 1; ++i) {
    difference[i] = mix.fraction[kept[i]] - mc.fraction[kept[i]];
    for (int j = 0; j < n - 1; ++j)
      covariance(i, j) = mix.covariance(kept[i], kept[j]) + mc.covariance(kept[i], kept[j]);
  }
  TMatrixDSymEigen eigen(covariance);
  const TVectorD values = eigen.GetEigenValues(); const TMatrixD vectors = eigen.GetEigenVectors();
  double maximum = 0;
  for (int k = 0; k < n - 1; ++k) maximum = std::max(maximum, values[k]);
  const double threshold = std::max(1e-18, maximum * 1e-10);
  for (int k = 0; k < n - 1; ++k) {
    if (values[k] <= threshold) continue;
    double projection = 0;
    for (int i = 0; i < n - 1; ++i) projection += vectors(i, k) * difference[i];
    out.chi2 += projection * projection / values[k]; ++out.covarianceRank;
  }
  out.pValue = out.covarianceRank > 0 ? TMath::Prob(out.chi2, out.covarianceRank) : 0;
  for (int i = 0; i < n; ++i)
    if (mc.fraction[i] > 0)
      out.maxRatioDeviation = std::max(out.maxRatioDeviation,
                                       std::abs(mix.fraction[i] / mc.fraction[i] - 1));
  return out;
}

std::string format(double value, int precision = 3) {
  std::ostringstream out; out << std::setprecision(precision) << std::fixed << value; return out.str();
}

std::string formatProbability(double value) {
  std::ostringstream out;
  if (value < 0.001) out << std::scientific << std::setprecision(2) << value;
  else out << std::fixed << std::setprecision(3) << value;
  return out.str();
}

void draw(const std::string &base, const Variable &var, const Shape &mix,
          const Shape &mc, const Metrics &metric, TFile &rootOutput,
          const char *mixLegend, const char *referenceLegend,
          const char *analysisTitle, bool diagnosticCovariance) {
  const int n = var.edges.size() - 1;
  TH1D hMix(("h_mixing_" + var.key).c_str(), "", n, var.edges.data());
  TH1D hMc(("h_dps_mc_" + var.key).c_str(), "", n, var.edges.data());
  TH1D ratio(("h_ratio_" + var.key).c_str(), "", n, var.edges.data());
  TH1D band(("h_mc_ratio_band_" + var.key).c_str(), "", n, var.edges.data());
  hMix.SetDirectory(nullptr); hMc.SetDirectory(nullptr); ratio.SetDirectory(nullptr); band.SetDirectory(nullptr);
  double ymax = 0, yminPositive = std::numeric_limits<double>::max();
  double ratioLow = 1, ratioHigh = 1;
  for (int i = 0; i < n; ++i) {
    const double width = var.edges[i + 1] - var.edges[i];
    const double emix = std::sqrt(std::max(0.0, mix.covariance(i, i)));
    const double emc = std::sqrt(std::max(0.0, mc.covariance(i, i)));
    hMix.SetBinContent(i + 1, mix.fraction[i] / width); hMix.SetBinError(i + 1, emix / width);
    hMc.SetBinContent(i + 1, mc.fraction[i] / width); hMc.SetBinError(i + 1, emc / width);
    ymax = std::max({ymax, (mix.fraction[i] + emix) / width, (mc.fraction[i] + emc) / width});
    if (mix.fraction[i] > 0) yminPositive = std::min(yminPositive, mix.fraction[i] / width);
    if (mc.fraction[i] > 0) yminPositive = std::min(yminPositive, mc.fraction[i] / width);
    if (mc.fraction[i] > 0) {
      const double r = mix.fraction[i] / mc.fraction[i];
      const double erMix = emix / mc.fraction[i], erMc = emc / mc.fraction[i];
      ratio.SetBinContent(i + 1, r); ratio.SetBinError(i + 1, erMix);
      band.SetBinContent(i + 1, 1); band.SetBinError(i + 1, erMc);
      const double total = std::sqrt(erMix * erMix + erMc * erMc);
      ratioLow = std::min(ratioLow, r - total); ratioHigh = std::max(ratioHigh, r + total);
    }
  }

  TCanvas canvas(("c_" + var.key).c_str(), "", 850, 820);
  auto *top = new TPad(("top_" + var.key).c_str(), "", 0, .30, 1, 1);
  auto *bottom = new TPad(("bottom_" + var.key).c_str(), "", 0, 0, 1, .30);
  top->SetLeftMargin(.14); top->SetRightMargin(.04);
  top->SetTopMargin(.08); top->SetBottomMargin(.02);
  bottom->SetLeftMargin(.14); bottom->SetRightMargin(.04);
  bottom->SetTopMargin(.03); bottom->SetBottomMargin(.35);
  top->Draw(); bottom->Draw(); top->cd(); if (var.logY) top->SetLogy();

  hMc.SetLineColor(kAzure + 2); hMc.SetFillColorAlpha(kAzure - 9, .55); hMc.SetLineWidth(2);
  hMix.SetLineColor(kBlack); hMix.SetMarkerColor(kBlack); hMix.SetMarkerStyle(20); hMix.SetMarkerSize(.9);
  hMc.GetYaxis()->SetTitle("Normalized density"); hMc.GetYaxis()->SetTitleSize(.055);
  hMc.GetYaxis()->SetLabelSize(.045); hMc.GetYaxis()->SetTitleOffset(1.15);
  hMc.GetXaxis()->SetLabelSize(0);
  if (var.logY) {
    hMc.SetMinimum(std::max(1e-10, yminPositive * .2)); hMc.SetMaximum(ymax * 25);
  } else { hMc.SetMinimum(0); hMc.SetMaximum(ymax * 1.55); }
  hMc.Draw("E2"); hMc.Draw("HIST SAME"); hMix.Draw("E1 SAME");

  TLegend legend(.57, .70, .94, .88); legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextSize(.040);
  legend.AddEntry(&hMix, mixLegend, "lep");
  legend.AddEntry(&hMc, referenceLegend, "lf"); legend.Draw();
  TLatex latex; latex.SetNDC(); latex.SetTextFont(42); latex.SetTextSize(.047);
  latex.DrawLatex(.14, .94, analysisTitle); latex.SetTextSize(.035);
  latex.DrawLatex(.14, .875, ("max CDF distance = " + format(metric.maxCdf, 4)).c_str());
  latex.DrawLatex(.14, .825,
                  ((diagnosticCovariance ? "diagnostic #chi^{2}/rank = "
                                         : "covariance #chi^{2}/rank = ") +
                   format(metric.chi2, 2) + "/" +
                   std::to_string(metric.covarianceRank) +
                   ", p = " + formatProbability(metric.pValue)).c_str());
  if (std::isfinite(metric.jsDistance))
    latex.DrawLatex(.14, .775, ("Jensen-Shannon distance = " + format(metric.jsDistance, 4)).c_str());

  bottom->cd();
  band.SetFillColorAlpha(kAzure - 9, .65); band.SetLineColor(kAzure + 2);
  band.GetYaxis()->SetTitle("Mix / MC"); band.GetXaxis()->SetTitle(var.label.c_str());
  band.GetYaxis()->SetNdivisions(505); band.GetYaxis()->SetTitleSize(.12); band.GetYaxis()->SetLabelSize(.10);
  band.GetYaxis()->SetTitleOffset(.48); band.GetXaxis()->SetTitleSize(.13); band.GetXaxis()->SetLabelSize(.105);
  band.GetXaxis()->SetTitleOffset(1.05);
  const double low = std::max(0.0, std::min(0.65, ratioLow - .10));
  const double high = std::max(1.35, ratioHigh + .10);
  band.SetMinimum(low); band.SetMaximum(high); band.Draw("E2");
  ratio.SetMarkerStyle(20); ratio.SetMarkerSize(.8); ratio.SetLineColor(kBlack); ratio.Draw("E1 SAME");
  TLine line(var.edges.front(), 1, var.edges.back(), 1); line.SetLineStyle(2); line.Draw();

  canvas.SaveAs((base + "/plots/" + var.key + ".pdf").c_str());
  rootOutput.cd(); auto *directory = rootOutput.mkdir(var.key.c_str()); directory->cd();
  hMix.Write(); hMc.Write(); ratio.Write(); band.Write();
  TMatrixDSym mixCov(mix.covariance), mcCov(mc.covariance);
  mixCov.Write("cov_mixing_source_event_jackknife"); mcCov.Write("cov_dps_mc_weighted_normalized");
}

}  // namespace

void compare_mixed_dps_to_mc_1d(
    const char *output,
    const char *mixing = "Data_driven/results/nominal_root640_corrected_combcombfix0_20260902_v1/mixed_dps_allpairs.root",
    const char *dpsMc = "Data_driven/results/dps_reference/nominal_dps_template_input.root",
    const char *mixSource = "original_data",
    const char *mixLegend = "single-J/#psi event mixing",
    const char *referenceLegend = "DPS MC",
    const char *analysisTitle = "DPS shape validation",
    const char *sampleRelation = "independent_samples") {
  try {
    const bool diagnosticCovariance =
        std::string(sampleRelation) != "independent_samples";
    gStyle->SetOptStat(0); gStyle->SetEndErrorSize(4);
    const std::vector<Variable> vars = variables();
    const RawSample mixedRaw = readSample(mixing, vars, true);
    const RawSample mcRaw = readSample(dpsMc, vars, false);
    gSystem->mkdir(output, true); gSystem->mkdir((std::string(output) + "/plots").c_str(), true);
    TFile rootOutput((std::string(output) + "/comparison.root").c_str(), "RECREATE");
    std::ofstream metricCsv(std::string(output) + "/shape_metrics.csv");
    std::ofstream binCsv(std::string(output) + "/bin_values.csv");
    std::ofstream summary(std::string(output) + "/summary.txt");
    metricCsv << "variable,bins,mix_input_entries,mix_entries_used,mix_outside_entries,mix_sum_weights,mix_sum_weights2,mix_source_events,mix_pair_naive_neff,mc_input_entries,mc_entries_used,mc_outside_entries,mc_sum_weights,mc_sum_weights2,mc_effective_entries,max_cdf_distance,max_cdf_edge,js_distance,covariance_chi2,covariance_rank,asymptotic_p_value,dropped_covariance_bin,max_abs_ratio_minus_one\n";
    binCsv << "variable,bin,low,high,mix_fraction,mix_error,mc_fraction,mc_error,ratio,mix_ratio_error,mc_ratio_band_error\n";
    summary << std::setprecision(12)
            << "status=complete_pending_user_confirmation\n"
            << "artifact=normalized_1d_event_mixing_vs_dps_mc_shape_validation\n"
            << "mixing_input=" << mixing << "\n"
            << "dps_mc_input=" << dpsMc << "\n"
            << "variables=jpsi1_pt,jpsi1_abs_y,jpsi2_pt,jpsi2_abs_y,delta_y,delta_phi,evt_mass,evt_abs_y,evt_pt\n"
            << "plots=9\n"
            << "mix_source_type=" << mixSource << "\n"
            << "mix_legend=" << mixLegend << "\n"
            << "reference_legend=" << referenceLegend << "\n"
            << "analysis_title=" << analysisTitle << "\n"
            << "sample_relation=" << sampleRelation << "\n"
            << "covariance_chi2_interpretation="
            << (diagnosticCovariance
                    ? "diagnostic_only_shared_source_cross_covariance_not_propagated"
                    : "standard_independent_sample_combination")
            << "\n"
            << "mix_uncertainty=delete_one_" << mixSource
            << "_source_event_cluster_jackknife\n"
            << "dps_mc_uncertainty=independent_weighted_event_normalized_delta_covariance\n"
            << "ratio_points=mixing_uncertainty_only\nratio_band=dps_mc_uncertainty_only\n"
            << "chi2=combined_full_normalized_covariance_eigen_pseudoinverse\n"
            << "max_cdf=binned_maximum_absolute_cumulative_fraction_difference\n"
            << "js_distance=natural_log_jensen_shannon_distance\n"
            << "mix_source_events=" << mixedRaw.incident.size() << "\n"
            << "mix_invalid_weights=" << mixedRaw.invalidWeights << "\n"
            << "mc_invalid_weights=" << mcRaw.invalidWeights << "\n";

    for (std::size_t iv = 0; iv < vars.size(); ++iv) {
      const Shape mix = normalizedClusterJackknife(mixedRaw, vars, iv);
      const Shape mc = normalizedIndependent(mcRaw, iv);
      const Metrics result = metrics(mix, mc, vars[iv]);
      draw(output, vars[iv], mix, mc, result, rootOutput, mixLegend,
           referenceLegend, analysisTitle, diagnosticCovariance);
      metricCsv << std::setprecision(12) << vars[iv].key << ',' << mix.fraction.size() << ','
                << mixedRaw.inputEntries << ',' << mixedRaw.accepted[iv] << ',' << mixedRaw.outside[iv] << ','
                << mix.totalWeight << ',' << mix.totalWeight2 << ',' << mixedRaw.incident.size() << ','
                << mix.effectiveEntries << ',' << mcRaw.inputEntries << ',' << mcRaw.accepted[iv] << ','
                << mcRaw.outside[iv] << ',' << mc.totalWeight << ',' << mc.totalWeight2 << ','
                << mc.effectiveEntries << ',' << result.maxCdf << ',' << result.maxCdfEdge << ','
                << result.jsDistance << ',' << result.chi2 << ',' << result.covarianceRank << ','
                << result.pValue << ',' << result.droppedBin << ',' << result.maxRatioDeviation << '\n';
      summary << vars[iv].key << "_max_cdf_distance=" << result.maxCdf << '\n'
              << vars[iv].key << "_max_cdf_edge=" << result.maxCdfEdge << '\n'
              << vars[iv].key << "_js_distance=" << result.jsDistance << '\n'
              << vars[iv].key << "_covariance_chi2=" << result.chi2 << '\n'
              << vars[iv].key << "_covariance_rank=" << result.covarianceRank << '\n'
              << vars[iv].key << "_asymptotic_p_value=" << result.pValue << '\n';
      for (std::size_t i = 0; i < mix.fraction.size(); ++i) {
        const double emix = std::sqrt(std::max(0.0, mix.covariance(i, i)));
        const double emc = std::sqrt(std::max(0.0, mc.covariance(i, i)));
        const double ratio = mc.fraction[i] > 0 ? mix.fraction[i] / mc.fraction[i] : std::numeric_limits<double>::quiet_NaN();
        binCsv << std::setprecision(12) << vars[iv].key << ',' << i << ',' << vars[iv].edges[i] << ','
               << vars[iv].edges[i + 1] << ',' << mix.fraction[i] << ',' << emix << ','
               << mc.fraction[i] << ',' << emc << ',' << ratio << ','
               << (mc.fraction[i] > 0 ? emix / mc.fraction[i] : 0) << ','
               << (mc.fraction[i] > 0 ? emc / mc.fraction[i] : 0) << '\n';
      }
    }
    summary << "metric_scope=shape_only_unit_normalized\n"
            << "sweight_fit_parameter_uncertainty="
            << (std::string(mixSource) == "original_data"
                    ? "not_propagated_conditional_on_saved_sweights"
                    : "not_applicable")
            << "\n"
            << "systematic_interpretation=false\n";
    TNamed metadata(
        "analysis_definition",
        "unit-normalized one-dimensional event-mixing versus DPS-MC shape "
        "comparison; mixing errors use delete-one-source-event cluster "
        "jackknife; reference errors use normalized weighted-event covariance");
    metadata.Write();
    rootOutput.Close(); metricCsv.close(); binCsv.close(); summary.close();
    std::cout << "DPS_SHAPE_COMPARISON_COMPLETE output=" << output << std::endl;
  } catch (const std::exception &error) {
    std::cerr << "DPS_SHAPE_COMPARISON_FAILED " << error.what() << std::endl;
  }
}
