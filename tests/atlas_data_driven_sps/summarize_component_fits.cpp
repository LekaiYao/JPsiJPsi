#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
namespace {
struct V {
  std::string n, l;
  std::vector<double> e;
};
struct R {
  std::string var, comp, weight;
  int bin;
  double lo, hi, raw, entries, neg, sw, sw2, negw, y, ey, status, cov, edm, nll,
      attempts, ok;
};
std::vector<std::string> split(const std::string &s) {
  std::vector<std::string> v;
  std::stringstream q(s);
  std::string x;
  while (std::getline(q, x, ","[0]))
    v.push_back(x);
  return v;
}
R read(const std::string &p) {
  std::ifstream f(p);
  std::string h, l;
  std::getline(f, h);
  std::getline(f, l);
  auto x = split(l);
  return {x[0],
          x[1],
          x[2],
          std::stoi(x[3]),
          std::stod(x[4]),
          std::stod(x[5]),
          std::stod(x[6]),
          std::stod(x[7]),
          std::stod(x[8]),
          std::stod(x[9]),
          std::stod(x[10]),
          std::stod(x[11]),
          std::stod(x[12]),
          std::stod(x[13]),
          std::stod(x[14]),
          std::stod(x[15]),
          std::stod(x[16]),
          std::stod(x[17]),
          std::stod(x[18]),
          std::stod(x[19])};
}
R *get(std::vector<R> &r, const std::string &v, int b, const std::string &c) {
  for (auto &x : r)
    if (x.var == v && x.bin == b && x.comp == c)
      return &x;
  return nullptr;
}
void draw(TH1D &t, TH1D &s, TH1D &d, const std::string &p, bool norm) {
  TH1D a(t), b(s), c(d);
  if (norm) {
    if (a.Integral())
      a.Scale(1 / a.Integral());
    if (b.Integral())
      b.Scale(1 / b.Integral());
    if (c.Integral())
      c.Scale(1 / c.Integral());
    a.GetYaxis()->SetTitle("Normalized shape");
  }
  a.SetLineColor(kBlack);
  a.SetMarkerColor(kBlack);
  a.SetMarkerStyle(20);
  b.SetLineColor(kRed + 1);
  b.SetMarkerColor(kRed + 1);
  b.SetMarkerStyle(21);
  c.SetLineColor(kBlue + 1);
  c.SetMarkerColor(kBlue + 1);
  c.SetMarkerStyle(22);
  double mx = std::max({a.GetMaximum(), b.GetMaximum(), c.GetMaximum()});
  a.SetMaximum(1.3 * mx);
  a.SetMinimum(0);
  TCanvas cv("cv", "", 850, 700);
  cv.SetLeftMargin(.13);
  cv.SetBottomMargin(.12);
  a.Draw("E1");
  b.Draw("E1 SAME");
  c.Draw("E1 SAME");
  TLegend leg(.62, .72, .86, .87);
  leg.SetBorderSize(0);
  leg.AddEntry(&a, "Nominal PP", "lep");
  leg.AddEntry(&b, "SPS-weighted PP", "lep");
  leg.AddEntry(&c, "DPS-weighted PP", "lep");
  leg.Draw();
  cv.SaveAs(p.c_str());
}
} // namespace
void summarize_component_fits(
    const char *outdir =
        "tests/atlas_data_driven_sps/results/component_fits_adaptive14",
    const char *onlyVar = "") {
  gSystem->mkdir(outdir, true);
  std::string base(outdir), jobs = base + "/jobs/";
  const double pi = 3.14159265358979323846;
  std::vector<V> vs = {
      {"delta_y", "|#Delta y|", {0, .5, 1, 1.5, 2, 2.5, 4}},
      {"delta_phi",
       "|#Delta#phi|",
       {0, .3927, .7854, 1.1781, 1.5708, 1.9635, 2.3562, 2.7489, 3.1416}},
      {"evt_mass",
       "m(J/#psi J/#psi) [GeV]",
       {7.5, 17.5, 27.5, 37.5, 47.5, 57.5, 67.5, 107.5}},
      {"evt_y", "|y(J/#psi J/#psi)|", {0, .4, .8, 1.2, 1.6, 2}}};
  if (std::string(onlyVar).size()) {
    std::vector<V> selected;
    for (const auto &v : vs)
      if (v.n == onlyVar)
        selected.push_back(v);
    vs.swap(selected);
  }
  std::vector<std::string> cs = {"total", "sps", "dps"};
  std::vector<R> r;
  for (auto &v : vs)
    for (int b = 0; b + 1 < (int)v.e.size(); b++)
      for (auto &c : cs)
        r.push_back(read(jobs + v.n + "_bin" + std::to_string(b) + "_" + c +
                         "/fit_results.csv"));
  std::ofstream all(base + "/fit_results.csv");
  all << "variable,component,weight,bin,low,high,raw_entries,fit_entries,"
         "negative_weight_entries,sum_weights,sum_weights2,sum_negative_"
         "weights,pp_yield,pp_error,status,covQual,edm,minNll,attempts,"
         "accepted\n"
      << std::setprecision(12);
  for (auto &x : r)
    all << x.var << "," << x.comp << "," << x.weight << "," << x.bin << ","
        << x.lo << "," << x.hi << "," << x.raw << "," << x.entries << ","
        << x.neg << "," << x.sw << "," << x.sw2 << "," << x.negw << "," << x.y
        << "," << x.ey << "," << x.status << "," << x.cov << "," << x.edm << ","
        << x.nll << "," << x.attempts << "," << x.ok << "\n";
  std::ofstream cl(base + "/closure.csv"), fr(base + "/fraction_summary.csv");
  cl << "variable,bin,low,high,total,total_error,sps,sps_error,dps,dps_error,"
        "sum_components,difference,relative_difference\n"
     << std::setprecision(12);
  fr << "variable,sum_total,sum_sps,sum_dps,sum_components,f_sps_components,f_"
        "dps_components,components_minus_total,relative_closure,max_abs_bin_"
        "relative_closure\n"
     << std::setprecision(12);
  TFile fo((base + "/component_yields.root").c_str(), "RECREATE");
  double globalMax = 0;
  int accepted = 0;
  for (auto &x : r)
    accepted += x.ok;
  for (auto &v : vs) {
    int nb = v.e.size() - 1;
    TH1D ht(("h_" + v.n + "_total").c_str(),
            (";" + v.l + ";PP corrected yield").c_str(), nb, v.e.data()),
        hs(("h_" + v.n + "_sps").c_str(),
           (";" + v.l + ";PP corrected yield").c_str(), nb, v.e.data()),
        hd(("h_" + v.n + "_dps").c_str(),
           (";" + v.l + ";PP corrected yield").c_str(), nb, v.e.data());
    double st = 0, ss = 0, sd = 0, mx = 0;
    for (int b = 0; b < nb; b++) {
      auto *t = get(r, v.n, b, "total"), *s = get(r, v.n, b, "sps"),
           *d = get(r, v.n, b, "dps");
      ht.SetBinContent(b + 1, t->y);
      ht.SetBinError(b + 1, t->ey);
      hs.SetBinContent(b + 1, s->y);
      hs.SetBinError(b + 1, s->ey);
      hd.SetBinContent(b + 1, d->y);
      hd.SetBinError(b + 1, d->ey);
      double sum = s->y + d->y, diff = sum - t->y,
             rel = t->y ? diff / t->y : NAN;
      mx = std::max(mx, std::abs(rel));
      cl << v.n << "," << b << "," << v.e[b] << "," << v.e[b + 1] << "," << t->y
         << "," << t->ey << "," << s->y << "," << s->ey << "," << d->y << ","
         << d->ey << "," << sum << "," << diff << "," << rel << "\n";
      st += t->y;
      ss += s->y;
      sd += d->y;
    }
    globalMax = std::max(globalMax, mx);
    fr << v.n << "," << st << "," << ss << "," << sd << "," << ss + sd << ","
       << ss / (ss + sd) << "," << sd / (ss + sd) << "," << ss + sd - st << ","
       << (ss + sd - st) / st << "," << mx << "\n";
    fo.cd();
    ht.Write();
    hs.Write();
    hd.Write();
    draw(ht, hs, hd, base + "/" + v.n + "_component_yields.pdf", false);
    draw(ht, hs, hd, base + "/" + v.n + "_component_shapes_normalized.pdf",
         true);
  }
  fo.Close();
  std::ofstream s(base + "/summary.txt");
  s << "fits_expected=" << r.size() << "\nfits_loaded=" << r.size()
    << "\nfits_accepted=" << accepted << "\nfits_failed=" << r.size() - accepted
    << "\nvariables=" << (std::string(onlyVar).size() ? onlyVar
                                                     : "delta_y,delta_phi,evt_mass,evt_y")
    << "\nbinning=mainline_"
       "Makefile\nmax_abs_bin_relative_sps_plus_dps_minus_total="
    << std::setprecision(12) << globalMax
    << "\nwarning=SPS_and_DPS_fit_errors_are_correlated_and_component_fraction_"
       "uncertainty_is_not_included\nwarning=signed_SPS_weights_are_accepted_"
       "by_RooFit_but_require_toy_validation_for_coverage\n";
  std::cout << "loaded=" << r.size() << " accepted=" << accepted
            << " max closure=" << globalMax << "\noutput=" << outdir
            << std::endl;
}
