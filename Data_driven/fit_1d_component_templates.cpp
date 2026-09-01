#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
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
struct Y {
  double t = 0, et = 0, s = 0, es = 0, d = 0, ed = 0;
};
struct F {
  double a = 0, b = 0, ea = 0, eb = 0, corr = 0, f = 0, ef = 0, chi2 = 0;
  int ndf = 0;
  bool ok = false;
};
std::vector<std::string> split(const std::string &s) {
  std::vector<std::string> v;
  std::stringstream q(s);
  std::string x;
  while (std::getline(q, x, ","[0]))
    v.push_back(x);
  return v;
}
std::map<std::string, std::vector<Y>>
load(const std::vector<V> &vs, const std::string &fitInput) {
  std::map<std::string, std::vector<Y>> m;
  for (auto &v : vs)
    m[v.n].resize(v.e.size() - 1);
  std::ifstream f(fitInput.c_str());
  std::string l;
  std::getline(f, l);
  while (std::getline(f, l)) {
    auto x = split(l);
    if (x.size() < 20)
      continue;
    auto &y = m[x[0]][std::stoi(x[3])];
    double z = std::stod(x[12]), e = std::stod(x[13]);
    if (x[1] == "total") {
      y.t = z;
      y.et = e;
    } else if (x[1] == "sps") {
      y.s = z;
      y.es = e;
    } else {
      y.d = z;
      y.ed = e;
    }
  }
  return m;
}
F fit(const std::vector<Y> &y, int skip = -1) {
  double ss = 0, dd = 0;
  for (auto &x : y) {
    ss += x.s;
    dd += x.d;
  }
  double m11 = 0, m12 = 0, m22 = 0, v1 = 0, v2 = 0;
  int n = 0;
  for (size_t i = 0; i < y.size(); i++) {
    if ((int)i == skip)
      continue;
    double s = y[i].s / ss, d = y[i].d / dd, w = 1 / (y[i].et * y[i].et);
    m11 += w * s * s;
    m12 += w * s * d;
    m22 += w * d * d;
    v1 += w * s * y[i].t;
    v2 += w * d * y[i].t;
    n++;
  }
  double det = m11 * m22 - m12 * m12;
  F r;
  if (det <= 0)
    return r;
  r.a = (v1 * m22 - v2 * m12) / det;
  r.b = (v2 * m11 - v1 * m12) / det;
  double va = m22 / det, vb = m11 / det, cab = -m12 / det;
  r.ea = std::sqrt(va);
  r.eb = std::sqrt(vb);
  r.corr = cab / std::sqrt(va * vb);
  double den = r.a + r.b;
  r.f = r.a / den;
  double ga = r.b / (den * den), gb = -r.a / (den * den);
  r.ef = std::sqrt(ga * ga * va + gb * gb * vb + 2 * ga * gb * cab);
  for (size_t i = 0; i < y.size(); i++) {
    if ((int)i == skip)
      continue;
    double pred = r.a * y[i].s / ss + r.b * y[i].d / dd;
    r.chi2 += std::pow((y[i].t - pred) / y[i].et, 2);
  }
  r.ndf = n - 2;
  r.ok = r.a >= 0 && r.b >= 0;
  return r;
}
void draw(const V &v, const std::vector<Y> &y, const F &r,
          const std::string &p) {
  double ss = 0, dd = 0;
  for (auto &x : y) {
    ss += x.s;
    dd += x.d;
  }
  TH1D h("hdata", "", v.e.size() - 1, v.e.data()),
      hs("hs", "", v.e.size() - 1, v.e.data()),
      hd("hd", "", v.e.size() - 1, v.e.data()),
      hp("hp", "", v.e.size() - 1, v.e.data());
  for (size_t i = 0; i < y.size(); i++) {
    h.SetBinContent(i + 1, y[i].t);
    h.SetBinError(i + 1, y[i].et);
    hs.SetBinContent(i + 1, r.a * y[i].s / ss);
    hd.SetBinContent(i + 1, r.b * y[i].d / dd);
    hp.SetBinContent(i + 1, hs.GetBinContent(i + 1) + hd.GetBinContent(i + 1));
  }
  h.SetMarkerStyle(20);
  h.SetLineColor(kBlack);
  hs.SetLineColor(kRed + 1);
  hs.SetLineWidth(2);
  hd.SetLineColor(kBlue + 1);
  hd.SetLineWidth(2);
  hp.SetLineColor(kGreen + 2);
  hp.SetLineWidth(3);
  h.SetTitle("");
  h.GetXaxis()->SetTitle(v.l.c_str());
  h.GetYaxis()->SetTitle("PP corrected yield");
  h.SetMaximum(1.3 * std::max(h.GetMaximum(), hp.GetMaximum()));
  TCanvas c("ct", "", 850, 700);
  c.SetLeftMargin(.13);
  c.SetBottomMargin(.12);
  h.Draw("E1");
  hp.Draw("HIST SAME");
  hs.Draw("HIST SAME");
  hd.Draw("HIST SAME");
  h.Draw("E1 SAME");
  TLegend l(.57, .66, .87, .87);
  l.SetBorderSize(0);
  l.AddEntry(&h, "Nominal PP", "lep");
  l.AddEntry(&hp, "Template fit", "l");
  l.AddEntry(&hs, "SPS component", "l");
  l.AddEntry(&hd, "DPS component", "l");
  l.Draw();
  c.SaveAs(p.c_str());
}
} // namespace
void fit_1d_component_templates(
    const char *outdir = "Data_driven/results/"
                         "template_fraction_fits_adaptive14",
    const char *fitInput = "Data_driven/results/"
                           "component_fits_adaptive14/fit_results.csv") {
  gSystem->mkdir(outdir, true);
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
  auto m = load(vs, fitInput);
  std::string base(outdir);
  std::ofstream o(base + "/template_fit_results.csv"),
      l(base + "/leave_one_bin_out.csv");
  o << "variable,n_bins,sps_yield,sps_yield_error,dps_yield,dps_yield_error,"
       "yield_correlation,f_sps,f_sps_error,f_dps,chi2,ndf,chi2_ndf,accepted\n"
    << std::setprecision(12);
  l << "variable,excluded_bin,low,high,f_sps,f_sps_error,chi2,ndf\n"
    << std::setprecision(12);
  for (auto &v : vs) {
    F r = fit(m[v.n]);
    o << v.n << "," << m[v.n].size() << "," << r.a << "," << r.ea << "," << r.b
      << "," << r.eb << "," << r.corr << "," << r.f << "," << r.ef << ","
      << 1 - r.f << "," << r.chi2 << "," << r.ndf << ","
      << (r.ndf ? r.chi2 / r.ndf : NAN) << "," << r.ok << "\n";
    for (size_t i = 0; i < m[v.n].size(); i++) {
      F q = fit(m[v.n], i);
      l << v.n << "," << i << "," << v.e[i] << "," << v.e[i + 1] << "," << q.f
        << "," << q.ef << "," << q.chi2 << "," << q.ndf << "\n";
    }
    draw(v, m[v.n], r, base + "/" + v.n + "_template_fit.pdf");
  }
  std::cout << "fit_input=" << fitInput << "\noutput=" << outdir << std::endl;
}
