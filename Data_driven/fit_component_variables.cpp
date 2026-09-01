#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "TFile.h"
#include "TH1D.h"
#include "TSystem.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
using namespace RooFit;
namespace {
std::string input = "Data_driven/results/"
                    "component_trees_adaptive14/component_data.root";
const char *model =
    "Data_driven/inputs/Model_4D_tot.root";
struct V {
  std::string name, label;
  std::vector<double> e;
};
struct R {
  std::string var, component, weight;
  int bin = -1;
  double lo = 0, hi = 0;
  Long64_t raw = 0, negative = 0;
  int entries = 0, attempts = 0, status = 999, cov = -1;
  double sw = 0, sw2 = 0, negw = 0, pp = NAN, epp = NAN, edm = NAN, nll = NAN;
  bool ok = false;
};
std::string cut(const std::string &v, double a, double b, bool last) {
  std::ostringstream s;
  s << std::setprecision(17) << v << ">=" << a << "&&" << v
    << (last ? "<=" : "<") << b;
  return s.str();
}
void init(RooWorkspace &w, double scale, int a) {
  const char *n[] = {"n_P_P", "n_P_NP", "n_NP_NP", "n_Sig_Comb", "n_Comb_Comb"};
  double min[] = {1, 0, 0, 0, 0}, shift[] = {1, .8, 1.2, .6, 1.5, .4};
  for (int i = 0; i < 5; i++) {
    auto *x = w.var(n[i]);
    double g = x->getVal();
    x->setMin(min[i]);
    x->setMax(std::max(1000., g * std::max(2., 10 * std::abs(scale))));
    x->setVal(std::max(min[i] + 1e-3,
                       g * std::max(.001, std::abs(scale)) * shift[a % 6]));
    x->setConstant(false);
  }
}
R run(const V &v, int b, const std::string &comp, const std::string &weight) {
  R o;
  TFile inputFile(input.c_str());
  auto *tree = (TTree *)inputFile.Get("data_components");
  TTree &t = *tree;
  o.var = v.name;
  o.component = comp;
  o.weight = weight;
  o.bin = b;
  o.lo = v.e[b];
  o.hi = v.e[b + 1];
  std::string sel = cut(v.name, o.lo, o.hi, b + 2 == (int)v.e.size());
  o.raw = t.GetEntries(sel.c_str());
  RooRealVar m1("Jpsi_mass1", "Jpsi_mass1", 2.95, 3.25),
      m2("Jpsi_mass2", "Jpsi_mass2", 2.95, 3.25),
      c1("Jpsi_ctau1", "Jpsi_ctau1", -.03, .16),
      c2("Jpsi_ctau2", "Jpsi_ctau2", -.03, .16),
      w(weight.c_str(), weight.c_str(), comp == "sps" ? -1000.0 : 0.0, 1000.0),
      x(v.name.c_str(), v.name.c_str(), v.e.front(), v.e.back());
  RooArgSet all(m1, m2, c1, c2, w, x);
  RooDataSet data("data", "", all, Import(t), Cut(sel.c_str()),
                  WeightVar(weight.c_str()));
  std::unique_ptr<RooAbsData> d(data.reduce(RooArgSet(m1, m2, c1, c2)));
  o.entries = d ? d->numEntries() : 0;
  o.sw = d ? d->sumEntries() : 0;
  TTreeFormula f("f", sel.c_str(), &t);
  TTreeFormula wf("wf", weight.c_str(), &t);
  double ww = 0;
  for (Long64_t i = 0; i < t.GetEntries(); i++) {
    t.GetEntry(i);
    ww = wf.EvalInstance();
    if (f.EvalInstance()) {
      o.sw2 += ww * ww;
      if (ww < 0) {
        o.negative++;
        o.negw += ww;
      }
    }
  }
  if (!d || o.entries < 10 || o.sw <= 0)
    return o;
  TFile mf(model);
  auto *wsp = (RooWorkspace *)mf.Get("wsp");
  auto *pdf = wsp ? wsp->pdf("pdf_all") : nullptr;
  if (!pdf)
    return o;
  double global =
      wsp->var("n_P_P")->getVal() + 2 * wsp->var("n_P_NP")->getVal() +
      wsp->var("n_NP_NP")->getVal() + 2 * wsp->var("n_Sig_Comb")->getVal() +
      wsp->var("n_Comb_Comb")->getVal();
  double scale = global ? o.sw / global : .05, best = 1e99;
  for (int a = 0; a < 6; a++) {
    init(*wsp, scale, a);
    std::unique_ptr<RooFitResult> z(pdf->fitTo(
        *d, Extended(true), Save(true), SumW2Error(false), Offset(true),
        Optimize(false), PrintLevel(-1), Strategy(1),
        Minimizer("Minuit2", "migrad")));
    o.attempts++;
    if (!z)
      continue;
    double score = 1000 * std::abs(z->status()) +
                   100 * std::max(0, 3 - z->covQual()) +
                   std::min(100., std::abs(z->edm()));
    if (score < best) {
      best = score;
      o.status = z->status();
      o.cov = z->covQual();
      o.edm = z->edm();
      o.nll = z->minNll();
      auto *p = wsp->var("n_P_P");
      o.pp = p->getVal();
      o.epp = p->getError();
    }
    if (z->status() == 0 && z->covQual() >= 2 && z->edm() < .01) {
      o.ok = true;
      break;
    }
  }
  return o;
}
} // namespace
void fit_component_variables(
    const char *onlyVar = "", int onlyBin = -1, const char *onlyComp = "",
    const char *outdir =
        "Data_driven/results/component_fits_adaptive14",
    const char *inputPath =
        "Data_driven/results/component_trees_adaptive14/"
        "component_data.root") {
  input = inputPath;
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  gSystem->mkdir(outdir, true);
  const double pi = 3.14159265358979323846;
  std::vector<V> vars = {
      {"delta_y", "|#Delta y|", {0, .5, 1, 1.5, 2, 2.5, 4}},
      {"delta_phi",
       "|#Delta#phi|",
       {0, .3927, .7854, 1.1781, 1.5708, 1.9635, 2.3562, 2.7489, 3.1416}},
      {"evt_mass",
       "m(J/#psi J/#psi)",
       {7.5, 17.5, 27.5, 37.5, 47.5, 57.5, 67.5, 107.5}},
      {"evt_y", "|y(J/#psi J/#psi)|", {0, .4, .8, 1.2, 1.6, 2}}};
  std::vector<std::pair<std::string, std::string>> comps = {
      {"total", "evt_weight"},
      {"sps", "evt_weight_sps"},
      {"dps", "evt_weight_dps"}};
  TFile inputCheck(input.c_str());
  if (!inputCheck.Get("data_components")) {
    std::cerr << "missing data_components\n";
    return;
  }
  inputCheck.Close();
  std::vector<R> res;
  for (auto &v : vars) {
    if (std::string(onlyVar).size() && v.name != onlyVar)
      continue;
    for (int b = 0; b + 1 < (int)v.e.size(); b++) {
      if (onlyBin >= 0 && b != onlyBin)
        continue;
      for (auto &c : comps) {
        if (std::string(onlyComp).size() && c.first != onlyComp)
          continue;
        std::cout << v.name << " [" << v.e[b] << "," << v.e[b + 1] << "] "
                  << c.first << std::endl;
        R r = run(v, b, c.first, c.second);
        std::cout << " raw=" << r.raw << " entries=" << r.entries
                  << " sw=" << r.sw << " neg=" << r.negative << " PP=" << r.pp
                  << " +/- " << r.epp << " status=" << r.status
                  << " cov=" << r.cov << " edm=" << r.edm << " ok=" << r.ok
                  << std::endl;
        res.push_back(r);
      }
    }
  }
  std::string base(outdir);
  std::ofstream csv(base + "/fit_results.csv");
  csv << "variable,component,weight,bin,low,high,raw_entries,fit_entries,"
         "negative_weight_entries,sum_weights,sum_weights2,sum_negative_"
         "weights,pp_yield,pp_error,status,covQual,edm,minNll,attempts,"
         "accepted\n"
      << std::setprecision(12);
  for (auto &r : res)
    csv << r.var << "," << r.component << "," << r.weight << "," << r.bin << ","
        << r.lo << "," << r.hi << "," << r.raw << "," << r.entries << ","
        << r.negative << "," << r.sw << "," << r.sw2 << "," << r.negw << ","
        << r.pp << "," << r.epp << "," << r.status << "," << r.cov << ","
        << r.edm << "," << r.nll << "," << r.attempts << "," << r.ok << "\n";
  if (std::string(onlyVar).size()) {
    csv.close();
    std::cout << "single_fit_result_written=" << base << "/fit_results.csv"
              << std::endl;
    std::cout.flush();
    _exit(0);
  }
  TFile fo((base + "/component_yields.root").c_str(), "RECREATE");
  int accepted = 0, failed = 0;
  std::ofstream close(base + "/closure.csv");
  close << "variable,bin,low,high,total,total_error,sps,sps_error,dps,dps_"
           "error,sum_components,difference,relative_difference\n"
        << std::setprecision(12);
  double maxrel = 0;
  for (auto &v : vars) {
    int nb = v.e.size() - 1;
    TH1D ht(("h_" + v.name + "_total").c_str(),
            (";" + v.label + ";PP corrected yield").c_str(), nb, v.e.data()),
        hs(("h_" + v.name + "_sps").c_str(),
           (";" + v.label + ";PP corrected yield").c_str(), nb, v.e.data()),
        hd(("h_" + v.name + "_dps").c_str(),
           (";" + v.label + ";PP corrected yield").c_str(), nb, v.e.data());
    for (int b = 0; b < nb; b++) {
      R *rt = nullptr, *rs = nullptr, *rd = nullptr;
      for (auto &r : res)
        if (r.var == v.name && r.bin == b) {
          if (r.component == "total")
            rt = &r;
          if (r.component == "sps")
            rs = &r;
          if (r.component == "dps")
            rd = &r;
        }
      for (R *q : {rt, rs, rd}) {
        if (q && q->ok)
          accepted++;
        else
          failed++;
      }
      if (rt && rs && rd) {
        ht.SetBinContent(b + 1, rt->pp);
        ht.SetBinError(b + 1, rt->epp);
        hs.SetBinContent(b + 1, rs->pp);
        hs.SetBinError(b + 1, rs->epp);
        hd.SetBinContent(b + 1, rd->pp);
        hd.SetBinError(b + 1, rd->epp);
        double sum = rs->pp + rd->pp, diff = sum - rt->pp,
               rel = rt->pp ? diff / rt->pp : NAN;
        if (std::isfinite(rel))
          maxrel = std::max(maxrel, std::abs(rel));
        close << v.name << "," << b << "," << v.e[b] << "," << v.e[b + 1] << ","
              << rt->pp << "," << rt->epp << "," << rs->pp << "," << rs->epp
              << "," << rd->pp << "," << rd->epp << "," << sum << "," << diff
              << "," << rel << "\n";
      }
    }
    ht.Write();
    hs.Write();
    hd.Write();
  }
  fo.Close();
  std::ofstream s(base + "/summary.txt");
  s << "input=" << input << "\nmodel=" << model
    << "\nvariables=4\none_dimensional_bins=26\nfits=78\naccepted=" << accepted
    << "\nfailed=" << failed
    << "\nacceptance=status==0 && covQual>=2 && "
       "edm<0.01\nsumw2_error=false\nlikelihood_offset=true\nmax_abs_relative_"
       "sps_plus_dps_minus_total="
    << std::setprecision(12) << maxrel
    << "\nwarning=component_fraction_uncertainty_not_included\nwarning=sps_"
       "dataset_contains_small_negative_weights_from_control_region\n";
  std::cout << "accepted=" << accepted << " failed=" << failed
            << " max closure=" << maxrel << "\noutput=" << outdir << std::endl;
}
