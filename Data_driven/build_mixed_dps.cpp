#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TNamed.h"
#include "TSystem.h"
#include "TTree.h"

namespace {

constexpr double kPi = 3.14159265358979323846;
const char *kDataBase =
    "/eos/user/c/chensh/JPsiJPsi/Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer";
const char *kAccFile =
    "Data_driven/inputs/acceptance_sps_full10_v1.txt";
const char *kEffFile =
    "Data_driven/inputs/efficiency_sps0p8_dps0p2_dedup60_v1.txt";

struct Jpsi {
    double pt = 0;
    double eta = 0;
    double y = 0;
    double phi = 0;
    double mass = 0;
};

struct PoolEvent {
    ULong64_t run = 0;
    ULong64_t lumi = 0;
    ULong64_t event = 0;
    Jpsi jpsi[2];
};

bool sameEvent(const PoolEvent &a, const PoolEvent &b) {
    return a.run == b.run && a.lumi == b.lumi && a.event == b.event;
}

int findBin(const std::vector<double> &edges, double value) {
    if (edges.size() < 2 || value < edges.front() || value > edges.back())
        return -1;
    if (value == edges.back())
        return static_cast<int>(edges.size()) - 2;
    return static_cast<int>(
               std::upper_bound(edges.begin(), edges.end(), value) -
               edges.begin()) -
           1;
}

class Correction {
  public:
    bool load() {
        return loadAcceptance(kAccFile) && loadEfficiency(kEffFile);
    }

    double weight(const Jpsi &a, const Jpsi &b) const {
        const int ay1 = findBin(accY_, a.y);
        const int ap1 = findBin(accPt_, a.pt);
        const int ay2 = findBin(accY_, b.y);
        const int ap2 = findBin(accPt_, b.pt);
        const int ey1 = findBin(effY_, a.y);
        const int ep1 = findBin(effPt_, a.pt);
        const int ey2 = findBin(effY_, b.y);
        const int ep2 = findBin(effPt_, b.pt);
        if (ay1 < 0 || ap1 < 0 || ay2 < 0 || ap2 < 0 || ey1 < 0 ||
            ep1 < 0 || ey2 < 0 || ep2 < 0)
            return -1;
        const double nAcc1 = nAcc_[ay1][ap1];
        const double nAcc2 = nAcc_[ay2][ap2];
        const double nVtx1 = nVtxJpsi_[ey1][ep1];
        const double nVtx2 = nVtxJpsi_[ey2][ep2];
        const double nTrg = nTrgEvt_[ep2][ep1];
        if (nAcc1 <= 0 || nAcc2 <= 0 || nVtx1 <= 0 || nVtx2 <= 0 ||
            nTrg <= 0)
            return -1;
        return nGen_[ay1][ap1] / nAcc1 * nGen_[ay2][ap2] / nAcc2 *
               nBinJpsi_[ey1][ep1] / nVtx1 *
               nBinJpsi_[ey2][ep2] / nVtx2 *
               nVtxEvt_[ep2][ep1] / nTrg;
    }

  private:
    std::vector<double> accPt_, accY_, effPt_, effY_;
    std::vector<std::vector<double>> nGen_, nAcc_, nBinJpsi_, nVtxJpsi_;
    std::vector<std::vector<double>> nVtxEvt_, nTrgEvt_;

    bool loadAcceptance(const std::string &path) {
        std::ifstream in(path);
        int nPt = 0, nY = 0;
        double ignored1 = 0, ignored2 = 0;
        if (!(in >> nPt >> nY >> ignored1 >> ignored2))
            return false;
        accPt_.resize(nPt + 1);
        accY_.resize(nY + 1);
        for (double &x : accPt_)
            in >> x;
        for (double &x : accY_)
            in >> x;
        nAcc_.assign(nY, std::vector<double>(nPt));
        nGen_.assign(nY, std::vector<double>(nPt));
        for (int iy = 0; iy < nY; ++iy)
            for (int ip = 0; ip < nPt; ++ip)
                in >> nAcc_[iy][ip] >> nGen_[iy][ip];
        return static_cast<bool>(in);
    }

    bool loadEfficiency(const std::string &path) {
        std::ifstream in(path);
        int nPt = 0, nY = 0;
        if (!(in >> nPt >> nY))
            return false;
        effPt_.resize(nPt + 1);
        effY_.resize(nY + 1);
        for (double &x : effPt_)
            in >> x;
        for (double &x : effY_)
            in >> x;
        nVtxJpsi_.assign(nY, std::vector<double>(nPt));
        nBinJpsi_.assign(nY, std::vector<double>(nPt));
        for (int iy = 0; iy < nY; ++iy)
            for (int ip = 0; ip < nPt; ++ip)
                in >> nVtxJpsi_[iy][ip] >> nBinJpsi_[iy][ip];
        nTrgEvt_.assign(nPt, std::vector<double>(nPt));
        nVtxEvt_.assign(nPt, std::vector<double>(nPt));
        for (int iy = 0; iy < nPt; ++iy)
            for (int ip = 0; ip < nPt; ++ip)
                in >> nTrgEvt_[iy][ip] >> nVtxEvt_[iy][ip];
        return static_cast<bool>(in);
    }
};

void addNominalData(TChain &chain) {
    struct Era {
        const char *name;
        int files;
    };
    const Era eras[] = {{"B", 20}, {"C", 9},  {"D", 14}, {"E", 3},
                        {"F", 8},  {"G", 29}, {"H", 36}};
    for (const Era &era : eras)
        for (int i = 1; i <= era.files; ++i)
            chain.Add(Form("%s/%s/Ntuple_2016_%s_%d.root", kDataBase,
                           era.name, era.name, i));
}

std::vector<PoolEvent> buildPool(TChain &chain) {
    ULong64_t run = 0, lumi = 0, event = 0;
    std::vector<double> *pt = nullptr, *eta = nullptr, *y = nullptr;
    std::vector<double> *phi = nullptr, *mass = nullptr;
    std::vector<int> *id1 = nullptr, *id2 = nullptr;
    std::vector<bool> *passHLT = nullptr, *matchTrg = nullptr;
    std::vector<bool> *samePV = nullptr;

    chain.SetBranchAddress("run", &run);
    chain.SetBranchAddress("lumi", &lumi);
    chain.SetBranchAddress("event", &event);
    chain.SetBranchAddress("REJpsi_pt", &pt);
    chain.SetBranchAddress("REJpsi_eta", &eta);
    chain.SetBranchAddress("REJpsi_y", &y);
    chain.SetBranchAddress("REJpsi_phi", &phi);
    chain.SetBranchAddress("REJpsi_mass", &mass);
    chain.SetBranchAddress("REevt_JpsiId1", &id1);
    chain.SetBranchAddress("REevt_JpsiId2", &id2);
    chain.SetBranchAddress("REevt_passHLT", &passHLT);
    chain.SetBranchAddress("REevt_matchTrg", &matchTrg);
    chain.SetBranchAddress("REevt_samePV", &samePV);

    std::vector<PoolEvent> pool;
    for (Long64_t entry = 0; entry < chain.GetEntries(); ++entry) {
        chain.GetEntry(entry);
        for (int k = static_cast<int>(matchTrg->size()) - 1; k >= 0; --k) {
            if (!passHLT->at(k) || !matchTrg->at(k) || !samePV->at(k))
                continue;
            const int a = id1->at(k), b = id2->at(k);
            if (pt->at(a) < 10 || pt->at(a) > 40 || pt->at(b) < 10 ||
                pt->at(b) > 40)
                continue;
            TLorentzVector va, vb;
            va.SetPtEtaPhiM(pt->at(a), eta->at(a), phi->at(a), mass->at(a));
            vb.SetPtEtaPhiM(pt->at(b), eta->at(b), phi->at(b), mass->at(b));
            if ((va + vb).M() < 7.5)
                continue;
            PoolEvent out;
            out.run = run;
            out.lumi = lumi;
            out.event = event;
            const int ids[2] = {a, b};
            for (int j = 0; j < 2; ++j) {
                const int q = ids[j];
                out.jpsi[j] = {pt->at(q), eta->at(q), y->at(q),
                               phi->at(q), mass->at(q)};
            }
            pool.push_back(out);
            break;
        }
    }
    return pool;
}

} // namespace

void build_mixed_dps(Long64_t requested = 100000,
                     ULong64_t seed = 20260728,
                     const char *output =
                         "Data_driven/results/mixed_dps.root") {
    TChain chain("rootuple/oniaTree");
    addNominalData(chain);
    if (chain.GetNtrees() != 119) {
        std::cerr << "Expected 119 input files, got " << chain.GetNtrees()
                  << std::endl;
        return;
    }
    const std::vector<PoolEvent> pool = buildPool(chain);
    if (pool.size() < 2) {
        std::cerr << "Selected pool is too small: " << pool.size() << std::endl;
        return;
    }
    Correction correction;
    if (!correction.load()) {
        std::cerr << "Cannot read nominal correction tables." << std::endl;
        return;
    }

    TString outPath(output);
    gSystem->mkdir(gSystem->DirName(outPath), true);
    TFile out(output, "RECREATE");
    TTree tree("mix", "cross-event mixed J/psi pairs");
    ULong64_t run1, lumi1, event1, run2, lumi2, event2;
    double jpsi_pt1, jpsi_y1, jpsi_pt2, jpsi_y2;
    double evt_mass, evt_y, delta_y, delta_phi, evt_weight;
    tree.Branch("run1", &run1);
    tree.Branch("lumi1", &lumi1);
    tree.Branch("event1", &event1);
    tree.Branch("run2", &run2);
    tree.Branch("lumi2", &lumi2);
    tree.Branch("event2", &event2);
    tree.Branch("Jpsi_pt1", &jpsi_pt1);
    tree.Branch("Jpsi_y1", &jpsi_y1);
    tree.Branch("Jpsi_pt2", &jpsi_pt2);
    tree.Branch("Jpsi_y2", &jpsi_y2);
    tree.Branch("evt_mass", &evt_mass);
    tree.Branch("evt_y", &evt_y);
    tree.Branch("delta_y", &delta_y);
    tree.Branch("delta_phi", &delta_phi);
    tree.Branch("evt_weight", &evt_weight);

    TH2D hDyDphi("h_mix_dy_dphi", ";|#Delta y|;|#Delta#phi|",
                 8, 0, 4, 8, 0, kPi);
    TH1D hMass("h_mix_mass", ";m(J/#psi J/#psi) [GeV];Weighted pairs",
               20, 7.5, 107.5);
    TH1D hY("h_mix_y", ";|y(J/#psi J/#psi)|;Weighted pairs", 10, 0, 2);
    hDyDphi.Sumw2();
    hMass.Sumw2();
    hY.Sumw2();

    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<size_t> eventPick(0, pool.size() - 1);
    std::uniform_int_distribution<int> jpsiPick(0, 1);
    Long64_t attempts = 0, sameEventRejected = 0, massRejected = 0;
    Long64_t weightRejected = 0;
    while (tree.GetEntries() < requested) {
        ++attempts;
        const PoolEvent &a = pool[eventPick(rng)];
        const PoolEvent &b = pool[eventPick(rng)];
        if (sameEvent(a, b)) {
            ++sameEventRejected;
            continue;
        }
        const Jpsi &ja = a.jpsi[jpsiPick(rng)];
        const Jpsi &jb = b.jpsi[jpsiPick(rng)];
        TLorentzVector va, vb;
        va.SetPtEtaPhiM(ja.pt, ja.eta, ja.phi, ja.mass);
        vb.SetPtEtaPhiM(jb.pt, jb.eta, jb.phi, jb.mass);
        const TLorentzVector pair = va + vb;
        if (pair.M() < 7.5) {
            ++massRejected;
            continue;
        }
        const double weight = correction.weight(ja, jb);
        if (!(weight > 0) || !std::isfinite(weight)) {
            ++weightRejected;
            continue;
        }
        run1 = a.run;
        lumi1 = a.lumi;
        event1 = a.event;
        run2 = b.run;
        lumi2 = b.lumi;
        event2 = b.event;
        jpsi_pt1 = ja.pt;
        jpsi_y1 = ja.y;
        jpsi_pt2 = jb.pt;
        jpsi_y2 = jb.y;
        evt_mass = pair.M();
        evt_y = std::abs(pair.Rapidity());
        delta_y = std::abs(ja.y - jb.y);
        delta_phi = kPi - std::abs(std::abs(ja.phi - jb.phi) - kPi);
        evt_weight = weight;
        tree.Fill();
        hDyDphi.Fill(delta_y, delta_phi, evt_weight);
        hMass.Fill(evt_mass, evt_weight);
        hY.Fill(evt_y, evt_weight);
    }

    std::ostringstream meta;
    meta << "requested=" << requested << "\n"
         << "seed=" << seed << "\n"
         << "input_files=119\n"
         << "input_entries=" << chain.GetEntries() << "\n"
         << "selected_pool_events=" << pool.size() << "\n"
         << "attempts=" << attempts << "\n"
         << "same_event_rejected=" << sameEventRejected << "\n"
         << "mass_rejected=" << massRejected << "\n"
         << "weight_rejected=" << weightRejected << "\n";
    TNamed metadata("run_metadata", meta.str().c_str());
    tree.Write();
    hDyDphi.Write();
    hMass.Write();
    hY.Write();
    metadata.Write();
    out.Close();
    std::cout << meta.str();
    std::cout << "output=" << output << std::endl;
}
