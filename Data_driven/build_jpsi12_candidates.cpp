#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TNamed.h"
#include "TSystem.h"
#include "TTree.h"

namespace {
const char *kDataBase =
    "/eos/user/c/chensh/JPsiJPsi/Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer";

void addNominalData(TChain &chain, const char *base) {
  struct Era { const char *name; int files; };
  const Era eras[] = {{"B",20},{"C",9},{"D",14},{"E",3},
                      {"F",8},{"G",29},{"H",36}};
  for (const auto &era : eras)
    for (int i=1; i<=era.files; ++i)
      chain.Add(Form("%s/%s/Ntuple_2016_%s_%d.root", base,
                     era.name, era.name, i));
}

void addDpsMc(TChain &chain) {
  const char *base =
      "Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/DPS_ntuple";
  for (int i=1; i<=65; ++i)
    chain.Add(Form("%s/Ntuple_2016_DPS_%d.root", base, i));
}
}

void build_jpsi12_candidates(
    const char *output=
        "Data_driven/results/route9p3_crossslot_splot/jpsi12_candidates.root",
    const char *inputMode="data",
    const char *dataBase=kDataBase) {
  TChain chain("rootuple/oniaTree");
  const std::string mode(inputMode);
  const int expectedFiles = mode=="data" ? 119 : mode=="dps_mc" ? 65 : -1;
  if (mode=="data") addNominalData(chain, dataBase);
  else if (mode=="dps_mc") addDpsMc(chain);
  else {
    std::cerr << "inputMode must be data or dps_mc\n";
    return;
  }
  if (chain.GetNtrees()!=expectedFiles) {
    std::cerr << "Expected " << expectedFiles << " input files, got "
              << chain.GetNtrees() << "\n";
    return;
  }
  ULong64_t run=0,lumi=0,event=0;
  std::vector<double> *pt=nullptr,*eta=nullptr,*y=nullptr,*phi=nullptr;
  std::vector<double> *mass=nullptr,*ctau=nullptr;
  std::vector<int> *id1=nullptr,*id2=nullptr;
  std::vector<bool> *passHLT=nullptr,*matchTrg=nullptr,*samePV=nullptr;
  chain.SetBranchAddress("run",&run); chain.SetBranchAddress("lumi",&lumi);
  chain.SetBranchAddress("event",&event); chain.SetBranchAddress("REJpsi_pt",&pt);
  chain.SetBranchAddress("REJpsi_eta",&eta); chain.SetBranchAddress("REJpsi_y",&y);
  chain.SetBranchAddress("REJpsi_phi",&phi); chain.SetBranchAddress("REJpsi_mass",&mass);
  chain.SetBranchAddress("REJpsi_ctau",&ctau);
  chain.SetBranchAddress("REevt_JpsiId1",&id1); chain.SetBranchAddress("REevt_JpsiId2",&id2);
  chain.SetBranchAddress("REevt_passHLT",&passHLT);
  chain.SetBranchAddress("REevt_matchTrg",&matchTrg);
  chain.SetBranchAddress("REevt_samePV",&samePV);

  TString outDir=gSystem->DirName(output); gSystem->mkdir(outDir,true);
  TFile out(output,"RECREATE");
  TTree t1("jpsi1","selected event Jpsi1 candidates, unweighted");
  TTree t2("jpsi2","selected event Jpsi2 candidates, unweighted");
  Long64_t source_entry=-1;
  int candidate_index1=-1,partner_index1=-1,candidate_index2=-1,partner_index2=-1;
  double m1=0,c1=0,pt1=0,eta1=0,y1=0,phi1=0;
  double m2=0,c2=0,pt2=0,eta2=0,y2=0,phi2=0;
  double original_pair_mass=0;
  for (TTree *t : {&t1,&t2}) {
    t->Branch("run",&run); t->Branch("lumi",&lumi); t->Branch("event",&event);
    t->Branch("source_entry",&source_entry);
    t->Branch("original_pair_mass",&original_pair_mass);
  }
  t1.Branch("candidate_index",&candidate_index1); t1.Branch("partner_index",&partner_index1);
  t1.Branch("Jpsi_mass1",&m1); t1.Branch("Jpsi_ctau1",&c1);
  t1.Branch("Jpsi_pt1",&pt1); t1.Branch("Jpsi_eta1",&eta1);
  t1.Branch("Jpsi_y1",&y1); t1.Branch("Jpsi_phi1",&phi1);
  t2.Branch("candidate_index",&candidate_index2); t2.Branch("partner_index",&partner_index2);
  t2.Branch("Jpsi_mass2",&m2); t2.Branch("Jpsi_ctau2",&c2);
  t2.Branch("Jpsi_pt2",&pt2); t2.Branch("Jpsi_eta2",&eta2);
  t2.Branch("Jpsi_y2",&y2); t2.Branch("Jpsi_phi2",&phi2);

  Long64_t selected=0;
  for (Long64_t entry=0; entry<chain.GetEntries(); ++entry) {
    chain.GetEntry(entry);
    for (int k=static_cast<int>(matchTrg->size())-1; k>=0; --k) {
      if (!passHLT->at(k)||!matchTrg->at(k)||!samePV->at(k)) continue;
      const int a=id1->at(k),b=id2->at(k);
      if (pt->at(a)<10||pt->at(a)>40||pt->at(b)<10||pt->at(b)>40) continue;
      TLorentzVector va,vb;
      va.SetPtEtaPhiM(pt->at(a),eta->at(a),phi->at(a),mass->at(a));
      vb.SetPtEtaPhiM(pt->at(b),eta->at(b),phi->at(b),mass->at(b));
      if ((va+vb).M()<7.5) continue;
      source_entry=entry; original_pair_mass=(va+vb).M();
      candidate_index1=a; partner_index1=b;
      m1=mass->at(a); c1=ctau->at(a); pt1=pt->at(a);
      eta1=eta->at(a); y1=y->at(a); phi1=phi->at(a);
      candidate_index2=b; partner_index2=a;
      m2=mass->at(b); c2=ctau->at(b); pt2=pt->at(b);
      eta2=eta->at(b); y2=y->at(b); phi2=phi->at(b);
      t1.Fill(); t2.Fill(); ++selected; break;
    }
  }
  std::ostringstream meta;
  meta << "input_mode="<<mode<<"\ninput_files="<<chain.GetNtrees()
       << "\ninput_base="<<(mode=="data" ? dataBase : "Data/DPS_ntuple")
       <<"\ninput_entries="<<chain.GetEntries()
       << "\nselected_events="<<selected
       << "\njpsi1_candidates="<<t1.GetEntries()
       << "\njpsi2_candidates="<<t2.GetEntries()
       << "\ninput_weight=none"
       << "\nselection=same_as_build_mixed_dps_with_mJJ_ge_7p5\n";
  TNamed metadata("run_metadata",meta.str().c_str());
  t1.Write(); t2.Write(); metadata.Write(); out.Close();
  std::cout << meta.str() << "output="<<output<<"\n";
}
