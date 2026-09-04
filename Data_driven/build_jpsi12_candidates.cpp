#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TLorentzVector.h"
#include "TNamed.h"
#include "TSystem.h"
#include "TTree.h"

namespace {
const char *kDataBase =
    "/eos/user/c/chensh/JPsiJPsi/Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer";

void addNominalData(std::vector<std::string> &files, const char *base) {
  struct Era { const char *name; int files; };
  const Era eras[] = {{"B",20},{"C",9},{"D",14},{"E",3},
                      {"F",8},{"G",29},{"H",36}};
  for (const auto &era : eras)
    for (int i=1; i<=era.files; ++i)
      files.emplace_back(Form("%s/%s/Ntuple_2016_%s_%d.root", base,
                              era.name, era.name, i));
}

void addDpsMc(std::vector<std::string> &files) {
  const char *base =
      "Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/DPS_ntuple";
  for (int i=1; i<=65; ++i)
    files.emplace_back(Form("%s/Ntuple_2016_DPS_%d.root", base, i));
}

int addManifest(std::vector<std::string> &files, const char *manifestPath) {
  std::ifstream manifest(manifestPath);
  if (!manifest) {
    std::cerr << "Cannot open input manifest " << manifestPath << "\n";
    return -1;
  }
  int count = 0;
  std::string path;
  while (std::getline(manifest, path)) {
    if (path.empty() || path[0] == '#') continue;
    files.push_back(path);
    ++count;
  }
  return count;
}
}

void build_jpsi12_candidates(
    const char *output=
        "Data_driven/results/route9p3_crossslot_splot/jpsi12_candidates.root",
    const char *inputMode="data",
    const char *dataBase=kDataBase,
    const char *inputManifest="",
    bool symmetrizeJpsiLabels=false,
    unsigned int symmetrizationSeed=50) {
  std::vector<std::string> inputFiles;
  const std::string mode(inputMode);
  const bool useManifest = inputManifest && inputManifest[0];
  const int manifestFiles = useManifest ? addManifest(inputFiles, inputManifest) : 0;
  const int expectedFiles = mode=="data" ? (useManifest ? manifestFiles : 119)
                                          : mode=="dps_mc"
                                                ? (useManifest ? manifestFiles : 65)
                                                : -1;
  if (!useManifest && mode=="data") addNominalData(inputFiles, dataBase);
  else if (!useManifest && mode=="dps_mc") addDpsMc(inputFiles);
  else if (mode!="data" && mode!="dps_mc") {
    std::cerr << "inputMode must be data or dps_mc\n";
    return;
  }
  if (static_cast<int>(inputFiles.size())!=expectedFiles) {
    std::cerr << "Expected " << expectedFiles << " input files, got "
              << inputFiles.size() << "\n";
    return;
  }
  ULong64_t run=0,lumi=0,event=0;

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

  Long64_t selected=0, randomizationDraws=0, randomizedCandidateSwaps=0;
  Long64_t randomizedSelectedSwaps=0;
  Long64_t inputEntries=0;
  int reseededFiles=0;
  for (const std::string &inputPath : inputFiles) {
    TFile input(inputPath.c_str(), "READ");
    if (input.IsZombie()) {
      std::cerr << "Cannot open input ROOT file " << inputPath << "\n";
      return;
    }
    TTree *tree=nullptr;
    input.GetObject("rootuple/oniaTree", tree);
    if (!tree) {
      std::cerr << "Missing rootuple/oniaTree in " << inputPath << "\n";
      return;
    }
    std::vector<double> *pt=nullptr,*eta=nullptr,*y=nullptr,*phi=nullptr;
    std::vector<double> *mass=nullptr,*ctau=nullptr;
    std::vector<int> *id1=nullptr,*id2=nullptr;
    std::vector<bool> *passHLT=nullptr,*matchTrg=nullptr,*samePV=nullptr;
    tree->SetBranchAddress("run",&run); tree->SetBranchAddress("lumi",&lumi);
    tree->SetBranchAddress("event",&event); tree->SetBranchAddress("REJpsi_pt",&pt);
    tree->SetBranchAddress("REJpsi_eta",&eta); tree->SetBranchAddress("REJpsi_y",&y);
    tree->SetBranchAddress("REJpsi_phi",&phi); tree->SetBranchAddress("REJpsi_mass",&mass);
    tree->SetBranchAddress("REJpsi_ctau",&ctau);
    tree->SetBranchAddress("REevt_JpsiId1",&id1); tree->SetBranchAddress("REevt_JpsiId2",&id2);
    tree->SetBranchAddress("REevt_passHLT",&passHLT);
    tree->SetBranchAddress("REevt_matchTrg",&matchTrg);
    tree->SetBranchAddress("REevt_samePV",&samePV);
    if (symmetrizeJpsiLabels) {
      std::srand(symmetrizationSeed);
      ++reseededFiles;
    }
    const Long64_t fileEntries=tree->GetEntries();
    for (Long64_t localEntry=0; localEntry<fileEntries; ++localEntry) {
      tree->GetEntry(localEntry);
      for (int k=static_cast<int>(matchTrg->size())-1; k>=0; --k) {
      if (!passHLT->at(k)||!matchTrg->at(k)) continue;
      int a=id1->at(k),b=id2->at(k);
      bool labelsSwapped=false;
      if (symmetrizeJpsiLabels) {
        ++randomizationDraws;
        if (static_cast<double>(std::rand())/(RAND_MAX+1.0)<=0.5) {
          std::swap(a,b);
          labelsSwapped=true;
          ++randomizedCandidateSwaps;
        }
      }
      if (pt->at(a)<10||pt->at(a)>40||pt->at(b)<10||pt->at(b)>40) continue;
      if (!samePV->at(k)) continue;
      TLorentzVector va,vb;
      va.SetPtEtaPhiM(pt->at(a),eta->at(a),phi->at(a),mass->at(a));
      vb.SetPtEtaPhiM(pt->at(b),eta->at(b),phi->at(b),mass->at(b));
      if ((va+vb).M()<7.5) continue;
      source_entry=inputEntries+localEntry; original_pair_mass=(va+vb).M();
      candidate_index1=a; partner_index1=b;
      m1=mass->at(a); c1=ctau->at(a); pt1=pt->at(a);
      eta1=eta->at(a); y1=y->at(a); phi1=phi->at(a);
      candidate_index2=b; partner_index2=a;
      m2=mass->at(b); c2=ctau->at(b); pt2=pt->at(b);
      eta2=eta->at(b); y2=y->at(b); phi2=phi->at(b);
      if (labelsSwapped) ++randomizedSelectedSwaps;
      t1.Fill(); t2.Fill(); ++selected; break;
      }
    }
    inputEntries+=fileEntries;
  }
  std::ostringstream meta;
  meta << "input_mode="<<mode<<"\ninput_files="<<inputFiles.size()
       << "\ninput_base="<<(mode=="data" ? dataBase : "Data/DPS_ntuple")
       << "\ninput_manifest="<<(useManifest ? inputManifest : "none")
       <<"\ninput_entries="<<inputEntries
       << "\nselected_events="<<selected
       << "\njpsi1_candidates="<<t1.GetEntries()
       << "\njpsi2_candidates="<<t2.GetEntries()
       << "\ninput_weight=none"
       << "\nlabel_randomization="<<(symmetrizeJpsiLabels ? "enabled" : "disabled")
       << "\nlabel_randomization_seed="<<symmetrizationSeed
       << "\nlabel_randomization_reseeded_files="<<reseededFiles
       << "\nrandomization_draws="<<randomizationDraws
       << "\nrandomized_candidate_swaps="<<randomizedCandidateSwaps
       << "\nrandomized_selected_event_swaps="<<randomizedSelectedSwaps
       << "\ninput_iteration=direct_tfile_rephrase_compatible"
       << "\nselection=same_as_build_mixed_dps_with_mJJ_ge_7p5\n";
  TNamed metadata("run_metadata",meta.str().c_str());
  out.cd();
  t1.Write(); t2.Write(); metadata.Write(); out.Close();
  std::cout << meta.str() << "output="<<output<<"\n";
}
