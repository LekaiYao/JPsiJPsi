#include <iomanip>

#include "build_mixed_dps.cpp"

namespace {
struct SlotCandidate {
  ULong64_t run=0,lumi=0,event=0;
  Long64_t sourceEntry=-1;
  Jpsi jpsi;
  double sw=0;
};

std::vector<SlotCandidate> loadSlot(const char *path,int slot) {
  TFile f(path,"READ"); auto *t=dynamic_cast<TTree *>(f.Get("candidates"));
  std::vector<SlotCandidate> pool;
  if (!t) return pool;
  ULong64_t run=0,lumi=0,event=0; Long64_t sourceEntry=-1;
  double pt=0,eta=0,y=0,phi=0,mass=0,sw=0;
  const std::string s=std::to_string(slot);
  t->SetBranchAddress("run",&run); t->SetBranchAddress("lumi",&lumi);
  t->SetBranchAddress("event",&event); t->SetBranchAddress(("Jpsi_pt"+s).c_str(),&pt);
  if (t->GetBranch("source_entry"))
    t->SetBranchAddress("source_entry",&sourceEntry);
  t->SetBranchAddress(("Jpsi_eta"+s).c_str(),&eta); t->SetBranchAddress(("Jpsi_y"+s).c_str(),&y);
  t->SetBranchAddress(("Jpsi_phi"+s).c_str(),&phi); t->SetBranchAddress(("Jpsi_mass"+s).c_str(),&mass);
  t->SetBranchAddress("prompt_sweight",&sw);
  for (Long64_t i=0;i<t->GetEntries();++i) {
    t->GetEntry(i); if (std::isfinite(sw))
      pool.push_back({run,lumi,event,sourceEntry,{pt,eta,y,phi,mass},sw});
  }
  return pool;
}
}

void build_crossslot_splot_mixed(
    Long64_t requested=0, ULong64_t seed=20260807,
    const char *input1=
        "Data_driven/results/route9p3_crossslot_splot/jpsi1_sweights.root",
    const char *input2=
        "Data_driven/results/route9p3_crossslot_splot/jpsi2_sweights.root",
    const char *output=
        "Data_driven/results/route9p3_crossslot_allpairs/mixed_dps_allpairs.root") {
  const auto pool1=loadSlot(input1,1),pool2=loadSlot(input2,2);
  if (pool1.empty()||pool2.empty()) { std::cerr << "Missing slot candidate pool\n"; return; }
  Correction correction; if (!correction.load()) { std::cerr << "Cannot load correction tables\n"; return; }
  TString outDir=gSystem->DirName(output); gSystem->mkdir(outDir,true);
  TFile fo(output,"RECREATE"); TTree tree("mix","cross-event Jpsi1 x Jpsi2 prompt-sWeighted pairs");
  ULong64_t run1,lumi1,event1,run2,lumi2,event2;
  double jpsi_pt1,jpsi_y1,jpsi_pt2,jpsi_y2,evt_mass,evt_y,delta_y,delta_phi;
  double evt_weight,sweight_product,correction_weight;
  tree.Branch("run1",&run1); tree.Branch("lumi1",&lumi1); tree.Branch("event1",&event1);
  tree.Branch("run2",&run2); tree.Branch("lumi2",&lumi2); tree.Branch("event2",&event2);
  tree.Branch("Jpsi_pt1",&jpsi_pt1); tree.Branch("Jpsi_y1",&jpsi_y1);
  tree.Branch("Jpsi_pt2",&jpsi_pt2); tree.Branch("Jpsi_y2",&jpsi_y2);
  tree.Branch("evt_mass",&evt_mass); tree.Branch("evt_y",&evt_y);
  tree.Branch("delta_y",&delta_y); tree.Branch("delta_phi",&delta_phi);
  tree.Branch("evt_weight",&evt_weight); tree.Branch("sweight_product",&sweight_product);
  tree.Branch("correction_weight",&correction_weight);
  Long64_t attempts=0,sameRejected=0,massRejected=0,correctionRejected=0;
  Long64_t positive=0,negative=0,zero=0; double sumw=0,sumw2=0,maxAbs=0;
  auto fillPair = [&](const SlotCandidate &a, const SlotCandidate &b) {
    ++attempts;
    const bool sameSource = a.sourceEntry>=0 && b.sourceEntry>=0
        ? a.sourceEntry==b.sourceEntry
        : a.run==b.run && a.lumi==b.lumi && a.event==b.event;
    if (sameSource) { ++sameRejected; return; }
    TLorentzVector va,vb; va.SetPtEtaPhiM(a.jpsi.pt,a.jpsi.eta,a.jpsi.phi,a.jpsi.mass);
    vb.SetPtEtaPhiM(b.jpsi.pt,b.jpsi.eta,b.jpsi.phi,b.jpsi.mass);
    const TLorentzVector pair=va+vb; if (pair.M()<7.5) { ++massRejected; return; }
    correction_weight=correction.weight(a.jpsi,b.jpsi);
    if (!(correction_weight>0)||!std::isfinite(correction_weight)) { ++correctionRejected; return; }
    sweight_product=a.sw*b.sw; evt_weight=sweight_product*correction_weight;
    if (!std::isfinite(evt_weight)) { ++correctionRejected; return; }
    run1=a.run;lumi1=a.lumi;event1=a.event;run2=b.run;lumi2=b.lumi;event2=b.event;
    jpsi_pt1=a.jpsi.pt;jpsi_y1=a.jpsi.y;jpsi_pt2=b.jpsi.pt;jpsi_y2=b.jpsi.y;
    evt_mass=pair.M();evt_y=std::abs(pair.Rapidity());delta_y=std::abs(a.jpsi.y-b.jpsi.y);
    delta_phi=kPi-std::abs(std::abs(a.jpsi.phi-b.jpsi.phi)-kPi);
    if (evt_weight>0) ++positive; else if (evt_weight<0) ++negative; else ++zero;
    sumw+=evt_weight;sumw2+=evt_weight*evt_weight;maxAbs=std::max(maxAbs,std::abs(evt_weight));tree.Fill();
  };
  if (requested<=0) {
    for (const auto &a : pool1)
      for (const auto &b : pool2)
        fillPair(a,b);
  } else {
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<size_t> pick1(0,pool1.size()-1),pick2(0,pool2.size()-1);
    while (tree.GetEntries()<requested) fillPair(pool1[pick1(rng)],pool2[pick2(rng)]);
  }
  std::ostringstream meta;
  meta << std::setprecision(12)<<"mixing_mode="<<(requested<=0 ? "all_cross_event_pairs" : "random_pairs")
       <<"\nrequested="<<requested<<"\nseed="<<(requested<=0 ? 0 : seed)
       <<"\npool1_candidates="<<pool1.size()<<"\npool2_candidates="<<pool2.size()
       <<"\ncartesian_pairs="<<Long64_t(pool1.size())*Long64_t(pool2.size())
       <<"\nmixing=Jpsi1_from_source_A_times_Jpsi2_from_source_B"
       <<"\nsame_event_identity=source_entry_with_run_lumi_event_fallback"
       <<"\nattempts="<<attempts<<"\nsame_event_rejected="<<sameRejected
       <<"\nmass_rejected="<<massRejected<<"\ncorrection_rejected="<<correctionRejected
       <<"\npositive_pairs="<<positive<<"\nnegative_pairs="<<negative<<"\nzero_pairs="<<zero
       <<"\naccepted_pairs="<<tree.GetEntries()
       <<"\nsumw="<<sumw<<"\nsumw2="<<sumw2<<"\nmax_abs_weight="<<maxAbs
       <<"\npair_weight=prompt_sweight_Jpsi1_times_prompt_sweight_Jpsi2_times_nominal_correction\n";
  TNamed metadata("run_metadata",meta.str().c_str()); tree.Write();metadata.Write();fo.Close();
  std::cout << meta.str()<<"output="<<output<<"\n";
}
