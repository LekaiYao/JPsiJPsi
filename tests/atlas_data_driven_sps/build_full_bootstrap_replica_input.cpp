#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TLeaf.h"
#include "TNamed.h"
#include "TSystem.h"
#include "TTree.h"

namespace {
bool closeEnough(double a,double b) {
  return std::abs(a-b)<=1e-9*std::max(1.0,std::max(std::abs(a),std::abs(b)));
}

bool aligned(TTree &jpsi1,TTree &jpsi2,TTree &data,Long64_t entry) {
  jpsi1.GetEntry(entry); jpsi2.GetEntry(entry); data.GetEntry(entry);
  const char *slot1[]={"Jpsi_mass1","Jpsi_ctau1","Jpsi_pt1","Jpsi_y1"};
  const char *slot2[]={"Jpsi_mass2","Jpsi_ctau2","Jpsi_pt2","Jpsi_y2"};
  for (const char *name:slot1)
    if (!jpsi1.GetLeaf(name)||!data.GetLeaf(name)||
        !closeEnough(jpsi1.GetLeaf(name)->GetValue(),data.GetLeaf(name)->GetValue())) return false;
  for (const char *name:slot2)
    if (!jpsi2.GetLeaf(name)||!data.GetLeaf(name)||
        !closeEnough(jpsi2.GetLeaf(name)->GetValue(),data.GetLeaf(name)->GetValue())) return false;
  return true;
}
}

void build_full_bootstrap_replica_input(
    int replicaId=0,ULong64_t baseSeed=20260806,
    const char *outputDir=
      "tests/atlas_data_driven_sps/results/route9p3_crossslot_category_absy1p2_allpairs/bootstrap_full/replica_000000/input",
    const char *candidateInput=
      "tests/atlas_data_driven_sps/results/route9p3_crossslot_splot/jpsi12_candidates.root",
    const char *dataInput=
      "tests/atlas_data_driven_sps/results/unified_input/WeightData_mJJ7p5.root") {
  TFile fc(candidateInput,"READ"),fd(dataInput,"READ");
  TTree *jpsi1=dynamic_cast<TTree *>(fc.Get("jpsi1"));
  TTree *jpsi2=dynamic_cast<TTree *>(fc.Get("jpsi2"));
  TTree *data=dynamic_cast<TTree *>(fd.Get("data"));
  if (!jpsi1||!jpsi2||!data||jpsi1->GetEntries()!=jpsi2->GetEntries()||
      jpsi1->GetEntries()!=data->GetEntries()||jpsi1->GetEntries()<=0) {
    Error("build_full_bootstrap_replica_input","missing or unequal input trees");
    gSystem->Exit(2); return;
  }
  const Long64_t n=jpsi1->GetEntries();
  for (Long64_t i=0;i<n;++i) if (!aligned(*jpsi1,*jpsi2,*data,i)) {
    Error("build_full_bootstrap_replica_input","candidate/Data mismatch at entry %lld",i);
    gSystem->Exit(3); return;
  }

  const ULong64_t seed=baseSeed+1000003ULL*static_cast<ULong64_t>(replicaId);
  std::mt19937_64 rng(seed); std::uniform_int_distribution<Long64_t> pick(0,n-1);
  std::vector<Long64_t> draws(n); std::map<Long64_t,int> multiplicity;
  for (Long64_t i=0;i<n;++i) { draws[i]=pick(rng); ++multiplicity[draws[i]]; }

  gSystem->mkdir(outputDir,true); const std::string base(outputDir);
  TFile candidateOutput((base+"/jpsi12_candidates.root").c_str(),"RECREATE");
  TTree *out1=jpsi1->CloneTree(0),*out2=jpsi2->CloneTree(0);
  candidateOutput.cd();
  for (Long64_t index:draws) { jpsi1->GetEntry(index); out1->Fill(); jpsi2->GetEntry(index); out2->Fill(); }
  out1->Write();out2->Write();candidateOutput.Close();

  TFile dataOutput((base+"/WeightData_mJJ7p5.root").c_str(),"RECREATE");
  TTree *outData=data->CloneTree(0);
  for (Long64_t index:draws) { data->GetEntry(index); outData->Fill(); }
  outData->Write();dataOutput.Close();

  long long sum=0,sum2=0; int maximum=0;
  std::ofstream csv(base+"/source_event_multiplicities.csv");
  csv<<"source_entry,multiplicity\n";
  for (const auto &item:multiplicity) {
    csv<<item.first<<","<<item.second<<"\n"; sum+=item.second;
    sum2+=static_cast<long long>(item.second)*item.second; maximum=std::max(maximum,item.second);
  }
  std::ostringstream text;
  text<<std::setprecision(12)<<"scope=full_source_event_bootstrap_primary_2d_fdps\n"
      <<"replica_id="<<replicaId<<"\nbase_seed="<<baseSeed<<"\nseed="<<seed
      <<"\nsource_events="<<n<<"\ndraws="<<draws.size()
      <<"\nunique_source_events="<<multiplicity.size()<<"\nsum_multiplicity="<<sum
      <<"\nsum_multiplicity2="<<sum2<<"\nmax_multiplicity="<<maximum
      <<"\ncandidate_data_entry_alignment_verified=true\n"
      <<"bootstrap_unit=selected_source_event_pair\nslot_pairing=Jpsi1_and_Jpsi2_resampled_together\n";
  std::ofstream meta(base+"/metadata.txt");meta<<text.str();meta.close();csv.close();
  if (!meta||!csv||sum!=n) { Error("build_full_bootstrap_replica_input","output invariant failed"); gSystem->Exit(4); return; }
  Printf("full bootstrap input replica %d complete: draws=%lld unique=%zu",replicaId,n,multiplicity.size());
}
