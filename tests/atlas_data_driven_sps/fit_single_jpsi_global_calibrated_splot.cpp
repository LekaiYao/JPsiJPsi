#include "fit_single_jpsi_category_splot.cpp"

void fit_single_jpsi_global_calibrated_splot(
    int slot=1,
    const char *input=
        "tests/atlas_data_driven_sps/results/route9p3_crossslot_splot/jpsi12_candidates.root",
    const char *modelFile=
        "Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/Model_4D_tot.root",
    const char *output=
        "tests/atlas_data_driven_sps/results/fdps_systematics_final_v1/category/global_calibrated/jpsi1_sweights.root") {
  if (slot!=1 && slot!=2) { std::cerr<<"slot must be 1 or 2\n"; return; }
  const std::string suffix=std::to_string(slot), treeName="jpsi"+suffix;
  const std::string massName="Jpsi_mass"+suffix, ctauName="Jpsi_ctau"+suffix;
  TFile fi(input,"READ"); TTree *tree=dynamic_cast<TTree *>(fi.Get(treeName.c_str()));
  TFile fm(modelFile,"READ"); RooWorkspace *wsp=dynamic_cast<RooWorkspace *>(fm.Get("wsp"));
  if (!tree||!wsp) { std::cerr<<"Missing candidate tree or workspace\n"; return; }
  RooRealVar *mass=wsp->var(massName.c_str()), *ctau=wsp->var(ctauName.c_str());
  RooAbsPdf *massSig=wsp->pdf(("JpsiMassSig"+suffix).c_str());
  RooAbsPdf *massComb=wsp->pdf(("JpsiMassComb"+suffix).c_str());
  RooAbsPdf *ctauPrompt=wsp->pdf(("JpsiCtauSig"+suffix).c_str());
  RooAbsPdf *ctauNonprompt=wsp->pdf(("JpsiCtauBkg"+suffix).c_str());
  RooAbsPdf *ctauComb=wsp->pdf(("JpsiCtauCombBkg"+suffix).c_str());
  if (!mass||!ctau||!massSig||!massComb||!ctauPrompt||!ctauNonprompt||!ctauComb||
      !wsp->var("Jpsi_devia1")||!wsp->var("Jpsi_devia2")||
      !wsp->var("Jpsi_sigma1")||!wsp->var("Jpsi_sigma2")) {
    std::cerr<<"Missing required global-fit object\n"; return;
  }
  mass->setRange(2.95,3.25); ctau->setRange(-0.03,0.16);
  RooDataSet data("data","unweighted selected global candidates",tree,RooArgSet(*mass,*ctau));
  if (data.numEntries()<100) { std::cerr<<"Global fit has too few entries\n"; return; }

  CategoryModel model=makeCategoryModel("s"+suffix+"_global_calibrated",data.numEntries(),
      *wsp,*massSig,*ctauPrompt,*massComb,*ctauNonprompt,*ctauComb);
  CalibrationSeed seed;
  int calibrationAttempts=0;
  std::string calibrationPath;
  RooFitResult *calibration=calibrateCategory(
      model,data,seed,calibrationAttempts,calibrationPath);
  if (!goodFit(calibration)) {
    std::cerr<<"Global calibration failure slot="<<slot<<"\n";
    delete calibration; return;
  }
  const double massScale=model.massScale->getVal();
  const double massScaleError=model.massScale->getError();
  const double ctauScale=model.ctauScale->getVal();
  const double ctauScaleError=model.ctauScale->getError();
  model.massScale->setConstant(true); model.ctauScale->setConstant(true);
  RooFitResult *yieldFit=model.model->fitTo(data,Extended(true),Save(true),Offset(true),
      Optimize(false),Strategy(1),Minimizer("Minuit2","migrad"),PrintLevel(-1));
  int yieldAttempts=1;
  if (!goodFit(yieldFit)) {
    delete yieldFit; ++yieldAttempts;
    yieldFit=model.model->fitTo(data,Extended(true),Save(true),Offset(true),
        Optimize(false),Strategy(2),Minimizer("Minuit2","migrad"),PrintLevel(-1));
  }
  if (!goodFit(yieldFit)) {
    std::cerr<<"Global yield-only fit failure slot="<<slot<<"\n";
    delete calibration; delete yieldFit; return;
  }
  RooStats::SPlot splot(("splot_s"+suffix+"_global_calibrated").c_str(),
      "global calibrated sPlot",data,model.model,
      RooArgList(*model.nPrompt,*model.nNonprompt,*model.nComb));

  TString outDir=gSystem->DirName(output); gSystem->mkdir(outDir,true);
  TFile fo(output,"RECREATE"); fo.cd();
  std::ostringstream fitRange;
  fitRange<<massName<<">=2.95&&"<<massName<<"<=3.25&&"
          <<ctauName<<">=-0.03&&"<<ctauName<<"<=0.16";
  TTree *out=tree->CopyTree(fitRange.str().c_str());
  out->SetName("candidates");
  out->SetTitle(("Jpsi"+suffix+" with global-calibrated 2D sWeights").c_str());
  double promptWeight=0,nonpromptWeight=0,combWeight=0;
  TBranch *bp=out->Branch("prompt_sweight",&promptWeight);
  TBranch *bn=out->Branch("nonprompt_sweight",&nonpromptWeight);
  TBranch *bc=out->Branch("comb_sweight",&combWeight);
  Long64_t negativePrompt=0;
  double sumPrompt=0,sumPrompt2=0;
  for (Long64_t i=0;i<out->GetEntries();++i) {
    out->GetEntry(i);
    promptWeight=splot.GetSWeight(i,model.nPrompt->GetName());
    nonpromptWeight=splot.GetSWeight(i,model.nNonprompt->GetName());
    combWeight=splot.GetSWeight(i,model.nComb->GetName());
    if (promptWeight<0) ++negativePrompt;
    sumPrompt+=promptWeight; sumPrompt2+=promptWeight*promptWeight;
    bp->Fill(); bn->Fill(); bc->Fill();
  }
  if (out->GetEntries()!=data.numEntries()) {
    std::cerr<<"Tree/global dataset alignment failure\n"; fo.Close(); return;
  }
  calibration->Write("global_calibration_fit_result");
  yieldFit->Write("global_yield_fit_result");
  out->Write("",TObject::kOverwrite);
  std::ostringstream metadata;
  metadata<<std::setprecision(12)<<"slot="<<slot
      <<"\nfit_mode=single_global_calibrated_two_stage"
      <<"\nentries="<<data.numEntries()
      <<"\ncalibration_optimizer=central_and_legacy_yield_multistart_keep_best_with_dual_fallback"
      <<"\ncalibration_status="<<calibration->status()
      <<"\ncalibration_covQual="<<calibration->covQual()
      <<"\ncalibration_edm="<<calibration->edm()
      <<"\ncalibration_attempts="<<calibrationAttempts
      <<"\ncalibration_best_path="<<calibrationPath
      <<"\nmass_scale="<<massScale<<"\nmass_scale_error="<<massScaleError
      <<"\nctau_scale="<<ctauScale<<"\nctau_scale_error="<<ctauScaleError
      <<"\nyield_fit_status="<<yieldFit->status()
      <<"\nyield_fit_covQual="<<yieldFit->covQual()
      <<"\nyield_fit_edm="<<yieldFit->edm()
      <<"\nyield_fit_attempts="<<yieldAttempts
      <<"\nprompt_yield="<<model.nPrompt->getVal()
      <<"\nnonprompt_yield="<<model.nNonprompt->getVal()
      <<"\ncomb_yield="<<model.nComb->getVal()
      <<"\nsum_prompt_sweight="<<sumPrompt
      <<"\nsum_prompt_sweight2="<<sumPrompt2
      <<"\nnegative_prompt_sweights="<<negativePrompt
      <<"\nfit_input_weight=none"
      <<"\nresolution_scales_frozen_before_splot=true\n";
  TNamed runMetadata("run_metadata",metadata.str().c_str()); runMetadata.Write(); fo.Close();
  std::ofstream summary(std::string(outDir.Data())+"/global_calibrated_splot_summary_jpsi"+suffix+".txt");
  summary<<metadata.str();
  std::cout<<metadata.str()<<"output="<<output<<"\n";
  delete calibration; delete yieldFit;
}
