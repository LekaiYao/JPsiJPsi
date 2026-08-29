#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooCategory.h"
#include "RooConstVar.h"
#include "RooCustomizer.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "RooStats/SPlot.h"
#include "RooWorkspace.h"
#include "TFile.h"
#include "TNamed.h"
#include "TSystem.h"
#include "TTree.h"

using namespace RooFit;

namespace {
struct CategoryModel {
  RooRealVar *massScale=0, *ctauScale=0;
  RooConstVar *baseDevia1=0, *baseDevia2=0, *baseSigma1=0, *baseSigma2=0;
  RooFormulaVar *devia1=0, *devia2=0, *sigma1=0, *sigma2=0;
  RooAbsPdf *massSig=0, *ctauPrompt=0;
  RooProdPdf *prompt=0, *nonprompt=0, *comb=0;
  RooRealVar *nPrompt=0, *nNonprompt=0, *nComb=0;
  RooAddPdf *model=0;
};

struct CalibrationSeed {
  double massScale=1.0, ctauScale=1.0;
  double promptFraction=0.75, nonpromptFraction=0.20, combFraction=0.05;
};

bool goodFit(const RooFitResult *fit) {
  return fit && fit->status()==0 && fit->covQual()>=2 && fit->edm()<0.01;
}

void resetFitSteps(CategoryModel &c) {
  c.massScale->setError(0.05); c.ctauScale->setError(0.05);
  c.nPrompt->setError(std::sqrt(std::max(1.0,c.nPrompt->getVal())));
  c.nNonprompt->setError(std::sqrt(std::max(1.0,c.nNonprompt->getVal())));
  c.nComb->setError(std::sqrt(std::max(1.0,c.nComb->getVal())));
}

void resetFitStart(CategoryModel &c, const CalibrationSeed &seed) {
  const double entries=c.nPrompt->getMax()/1.5;
  const double fractionSum=seed.promptFraction+seed.nonpromptFraction+seed.combFraction;
  c.massScale->setVal(seed.massScale); c.ctauScale->setVal(seed.ctauScale);
  c.nPrompt->setVal(entries*seed.promptFraction/fractionSum);
  c.nNonprompt->setVal(entries*seed.nonpromptFraction/fractionSum);
  c.nComb->setVal(entries*seed.combFraction/fractionSum); resetFitSteps(c);
}

RooFitResult *calibrateCategory(CategoryModel &model, RooDataSet &data,
                                const CalibrationSeed &seed, int &attempts,
                                std::string &bestPath) {
  model.massScale->setConstant(false); model.ctauScale->setConstant(false);
  const CalibrationSeed legacySeed;
  const double centralStarts[4]={seed.massScale,0.75,1.0,1.5};
  const double legacyStarts[3]={0.75,1.0,1.5};
  bool found=false; double bestNll=0,bestMass=1,bestCtau=1;
  double bestPrompt=0,bestNonprompt=0,bestComb=0;
  double bestMassError=0.05,bestCtauError=0.05;
  RooFitResult *bestFit=0;
  attempts=0;
  bestPath="none";
  const auto considerTrial = [&](RooFitResult *trial, const char *path) {
    if (goodFit(trial) && (!found||trial->minNll()<bestNll)) {
      delete bestFit; bestFit=trial;
      found=true; bestNll=bestFit->minNll(); bestPath=path;
      bestMass=model.massScale->getVal(); bestCtau=model.ctauScale->getVal();
      bestMassError=model.massScale->getError(); bestCtauError=model.ctauScale->getError();
      bestPrompt=model.nPrompt->getVal(); bestNonprompt=model.nNonprompt->getVal();
      bestComb=model.nComb->getVal();
    } else {
      delete trial;
    }
  };
  for (int i=0;i<4;++i) {
    for (int strategy=1;strategy<=2;++strategy) {
      resetFitStart(model,seed); model.massScale->setVal(centralStarts[i]);
      RooFitResult *trial=model.model->fitTo(data,Extended(true),Save(true),Offset(true),
          Optimize(false),Strategy(strategy),Minimizer("Minuit2","migrad"),PrintLevel(-1));
      ++attempts;
      considerTrial(trial,"central_yields_primary");
    }
  }
  for (int i=0;i<3;++i) {
    for (int strategy=1;strategy<=2;++strategy) {
      resetFitStart(model,legacySeed); model.massScale->setVal(legacyStarts[i]);
      RooFitResult *trial=model.model->fitTo(data,Extended(true),Save(true),Offset(true),
          Optimize(false),Strategy(strategy),Minimizer("Minuit2","migrad"),PrintLevel(-1));
      ++attempts;
      considerTrial(trial,"legacy_yields_primary");
    }
  }
  if (!found) {
    resetFitStart(model,seed);
    model.massScale->setConstant(true); model.ctauScale->setConstant(true);
    RooFitResult *yieldPrefit=model.model->fitTo(data,Extended(true),Save(true),Offset(true),
        Optimize(false),Strategy(2),Minimizer("Minuit2","migrad"),PrintLevel(-1));
    ++attempts;
    const double warmPrompt=model.nPrompt->getVal();
    const double warmNonprompt=model.nNonprompt->getVal();
    const double warmComb=model.nComb->getVal();
    delete yieldPrefit;
    model.massScale->setConstant(false); model.ctauScale->setConstant(false);
    model.massScale->setError(0.05); model.ctauScale->setError(0.05);
    RooFitResult *warmTrial=model.model->fitTo(data,Extended(true),Save(true),Offset(true),
        Optimize(false),Strategy(2),Minimizer("Minuit2","migrad"),PrintLevel(-1));
    ++attempts;
    considerTrial(warmTrial,"fixed_scale_yield_prefit_release");
    const double fallbackMass[]={0.65,0.8,1.0,1.25,1.5,1.8};
    const double fallbackCtau[]={0.65,0.85,1.05,1.3,1.6};
    for (double massStart:fallbackMass) {
      for (double ctauStart:fallbackCtau) {
        resetFitStart(model,seed); model.massScale->setVal(massStart);
        model.ctauScale->setVal(ctauStart);
        model.nPrompt->setVal(warmPrompt); model.nNonprompt->setVal(warmNonprompt);
        model.nComb->setVal(warmComb); resetFitSteps(model);
        RooFitResult *trial=model.model->fitTo(data,Extended(true),Save(true),Offset(true),
            Optimize(false),Strategy(2),Minimizer("Minuit2","migrad"),PrintLevel(-1));
        ++attempts;
        considerTrial(trial,"prefit_yields_fallback_grid");
      }
    }
    for (double massStart:fallbackMass) {
      for (double ctauStart:fallbackCtau) {
        resetFitStart(model,legacySeed); model.massScale->setVal(massStart);
        model.ctauScale->setVal(ctauStart);
        RooFitResult *trial=model.model->fitTo(data,Extended(true),Save(true),Offset(true),
            Optimize(false),Strategy(2),Minimizer("Minuit2","migrad"),PrintLevel(-1));
        ++attempts;
        considerTrial(trial,"legacy_yields_fallback_grid");
      }
    }
  }
  if (!found) return 0;
  model.massScale->setVal(bestMass); model.ctauScale->setVal(bestCtau);
  model.nPrompt->setVal(bestPrompt); model.nNonprompt->setVal(bestNonprompt);
  model.nComb->setVal(bestComb); resetFitSteps(model);
  model.massScale->setError(bestMassError); model.ctauScale->setError(bestCtauError);
  return bestFit;
}

CalibrationSeed categorySeed(int slot, const std::string &variable, bool high) {
  CalibrationSeed seed;
  if (variable!="absy") return seed;
  if (slot==1 && !high) {
    seed.massScale=0.85634113833; seed.ctauScale=0.994608187154;
    seed.promptFraction=0.265440549116; seed.nonpromptFraction=0.692003329634;
    seed.combFraction=0.042557769145;
  } else if (slot==1) {
    seed.massScale=1.47874770968; seed.ctauScale=1.07446139402;
    seed.promptFraction=0.359792451461; seed.nonpromptFraction=0.617894187779;
    seed.combFraction=0.022312816692;
  } else if (!high) {
    seed.massScale=0.827714056645; seed.ctauScale=0.894987798789;
    seed.promptFraction=0.278188663791; seed.nonpromptFraction=0.665620267909;
    seed.combFraction=0.056189924287;
  } else {
    seed.massScale=1.41151480695; seed.ctauScale=1.08890105628;
    seed.promptFraction=0.328460576157; seed.nonpromptFraction=0.651795244386;
    seed.combFraction=0.019744253633;
  }
  return seed;
}

CategoryModel makeCategoryModel(
    const std::string &tag, double entries, RooWorkspace &wsp,
    RooAbsPdf &baseMassSig, RooAbsPdf &baseCtauPrompt,
    RooAbsPdf &massComb, RooAbsPdf &ctauNonprompt, RooAbsPdf &ctauComb) {
  CategoryModel c;
  c.massScale=new RooRealVar(("mass_scale_"+tag).c_str(),"mass resolution scale",1.0,0.5,2.0);
  c.ctauScale=new RooRealVar(("ctau_scale_"+tag).c_str(),"prompt ctau resolution scale",1.0,0.5,2.0);
  RooRealVar *d1=wsp.var("Jpsi_devia1"), *d2=wsp.var("Jpsi_devia2");
  RooRealVar *s1=wsp.var("Jpsi_sigma1"), *s2=wsp.var("Jpsi_sigma2");
  c.baseDevia1=new RooConstVar(("base_devia1_"+tag).c_str(),"",d1->getVal());
  c.baseDevia2=new RooConstVar(("base_devia2_"+tag).c_str(),"",d2->getVal());
  c.baseSigma1=new RooConstVar(("base_sigma1_"+tag).c_str(),"",s1->getVal());
  c.baseSigma2=new RooConstVar(("base_sigma2_"+tag).c_str(),"",s2->getVal());
  c.devia1=new RooFormulaVar(("devia1_"+tag).c_str(),"@0*@1",RooArgList(*c.baseDevia1,*c.massScale));
  c.devia2=new RooFormulaVar(("devia2_"+tag).c_str(),"@0*@1",RooArgList(*c.baseDevia2,*c.massScale));
  c.sigma1=new RooFormulaVar(("sigma1_"+tag).c_str(),"@0*@1",RooArgList(*c.baseSigma1,*c.ctauScale));
  c.sigma2=new RooFormulaVar(("sigma2_"+tag).c_str(),"@0*@1",RooArgList(*c.baseSigma2,*c.ctauScale));
  RooCustomizer massCustomizer(baseMassSig,("mass_"+tag).c_str());
  massCustomizer.replaceArg(*d1,*c.devia1); massCustomizer.replaceArg(*d2,*c.devia2);
  c.massSig=dynamic_cast<RooAbsPdf *>(massCustomizer.build());
  RooCustomizer ctauCustomizer(baseCtauPrompt,("ctau_"+tag).c_str());
  ctauCustomizer.replaceArg(*s1,*c.sigma1); ctauCustomizer.replaceArg(*s2,*c.sigma2);
  c.ctauPrompt=dynamic_cast<RooAbsPdf *>(ctauCustomizer.build());
  c.prompt=new RooProdPdf(("prompt_"+tag).c_str(),"prompt",RooArgSet(*c.massSig,*c.ctauPrompt));
  c.nonprompt=new RooProdPdf(("nonprompt_"+tag).c_str(),"nonprompt",RooArgSet(*c.massSig,ctauNonprompt));
  c.comb=new RooProdPdf(("comb_"+tag).c_str(),"combinatorial",RooArgSet(massComb,ctauComb));
  c.nPrompt=new RooRealVar(("n_prompt_"+tag).c_str(),"prompt yield",0.75*entries,0,1.5*entries);
  c.nNonprompt=new RooRealVar(("n_nonprompt_"+tag).c_str(),"nonprompt yield",0.20*entries,0,1.5*entries);
  c.nComb=new RooRealVar(("n_comb_"+tag).c_str(),"combinatorial yield",0.05*entries,0,1.5*entries);
  c.model=new RooAddPdf(("model_"+tag).c_str(),"three-component category model",
                        RooArgList(*c.prompt,*c.nonprompt,*c.comb),
                        RooArgList(*c.nPrompt,*c.nNonprompt,*c.nComb));
  return c;
}
}

void fit_single_jpsi_category_splot(
    int slot=1, const char *categoryVariable="pt", double split=15.0,
    const char *input=
        "tests/atlas_data_driven_sps/results/route9p3_crossslot_splot/jpsi12_candidates.root",
    const char *modelFile=
        "Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/Model_4D_tot.root",
    const char *output=
        "tests/atlas_data_driven_sps/results/route9p3_crossslot_category_pt15_allpairs/jpsi1_sweights.root") {
  if (slot!=1 && slot!=2) { std::cerr<<"slot must be 1 or 2\n"; return; }
  const std::string variable(categoryVariable);
  if (variable!="pt" && variable!="absy") { std::cerr<<"categoryVariable must be pt or absy\n"; return; }
  const std::string suffix=std::to_string(slot), treeName="jpsi"+suffix;
  const std::string massName="Jpsi_mass"+suffix, ctauName="Jpsi_ctau"+suffix;
  const std::string ptName="Jpsi_pt"+suffix, yName="Jpsi_y"+suffix;
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
    std::cerr<<"Missing required category-fit object\n"; return;
  }
  mass->setRange(2.95,3.25); ctau->setRange(-0.03,0.16);
  RooRealVar pt(ptName.c_str(),ptName.c_str(),10,40);
  RooRealVar rapidity(yName.c_str(),yName.c_str(),-3,3);
  RooArgSet allObservables(*mass,*ctau,pt,rapidity);
  RooDataSet all("all","unweighted selected candidates",tree,allObservables);
  std::ostringstream lowCut,highCut;
  if (variable=="pt") {
    lowCut<<ptName<<"<"<<std::setprecision(12)<<split;
    highCut<<ptName<<">="<<std::setprecision(12)<<split;
  } else {
    lowCut<<"abs("<<yName<<")<"<<std::setprecision(12)<<split;
    highCut<<"abs("<<yName<<")>="<<std::setprecision(12)<<split;
  }
  RooDataSet *low=dynamic_cast<RooDataSet *>(all.reduce(RooArgSet(*mass,*ctau),lowCut.str().c_str()));
  RooDataSet *high=dynamic_cast<RooDataSet *>(all.reduce(RooArgSet(*mass,*ctau),highCut.str().c_str()));
  if (!low||!high||low->numEntries()<100||high->numEntries()<100) {
    std::cerr<<"Category has too few entries\n"; return;
  }
  const std::string baseTag="s"+suffix+"_"+variable;
  CategoryModel lowModel=makeCategoryModel(baseTag+"_low",low->numEntries(),*wsp,
      *massSig,*ctauPrompt,*massComb,*ctauNonprompt,*ctauComb);
  CategoryModel highModel=makeCategoryModel(baseTag+"_high",high->numEntries(),*wsp,
      *massSig,*ctauPrompt,*massComb,*ctauNonprompt,*ctauComb);
  int lowCalibrationAttempts=0,highCalibrationAttempts=0;
  std::string lowCalibrationPath,highCalibrationPath;
  RooFitResult *lowCalibration=calibrateCategory(
      lowModel,*low,categorySeed(slot,variable,false),lowCalibrationAttempts,lowCalibrationPath);
  RooFitResult *highCalibration=calibrateCategory(
      highModel,*high,categorySeed(slot,variable,true),highCalibrationAttempts,highCalibrationPath);
  if (!goodFit(lowCalibration)||!goodFit(highCalibration)) {
    std::cerr<<"Independent category calibration failure slot="<<slot<<"\n";
    delete lowCalibration; delete highCalibration; return;
  }
  const int calibrationStatus=std::max(lowCalibration->status(),highCalibration->status());
  const int calibrationCovQual=std::min(lowCalibration->covQual(),highCalibration->covQual());
  const double calibrationEdm=std::max(lowCalibration->edm(),highCalibration->edm());
  const int calibrationAttempts=lowCalibrationAttempts+highCalibrationAttempts;
  const double lowMassScale=lowModel.massScale->getVal(), lowMassScaleError=lowModel.massScale->getError();
  const double lowCtauScale=lowModel.ctauScale->getVal(), lowCtauScaleError=lowModel.ctauScale->getError();
  const double highMassScale=highModel.massScale->getVal(), highMassScaleError=highModel.massScale->getError();
  const double highCtauScale=highModel.ctauScale->getVal(), highCtauScaleError=highModel.ctauScale->getError();
  lowModel.massScale->setConstant(true); lowModel.ctauScale->setConstant(true);
  highModel.massScale->setConstant(true); highModel.ctauScale->setConstant(true);
  RooFitResult *lowFit=lowModel.model->fitTo(*low,Extended(true),Save(true),Offset(true),
      Optimize(false),Strategy(1),Minimizer("Minuit2","migrad"),PrintLevel(-1));
  RooFitResult *highFit=highModel.model->fitTo(*high,Extended(true),Save(true),Offset(true),
      Optimize(false),Strategy(1),Minimizer("Minuit2","migrad"),PrintLevel(-1));
  int lowFitAttempts=1,highFitAttempts=1;
  if (!goodFit(lowFit)) {
    delete lowFit; ++lowFitAttempts;
    lowFit=lowModel.model->fitTo(*low,Extended(true),Save(true),Offset(true),
        Optimize(false),Strategy(2),Minimizer("Minuit2","migrad"),PrintLevel(-1));
  }
  if (!goodFit(highFit)) {
    delete highFit; ++highFitAttempts;
    highFit=highModel.model->fitTo(*high,Extended(true),Save(true),Offset(true),
        Optimize(false),Strategy(2),Minimizer("Minuit2","migrad"),PrintLevel(-1));
  }
  if (!goodFit(lowFit)||!goodFit(highFit)) {
    std::cerr<<"Yield-only category fit failure slot="<<slot<<"\n";
    delete lowCalibration; delete highCalibration; delete lowFit; delete highFit; return;
  }
  RooStats::SPlot lowSPlot(("splot_"+baseTag+"_low").c_str(),"low category sPlot",
      *low,lowModel.model,RooArgList(*lowModel.nPrompt,*lowModel.nNonprompt,*lowModel.nComb));
  RooStats::SPlot highSPlot(("splot_"+baseTag+"_high").c_str(),"high category sPlot",
      *high,highModel.model,RooArgList(*highModel.nPrompt,*highModel.nNonprompt,*highModel.nComb));

  TString outDir=gSystem->DirName(output); gSystem->mkdir(outDir,true);
  TFile fo(output,"RECREATE"); fo.cd();
  std::ostringstream fitRange;
  fitRange<<massName<<">=2.95&&"<<massName<<"<=3.25&&"
          <<ctauName<<">=-0.03&&"<<ctauName<<"<=0.16";
  TTree *out=tree->CopyTree(fitRange.str().c_str());
  out->SetName("candidates");
  out->SetTitle(("Jpsi"+suffix+" with category-calibrated 2D sWeights").c_str());
  double treePt=0,treeY=0,promptWeight=0,nonpromptWeight=0,combWeight=0,categoryValue=0;
  int categoryIndex=0;
  out->SetBranchAddress(ptName.c_str(),&treePt); out->SetBranchAddress(yName.c_str(),&treeY);
  TBranch *bp=out->Branch("prompt_sweight",&promptWeight);
  TBranch *bn=out->Branch("nonprompt_sweight",&nonpromptWeight);
  TBranch *bc=out->Branch("comb_sweight",&combWeight);
  TBranch *bi=out->Branch("category_index",&categoryIndex);
  TBranch *bv=out->Branch("category_value",&categoryValue);
  Long64_t index[2]={0,0},negative[2]={0,0}; double sumPrompt[2]={0,0},sumPrompt2[2]={0,0};
  for (Long64_t i=0;i<out->GetEntries();++i) {
    out->GetEntry(i); categoryValue=(variable=="pt" ? treePt : std::abs(treeY));
    categoryIndex=(categoryValue<split ? 0 : 1);
    RooStats::SPlot &sp=(categoryIndex==0 ? lowSPlot : highSPlot);
    CategoryModel &cm=(categoryIndex==0 ? lowModel : highModel);
    promptWeight=sp.GetSWeight(index[categoryIndex],cm.nPrompt->GetName());
    nonpromptWeight=sp.GetSWeight(index[categoryIndex],cm.nNonprompt->GetName());
    combWeight=sp.GetSWeight(index[categoryIndex],cm.nComb->GetName());
    if (promptWeight<0) ++negative[categoryIndex];
    sumPrompt[categoryIndex]+=promptWeight; sumPrompt2[categoryIndex]+=promptWeight*promptWeight;
    ++index[categoryIndex]; bp->Fill(); bn->Fill(); bc->Fill(); bi->Fill(); bv->Fill();
  }
  if (index[0]!=low->numEntries()||index[1]!=high->numEntries()) {
    std::cerr<<"Tree/category dataset alignment failure\n"; fo.Close(); return;
  }
  lowCalibration->Write("low_calibration_fit_result");
  highCalibration->Write("high_calibration_fit_result"); lowFit->Write("low_yield_fit_result");
  highFit->Write("high_yield_fit_result"); out->Write("",TObject::kOverwrite);
  std::ostringstream metadata;
  metadata<<std::setprecision(12)<<"slot="<<slot<<"\ncategory_variable="<<variable
      <<"\ncategory_split="<<split<<"\nlow_entries="<<low->numEntries()
      <<"\nhigh_entries="<<high->numEntries()<<"\ncalibration_mode=independent_factorized_categories"
      <<"\ncalibration_optimizer=central_and_legacy_yield_multistart_keep_best_with_dual_fallback"
      <<"\ncalibration_status="<<calibrationStatus
      <<"\ncalibration_attempts="<<calibrationAttempts
      <<"\ncalibration_covQual="<<calibrationCovQual<<"\ncalibration_edm="<<calibrationEdm
      <<"\nlow_calibration_attempts="<<lowCalibrationAttempts
      <<"\nhigh_calibration_attempts="<<highCalibrationAttempts
      <<"\nlow_calibration_best_path="<<lowCalibrationPath
      <<"\nhigh_calibration_best_path="<<highCalibrationPath
      <<"\nlow_mass_scale="<<lowMassScale<<"\nlow_mass_scale_error="<<lowMassScaleError
      <<"\nlow_ctau_scale="<<lowCtauScale<<"\nlow_ctau_scale_error="<<lowCtauScaleError
      <<"\nhigh_mass_scale="<<highMassScale<<"\nhigh_mass_scale_error="<<highMassScaleError
      <<"\nhigh_ctau_scale="<<highCtauScale<<"\nhigh_ctau_scale_error="<<highCtauScaleError
      <<"\nlow_fit_status="<<lowFit->status()<<"\nlow_fit_attempts="<<lowFitAttempts
      <<"\nlow_fit_covQual="<<lowFit->covQual()
      <<"\nlow_fit_edm="<<lowFit->edm()<<"\nhigh_fit_status="<<highFit->status()
      <<"\nhigh_fit_attempts="<<highFitAttempts
      <<"\nhigh_fit_covQual="<<highFit->covQual()<<"\nhigh_fit_edm="<<highFit->edm()
      <<"\nlow_prompt_yield="<<lowModel.nPrompt->getVal()
      <<"\nlow_sum_prompt_sweight="<<sumPrompt[0]
      <<"\nlow_sum_prompt_sweight2="<<sumPrompt2[0]
      <<"\nlow_negative_prompt_sweights="<<negative[0]
      <<"\nhigh_prompt_yield="<<highModel.nPrompt->getVal()
      <<"\nhigh_sum_prompt_sweight="<<sumPrompt[1]
      <<"\nhigh_sum_prompt_sweight2="<<sumPrompt2[1]
      <<"\nhigh_negative_prompt_sweights="<<negative[1]
      <<"\nfit_input_weight=none\nresolution_scales_frozen_before_splot=true\n";
  TNamed runMetadata("run_metadata",metadata.str().c_str()); runMetadata.Write(); fo.Close();
  std::ofstream summary(std::string(outDir.Data())+"/category_splot_summary_jpsi"+suffix+".txt");
  summary<<metadata.str();
  std::cout<<metadata.str()<<"output="<<output<<"\n";
  delete lowCalibration; delete highCalibration; delete lowFit; delete highFit;
}
