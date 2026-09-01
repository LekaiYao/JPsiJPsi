#include <iostream>
#include <memory>
#include <string>

#include "RooAbsArg.h"
#include "RooAbsCollection.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TIterator.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPad.h"
#include "RVersion.h"
#include "TStyle.h"

// Reuse the exact accepted model constructor.  This include only provides
// declarations and implementations; this plotting entry point never calls
// fitTo and never changes the saved fit result.
#include "../../fit_single_jpsi_category_splot.cpp"

namespace {

void assignValues(RooAbsCollection &destination, const RooAbsCollection &source) {
  const auto assignOne = [&](RooAbsArg *argument) {
    RooRealVar *sourceVariable=dynamic_cast<RooRealVar *>(argument);
    RooRealVar *destinationVariable=
        dynamic_cast<RooRealVar *>(destination.find(argument->GetName()));
    if (sourceVariable && destinationVariable) {
      destinationVariable->setVal(sourceVariable->getVal());
      destinationVariable->setError(sourceVariable->getError());
    }
  };
#if ROOT_VERSION_CODE < ROOT_VERSION(6,30,0)
  std::unique_ptr<TIterator> iterator(source.createIterator());
  RooAbsArg *argument=0;
  while ((argument=static_cast<RooAbsArg *>(iterator->Next()))) {
    assignOne(argument);
  }
#else
  for (RooAbsArg *argument : source) {
    assignOne(argument);
  }
#endif
}

void drawProjection(
    RooPlot *frame,
    RooDataSet &data,
    CategoryModel &model,
    const char *panelTitle,
    bool drawLegend) {
  data.plotOn(frame,Name("plot_data"),Binning(38),MarkerStyle(20),MarkerSize(0.75),LineColor(kBlack));
  model.model->plotOn(frame,Name("plot_total"),LineColor(kBlack),LineWidth(3));
  model.model->plotOn(frame,Components(*model.prompt),Name("plot_prompt"),
                      LineColor(kRed+1),LineWidth(2),LineStyle(kSolid));
  model.model->plotOn(frame,Components(*model.nonprompt),Name("plot_nonprompt"),
                      LineColor(kBlue+1),LineWidth(2),LineStyle(kDashed));
  model.model->plotOn(frame,Components(*model.comb),Name("plot_comb"),
                      LineColor(kGreen+2),LineWidth(2),LineStyle(kDotted));
  frame->SetTitle(panelTitle);
  frame->GetXaxis()->SetTitleSize(0.055);
  frame->GetXaxis()->SetLabelSize(0.045);
  frame->GetYaxis()->SetTitle("Candidates / bin");
  frame->GetYaxis()->SetTitleSize(0.052);
  frame->GetYaxis()->SetLabelSize(0.043);
  frame->Draw();
  if (drawLegend) {
    TLegend *legend=new TLegend(0.53,0.55,0.91,0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.040);
    legend->AddEntry(frame->findObject("plot_data"),"Data","lep");
    legend->AddEntry(frame->findObject("plot_total"),"Total model","l");
    legend->AddEntry(frame->findObject("plot_prompt"),"Prompt","l");
    legend->AddEntry(frame->findObject("plot_nonprompt"),"Nonprompt","l");
    legend->AddEntry(frame->findObject("plot_comb"),"Combinatorial","l");
    legend->Draw();
  }
}

bool drawSlot(
    int slot,
    const char *sourcePath,
    const char *modelPath,
    TPad *massPad,
    TPad *ctauPad) {
  const std::string suffix=std::to_string(slot);
  TFile *source=new TFile(sourcePath,"READ");
  TFile *modelFile=new TFile(modelPath,"READ");
  TTree *tree=dynamic_cast<TTree *>(source->Get("candidates"));
  RooFitResult *fit=dynamic_cast<RooFitResult *>(source->Get("global_yield_fit_result"));
  RooWorkspace *workspace=dynamic_cast<RooWorkspace *>(modelFile->Get("wsp"));
  if (!tree || !fit || !workspace) {
    std::cerr<<"Missing saved tree, yield fit result, or model workspace for J/psi "<<slot<<"\n";
    return false;
  }

  RooRealVar *mass=workspace->var(("Jpsi_mass"+suffix).c_str());
  RooRealVar *ctau=workspace->var(("Jpsi_ctau"+suffix).c_str());
  RooAbsPdf *massSignal=workspace->pdf(("JpsiMassSig"+suffix).c_str());
  RooAbsPdf *massComb=workspace->pdf(("JpsiMassComb"+suffix).c_str());
  RooAbsPdf *ctauPrompt=workspace->pdf(("JpsiCtauSig"+suffix).c_str());
  RooAbsPdf *ctauNonprompt=workspace->pdf(("JpsiCtauBkg"+suffix).c_str());
  RooAbsPdf *ctauComb=workspace->pdf(("JpsiCtauCombBkg"+suffix).c_str());
  if (!mass || !ctau || !massSignal || !massComb || !ctauPrompt || !ctauNonprompt || !ctauComb) {
    std::cerr<<"Missing model component for J/psi "<<slot<<"\n";
    return false;
  }
  mass->setRange(2.95,3.25);
  ctau->setRange(-0.03,0.16);
  RooArgSet observables(*mass,*ctau);
  RooDataSet data(("data_plot_s"+suffix).c_str(),"saved central candidates",
                  observables,Import(*tree));
  CategoryModel model=makeCategoryModel(
      "s"+suffix+"_global_calibrated",data.numEntries(),*workspace,
      *massSignal,*ctauPrompt,*massComb,*ctauNonprompt,*ctauComb);
  std::unique_ptr<RooArgSet> variables(model.model->getVariables());
  assignValues(*variables,fit->constPars());
  assignValues(*variables,fit->floatParsFinal());
  model.massScale->setConstant(true);
  model.ctauScale->setConstant(true);

  massPad->cd();
  massPad->SetLeftMargin(0.13);
  massPad->SetBottomMargin(0.13);
  RooPlot *massFrame=mass->frame();
  massFrame->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
  drawProjection(massFrame,data,model,("J/#psi "+suffix+": mass").c_str(),true);
  massPad->Modified(); massPad->Update();

  ctauPad->cd();
  ctauPad->SetLeftMargin(0.13);
  ctauPad->SetBottomMargin(0.13);
  ctauPad->SetLogy();
  RooPlot *ctauFrame=ctau->frame();
  ctauFrame->GetXaxis()->SetTitle("Decay length [cm]");
  drawProjection(ctauFrame,data,model,("J/#psi "+suffix+": decay length").c_str(),false);
  ctauFrame->SetMinimum(0.5);
  ctauPad->Modified(); ctauPad->Update();
  return true;
}

}  // namespace

void plot_single_jpsi_fit_summary(
    const char *outputPdf="Data_driven/presentation/figures/single_jpsi_fit_summary.pdf",
    const char *outputPng="Data_driven/presentation/figures/single_jpsi_fit_summary.png",
    const char *model="Data_driven/reference/preupdates_20260829/h015/Model_4D_tot.root",
    const char *slot1="Data_driven/reference/preupdates_20260829/h015/jpsi1_sweights.root",
    const char *slot2="Data_driven/reference/preupdates_20260829/h015/jpsi2_sweights.root") {
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.055);
  gStyle->SetPaperSize(15.0,26.666);
  TCanvas canvas("single_jpsi_summary","single-J/psi fit summary",1600,900);
  canvas.SetCanvasSize(1600,900);
  TPad p1("p1","",0.00,0.48,0.50,0.91);
  TPad p2("p2","",0.50,0.48,1.00,0.91);
  TPad p3("p3","",0.00,0.04,0.50,0.47);
  TPad p4("p4","",0.50,0.04,1.00,0.47);
  p1.Draw(); p2.Draw(); p3.Draw(); p4.Draw();
  if (!drawSlot(1,slot1,model,&p1,&p3) || !drawSlot(2,slot2,model,&p2,&p4)) {
    std::cerr<<"Unable to draw saved single-J/psi fit projections\n";
    return;
  }
  canvas.cd();
  TLatex title;
  title.SetNDC(true);
  title.SetTextAlign(22);
  title.SetTextFont(42);
  title.SetTextSize(0.038);
  title.DrawLatex(0.50,0.965,"Prompt single-J/#psi extraction in Data");
  canvas.Modified(); canvas.Update();
  canvas.SaveAs(outputPdf);
  canvas.SaveAs(outputPng);
}
