#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooGExpModel.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooProdPdf.h"
#include "RooExtendPdf.h"
#include "RooWorkspace.h"

using namespace std;
using namespace RooFit;

void plot_on(
    RooPlot *frame,
    const RooAddPdf &pdf_all,
    const RooAbsPdf &pdf_P_P, const RooAbsPdf &pdf_P_NP, const RooAbsPdf &pdf_NP_P, const RooAbsPdf &pdf_NP_NP,
    const RooAbsPdf &pdf_Sig_Comb, const RooAbsPdf &pdf_Comb_Sig, const RooAbsPdf &pdf_Comb_Comb,
    const string proj
) {
    pdf_all.plotOn(frame, LineColor(kBlack), LineWidth(2), Name("All"), Normalization(pdf_all.expectedEvents(RooArgSet()), RooAbsReal::NumEvent), ProjectionRange(proj.c_str()));
    pdf_P_P.plotOn(frame, LineColor(kBlue), LineStyle(kSolid), LineWidth(2), Name("P_P"), Normalization((dynamic_cast<RooRealVar &>(pdf_all.coefList()[0])).getVal(), RooAbsReal::NumEvent), ProjectionRange(proj.c_str()));
    pdf_P_NP.plotOn(frame, LineColor(kBlue), LineStyle(kDashed), LineWidth(2), Name("P_NP"), Normalization((dynamic_cast<RooRealVar &>(pdf_all.coefList()[1])).getVal(), RooAbsReal::NumEvent), ProjectionRange(proj.c_str()));
    pdf_NP_P.plotOn(frame, LineColor(kBlue), LineStyle(kDotted), LineWidth(2), Name("NP_P"), Normalization((dynamic_cast<RooRealVar &>(pdf_all.coefList()[2])).getVal(), RooAbsReal::NumEvent), ProjectionRange(proj.c_str()));
    pdf_NP_NP.plotOn(frame, LineColor(kBlue), LineStyle(kDashDotted), LineWidth(2), Name("NP_NP"), Normalization((dynamic_cast<RooRealVar &>(pdf_all.coefList()[3])).getVal(), RooAbsReal::NumEvent), ProjectionRange(proj.c_str()));
    pdf_Sig_Comb.plotOn(frame, LineColor(kRed), LineWidth(1), Name("Sig_Comb"), Normalization((dynamic_cast<RooRealVar &>(pdf_all.coefList()[4])).getVal(), RooAbsReal::NumEvent), ProjectionRange(proj.c_str()));
    pdf_Comb_Sig.plotOn(frame, LineColor(kGreen), LineWidth(1), Name("Comb_Sig"), Normalization((dynamic_cast<RooRealVar &>(pdf_all.coefList()[5])).getVal(), RooAbsReal::NumEvent), ProjectionRange(proj.c_str()));
    pdf_Comb_Comb.plotOn(frame, LineColor(kMagenta), LineWidth(1), Name("Comb_Comb"), Normalization((dynamic_cast<RooRealVar &>(pdf_all.coefList()[6])).getVal(), RooAbsReal::NumEvent), ProjectionRange(proj.c_str()));
    return;
}
void add_entry(RooPlot *frame, TLegend *legend) {
    legend->AddEntry(frame->findObject("All"), "Total p.d.f.", "L");
    legend->AddEntry(frame->findObject("P_P"), "prompt, prompt", "L");
    legend->AddEntry(frame->findObject("P_NP"), "prompt, non-prompt", "L");
    legend->AddEntry(frame->findObject("NP_P"), "non-prompt, prompt", "L");
    legend->AddEntry(frame->findObject("NP_NP"), "non-prompt, non-prompt", "L");
    legend->AddEntry(frame->findObject("Sig_Comb"), "J/#psi, #mu^{+}#mu^{-}", "L");
    legend->AddEntry(frame->findObject("Comb_Sig"), "#mu^{+}#mu^{-}, J/#psi", "L");
    legend->AddEntry(frame->findObject("Comb_Comb"), "#mu^{+}#mu^{-}, #mu^{+}#mu^{-}", "L");
    return;
}

void Plot_4D(
    const RooAbsData *data, const RooAddPdf &pdf_all,
    const RooAbsPdf &pdf_P_P, const RooAbsPdf &pdf_P_NP, const RooAbsPdf &pdf_NP_P, const RooAbsPdf &pdf_NP_NP,
    const RooAbsPdf &pdf_Sig_Comb, const RooAbsPdf &pdf_Comb_Sig, const RooAbsPdf &pdf_Comb_Comb,
    const string &prefix,
    const RooRealVar &Jpsi_mass, const RooRealVar &psi2S_mass, const RooRealVar &Jpsi_ctau, const RooRealVar &psi2S_ctau,
    const string proj = ""
) {
    Int_t BinNum = 100;
    TCanvas *canvas7 = new TCanvas("canvas7", "canvas7", 1500, 1000);
    RooPlot *frame7 = Jpsi_mass.frame(RooFit::Title("Jpsi Mass 4D"), Bins(BinNum));
    data->plotOn(frame7, DataError(RooAbsData::SumW2), Name("Data"));
    plot_on(frame7, pdf_all, pdf_P_P, pdf_P_NP, pdf_NP_P, pdf_NP_NP, pdf_Sig_Comb, pdf_Comb_Sig, pdf_Comb_Comb, proj);
    TLegend *legend7 = new TLegend(.65, .60, .85, .85);
    legend7->AddEntry(frame7->findObject("Data"), "RunII 2018", "L");
    add_entry(frame7, legend7);
    frame7->Draw();
    legend7->DrawClone();
    canvas7->SaveAs((prefix + "4D_JpsiMass.png").c_str());
    TCanvas *canvas8 = new TCanvas("canvas8", "canvas8", 1500, 1000);
    RooPlot *frame8 = psi2S_mass.frame(RooFit::Title("psi2S Mass 4D"), Bins(BinNum));
    data->plotOn(frame8, DataError(RooAbsData::SumW2), Name("Data"));
    plot_on(frame8, pdf_all, pdf_P_P, pdf_P_NP, pdf_NP_P, pdf_NP_NP, pdf_Sig_Comb, pdf_Comb_Sig, pdf_Comb_Comb, proj);
    TLegend *legend8 = new TLegend(.65, .60, .85, .85);
    legend8->AddEntry(frame8->findObject("Data"), "RunII 2018", "L");
    add_entry(frame8, legend8);
    frame8->Draw();
    legend8->DrawClone();
    canvas8->SaveAs((prefix + "4D_psi2SMass.png").c_str());
    TCanvas *canvas9 = new TCanvas("canvas9", "canvas9", 1500, 1000);
    RooPlot *frame9 = Jpsi_ctau.frame(RooFit::Title("Jpsi Ctau 4D"), Bins(BinNum));
    data->plotOn(frame9, DataError(RooAbsData::SumW2), Name("Data"));
    plot_on(frame9, pdf_all, pdf_P_P, pdf_P_NP, pdf_NP_P, pdf_NP_NP, pdf_Sig_Comb, pdf_Comb_Sig, pdf_Comb_Comb, proj);
    TLegend *legend9 = new TLegend(.65, .60, .85, .85);
    legend9->AddEntry(frame9->findObject("Data"), "RunII 2018", "L");
    add_entry(frame9, legend9);
    frame9->SetAxisRange(1, 5e3, "Y");
    gPad->SetLogy();
    frame9->Draw();
    legend9->DrawClone();
    canvas9->SaveAs((prefix + "4D_JpsiCtau.png").c_str());
    TCanvas *canvas0 = new TCanvas("canvas0", "canvas0", 1500, 1000);
    RooPlot *frame0 = psi2S_ctau.frame(RooFit::Title("psi2S Ctau 4D"), Bins(BinNum));
    data->plotOn(frame0, DataError(RooAbsData::SumW2), Name("Data"));
    plot_on(frame0, pdf_all, pdf_P_P, pdf_P_NP, pdf_NP_P, pdf_NP_NP, pdf_Sig_Comb, pdf_Comb_Sig, pdf_Comb_Comb, proj);
    TLegend *legend0 = new TLegend(.65, .60, .85, .85);
    legend0->AddEntry(frame0->findObject("Data"), "RunII 2018", "L");
    add_entry(frame0, legend0);
    frame0->SetAxisRange(1, 5e3, "Y");
    gPad->SetLogy();
    frame0->Draw();
    legend0->DrawClone();
    canvas0->SaveAs((prefix + "4D_psi2SCtau.png").c_str());
    return;
}