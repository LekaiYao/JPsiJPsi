#include "TCanvas.h"
#include "TAxis.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "RVersion.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooHist.h"
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
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 40, 0)
    RooArgSet normalizationSet;
    const double expectedEvents = pdf_all.expectedEvents(&normalizationSet);
#else
    const double expectedEvents = pdf_all.expectedEvents(RooArgSet());
#endif
    pdf_all.plotOn(frame, LineColor(kBlack), LineWidth(2), Name("All"), Normalization(expectedEvents, RooAbsReal::NumEvent), ProjectionRange(proj.c_str()));
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
    legend->AddEntry(frame->findObject("Data"), "RunII 2016", "L");
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
void add_latex(TLatex &latex) {
    latex.SetTextSize(0.05);
    latex.SetTextFont(62);
    latex.DrawLatex(0.18, 0.85, "CMS");
    latex.SetTextSize(0.04);
    latex.SetTextFont(52);
    latex.DrawLatex(0.27, 0.85, "Preliminary");
    latex.SetTextSize(0.03);
    latex.SetTextFont(42);
    latex.DrawLatex(0.74, 0.91, "36.31 fb^{-1} (13TeV)");
    return;
}
void set_pull_style(RooPlot* frame_pull) {
    frame_pull->SetTitle("");
    frame_pull->GetYaxis()->SetTitleSize(0.1);
    frame_pull->GetYaxis()->SetLabelSize(0.1);
    frame_pull->GetXaxis()->SetTitleSize(0.1);
    frame_pull->GetXaxis()->SetLabelSize(0.1);
    return;
}
void set_line_style(TLine* l) {
    l->SetLineWidth(2);
    l->SetLineColor(kRed);
    l->SetLineStyle(kDashed);
}
void handle_plot(TCanvas *canvas, RooPlot *frame, RooPlot *frame_pull, string title, bool logy=false) {
    RooHist* pull = frame->pullHist("Data", "All");
    frame_pull->addPlotable(pull);
    TLegend *legend = new TLegend(.65, .50, .9, .9);
    add_entry(frame, legend);
    canvas->Divide(1, 2);
    canvas->cd(1)->SetPad(0.01, 0.20, 0.99, 0.99);
    gPad->SetLeftMargin(0.15);
    if(logy) {
        frame->SetAxisRange(1, 5*frame->GetMaximum(), "Y");
        gPad->SetLogy();
    }
    frame->SetTitle("");
    frame->Draw();
    legend->DrawClone();
    TLatex latex;
    latex.SetNDC();
    add_latex(latex);
    canvas->cd(2)->SetPad(0.01, 0.03, 0.99, 0.25);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.25);
    set_pull_style(frame_pull);
    frame_pull->GetXaxis()->SetTitle(title.c_str());
    frame_pull->Draw();
    TLine* l = new TLine(frame->GetXaxis()->GetXmin(), 0, frame->GetXaxis()->GetXmax(), 0);
    set_line_style(l);
    l->Draw("same");
    return;
}

void Plot_4D(
    const RooAbsData *data, const RooAddPdf &pdf_all,
    const RooAbsPdf &pdf_P_P, const RooAbsPdf &pdf_P_NP, const RooAbsPdf &pdf_NP_P, const RooAbsPdf &pdf_NP_NP,
    const RooAbsPdf &pdf_Sig_Comb, const RooAbsPdf &pdf_Comb_Sig, const RooAbsPdf &pdf_Comb_Comb,
    const string &prefix,
    const RooRealVar &Jpsi_mass1, const RooRealVar &Jpsi_mass2, const RooRealVar &Jpsi_ctau1, const RooRealVar &Jpsi_ctau2,
    const string proj = "", const string &suffix = "png"
) {
    Int_t BinNum = 100;
    TCanvas *canvas7 = new TCanvas("canvas7", "canvas7", 1500, 1500);
    RooPlot *frame7 = Jpsi_mass1.frame(RooFit::Title("Jpsi Mass 4D"), Bins(BinNum)), *frame_pull7 = Jpsi_mass1.frame(RooFit::Title("Pull"), Bins(BinNum));
    data->plotOn(frame7, DataError(RooAbsData::SumW2), Name("Data"));
    plot_on(frame7, pdf_all, pdf_P_P, pdf_P_NP, pdf_NP_P, pdf_NP_NP, pdf_Sig_Comb, pdf_Comb_Sig, pdf_Comb_Comb, proj);
    handle_plot(canvas7, frame7, frame_pull7, "M(J/#psi_{1})[GeV]");
    canvas7->SaveAs((prefix + "4D_JpsiMass1." + suffix).c_str());
    TCanvas *canvas8 = new TCanvas("canvas8", "canvas8", 1500, 1500);
    RooPlot *frame8 = Jpsi_mass2.frame(RooFit::Title("psi2S Mass 4D"), Bins(BinNum)), *frame_pull8 = Jpsi_mass2.frame(RooFit::Title("Pull"), Bins(BinNum));
    data->plotOn(frame8, DataError(RooAbsData::SumW2), Name("Data"));
    plot_on(frame8, pdf_all, pdf_P_P, pdf_P_NP, pdf_NP_P, pdf_NP_NP, pdf_Sig_Comb, pdf_Comb_Sig, pdf_Comb_Comb, proj);
    handle_plot(canvas8, frame8, frame_pull8, "M(J/#psi_{2})[GeV]");
    canvas8->SaveAs((prefix + "4D_JpsiMass2." + suffix).c_str());
    TCanvas *canvas9 = new TCanvas("canvas9", "canvas9", 1500, 1500);
    RooPlot *frame9 = Jpsi_ctau1.frame(RooFit::Title("Jpsi Ctau 4D"), Bins(BinNum)), *frame_pull9 = Jpsi_ctau1.frame(RooFit::Title("Pull"), Bins(BinNum));
    data->plotOn(frame9, DataError(RooAbsData::SumW2), Name("Data"));
    plot_on(frame9, pdf_all, pdf_P_P, pdf_P_NP, pdf_NP_P, pdf_NP_NP, pdf_Sig_Comb, pdf_Comb_Sig, pdf_Comb_Comb, proj);
    handle_plot(canvas9, frame9, frame_pull9, "c#tau(J/#psi_{1})[GeV]", true);
    canvas9->SaveAs((prefix + "4D_JpsiCtau1." + suffix).c_str());
    TCanvas *canvas0 = new TCanvas("canvas0", "canvas0", 1500, 1500);
    RooPlot *frame0 = Jpsi_ctau2.frame(RooFit::Title("psi2S Ctau 4D"), Bins(BinNum)), *frame_pull0 = Jpsi_ctau2.frame(RooFit::Title("Pull"), Bins(BinNum));
    data->plotOn(frame0, DataError(RooAbsData::SumW2), Name("Data"));
    plot_on(frame0, pdf_all, pdf_P_P, pdf_P_NP, pdf_NP_P, pdf_NP_NP, pdf_Sig_Comb, pdf_Comb_Sig, pdf_Comb_Comb, proj);
    handle_plot(canvas0, frame0, frame_pull0, "c#tau(J/#psi_{2})[GeV]", true);
    canvas0->SaveAs((prefix + "4D_JpsiCtau2." + suffix).c_str());
    return;
}
