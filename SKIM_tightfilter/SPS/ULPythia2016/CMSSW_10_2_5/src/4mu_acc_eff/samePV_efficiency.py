#!/usr/bin/env python3
from array import array
from math import sqrt
import os

import ROOT
from ROOT import TFile, TH2D, TCanvas, TLorentzVector

from tdrStyle import setTDRStyle

# Use the same pT binning style as efficiencyPlot.py
ptBin = array("d", [i for i in range(10, 24)] + [24, 26, 28, 30, 35, 40])

# Input settings (same files as efficiencyPlot.py)
inFileDir = "/eos/home-l/leyao/26JJ/MC_Maker/HelacOnia2016/CMSSW_10_6_20/src/NTUPLE/NLO/"
inFileList = ["Ntuple_2016_SPSNLO_{}.root".format(i) for i in range(1, 2)]

outPlotDir = "plot/"
if not os.path.exists(outPlotDir):
    os.system("mkdir " + outPlotDir)


def plot_save(histo, titles, name, z_range_max=1.0):
    c = TCanvas(name, name, 2500, 2500)
    c.cd()
    c.SetTopMargin(0.10)
    c.SetRightMargin(0.20)
    histo.GetXaxis().SetTitle(titles[0])
    histo.GetYaxis().SetTitle(titles[1])
    histo.GetZaxis().SetTitle(titles[2])
    histo.GetXaxis().SetTitleOffset(1.0)
    histo.GetYaxis().SetTitleOffset(0.8)
    histo.GetZaxis().SetTitleOffset(0.6)
    histo.GetZaxis().SetRangeUser(0.0, z_range_max)
    histo.GetXaxis().SetNdivisions(505)
    histo.GetYaxis().SetNdivisions(505)
    histo.Draw("colzTEXT45")
    c.SaveAs(outPlotDir + name + ".png")
    del c


def in_pt_range(pt):
    return (pt > 10.0) and (pt < 40.0)


def main():
    # Numerator/denominator for samePV efficiency
    hDen = TH2D("hDen_samePV", "hDen_samePV", len(ptBin) - 1, ptBin, len(ptBin) - 1, ptBin)
    hNum = TH2D("hNum_samePV", "hNum_samePV", len(ptBin) - 1, ptBin, len(ptBin) - 1, ptBin)
    hEff = TH2D("hEff_samePV", "hEff_samePV", len(ptBin) - 1, ptBin, len(ptBin) - 1, ptBin)
    hErr = TH2D("hErr_samePV", "hErr_samePV", len(ptBin) - 1, ptBin, len(ptBin) - 1, ptBin)

    for inFileName in inFileList:
        f = TFile(inFileDir + inFileName, "READ")
        t = f.Get("rootuple/oniaTree")
        if not t:
            print("Skip missing tree in", inFileName)
            continue

        n_evt = t.GetEntries()
        print("Processing", inFileName, "entries =", n_evt)

        for i in range(n_evt):
            t.GetEntry(i)

            # Add all event-level fiducial/trigger cuts used in efficiencyPlot.py chain.
            if (not bool(t.GEevt_valid)) or (not bool(t.GEevt_passAcc)):
                continue
            if float(t.GEevt_fourMuMass) <= 7.5:
                continue
            if not (10.0 < float(t.GEJpsi1_pt) < 40.0 and 10.0 < float(t.GEJpsi2_pt) < 40.0):
                continue
            if not (-2.0 < float(t.GEJpsi1_y) < 2.0 and -2.0 < float(t.GEJpsi2_y) < 2.0):
                continue
            if (not bool(t.GEevt_passHLT)) or (not bool(t.GEevt_matchTrg)):
                continue
            if (not bool(t.GEJpsi1_passVtx)) or (not bool(t.GEJpsi2_passVtx)):
                continue

            evt_match = t.REevt_matchTrg
            evt_hlt = t.REevt_passHLT
            evt_samepv = t.REevt_samePV
            evt_id1 = t.REevt_JpsiId1
            evt_id2 = t.REevt_JpsiId2

            for j in range(len(evt_match) - 1, -1, -1):
                # all cuts minus samePV: HLT + matchTrg + J/psi pT window + 4mu mass > 7.5
                if not evt_hlt[j]:
                    continue
                if not evt_match[j]:
                    continue

                id1 = int(evt_id1[j])
                id2 = int(evt_id2[j])
                pt1 = float(t.REJpsi_pt[id1])
                pt2 = float(t.REJpsi_pt[id2])
                if (not in_pt_range(pt1)) or (not in_pt_range(pt2)):
                    continue

                lv1 = TLorentzVector()
                lv2 = TLorentzVector()
                lv1.SetPtEtaPhiM(float(t.REJpsi_pt[id1]), float(t.REJpsi_eta[id1]), float(t.REJpsi_phi[id1]), float(t.REJpsi_mass[id1]))
                lv2.SetPtEtaPhiM(float(t.REJpsi_pt[id2]), float(t.REJpsi_eta[id2]), float(t.REJpsi_phi[id2]), float(t.REJpsi_mass[id2]))
                if (lv1 + lv2).M() < 7.5:
                    continue

                hDen.Fill(pt1, pt2)
                if evt_samepv[j]:
                    hNum.Fill(pt1, pt2)

                # Match rephrase.cpp behavior: keep first accepted candidate per event.
                break

        f.Close()

    # Build efficiency and binomial uncertainty
    for ix in range(1, len(ptBin)):
        for iy in range(1, len(ptBin)):
            n_den = hDen.GetBinContent(ix, iy)
            n_num = hNum.GetBinContent(ix, iy)
            if n_den <= 0:
                continue
            eff = n_num / n_den
            err = sqrt(eff * (1.0 - eff) / n_den)
            hEff.SetBinContent(ix, iy, eff)
            hErr.SetBinContent(ix, iy, err)

    setTDRStyle()
    ROOT.gStyle.SetPaintTextFormat("1.3f")
    plot_save(hEff, ["pT(J/psi1) [GeV]", "pT(J/psi2) [GeV]", "eff_samePV"], "Eff_samePV", 1.0)
    plot_save(hErr, ["pT(J/psi1) [GeV]", "pT(J/psi2) [GeV]", "effErr_samePV"], "EffErr_samePV", 0.5)
    ROOT.gStyle.SetPaintTextFormat("1.0f")
    plot_save(hDen, ["pT(J/psi1) [GeV]", "pT(J/psi2) [GeV]", "N(all cuts - samePV)"], "N_denom_samePV", max(1.0, hDen.GetMaximum()))
    plot_save(hNum, ["pT(J/psi1) [GeV]", "pT(J/psi2) [GeV]", "N(all cuts)"], "N_numer_samePV", max(1.0, hNum.GetMaximum()))

    print("Saved plots to", outPlotDir)


if __name__ == "__main__":
    main()
