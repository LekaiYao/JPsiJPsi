import ROOT
from ROOT import TFile, TH2D, TCanvas
from array import array
import os
from tdrStyle import *
from math import sqrt
import numpy as np

# Specify bin settings for efficiency maps.
ptBin = array('d', [i for i in range(10, 24)] + [24, 26, 28, 30, 35, 40])
yBin = array('d', [-2, -1.75, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 1.75, 2])

# Specify output directory.
outPlotDir = "plot/"
if not os.path.exists(outPlotDir):
    os.system("mkdir " + outPlotDir)

# Handle input files and extract trees.
# inFileDir = "/eos/home-c/chensh/JPsiJPsi/SKIM_tightfilter/SPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc_eff/"
inFileDir = "/eos/home-l/leyao/26JJ/MC_Maker/HelacOnia2016/CMSSW_10_6_20/src/NTUPLE/NLO/"
inFileList = ["Ntuple_2016_SPSNLO_{}.root".format(i) for i in range(1, 2)]
inFiles = {}
inTrees = []
for inFileName in inFileList:
    inFiles[inFileName] = TFile(inFileDir + inFileName, "READ")
    inTrees.append(inFiles[inFileName].Get("rootuple/oniaTree"))
for inTree in inTrees:
    print(inTree)

# Define histograms, including efficiency(h) and error(d).
# x-axis: pT, y-axis: y
hReco_Jpsi = TH2D("hReco_Jpsi", "hReco_Jpsi", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
hId_Jpsi = TH2D("hId_Jpsi", "hId_Jpsi", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
hVtx_Jpsi = TH2D("hVtx_Jpsi", "hVtx_Jpsi", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
dReco_Jpsi = TH2D("dReco_Jpsi", "dReco_Jpsi", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
dId_Jpsi = TH2D("dId_Jpsi", "dId_Jpsi", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
dVtx_Jpsi = TH2D("dVtx_Jpsi", "dVtx_Jpsi", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
hBin_Jpsi = TH2D("hBin_Jpsi", "hBin_Jpsi", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)# data point distribution(J/psi)
# x-axis: J/psi pT, y-axis: psi(2S) pT
hHlt = TH2D("hHlt", "hHlt", len(ptBin) - 1, ptBin, len(ptBin) - 1, ptBin)
hTrg = TH2D("hTrg", "hTrg", len(ptBin) - 1, ptBin, len(ptBin) - 1, ptBin)
dHlt = TH2D("dHlt", "dHlt", len(ptBin) - 1, ptBin, len(ptBin) - 1, ptBin)
dTrg = TH2D("dTrg", "dTrg", len(ptBin) - 1, ptBin, len(ptBin) - 1, ptBin)
hBin_evt = TH2D("hBin_evt", "hBin_evt", len(ptBin) - 1, ptBin, len(ptBin) - 1, ptBin)# data point distribution

# Define matrices
mBin_Jpsi, mRec_Jpsi, mIdt_Jpsi, mVtx_Jpsi = tuple([np.zeros((len(yBin) - 1, len(ptBin) - 1), dtype=int) for i in range(4)])
mVtx_evt, mHlt_evt, mTrg_evt = tuple([np.zeros((len(ptBin) - 1, len(ptBin) - 1), dtype=int) for i in range(3)])

# Loop on all rectangular bins, calculate efficiency and fill histograms.
entryCut = "GEevt_valid && GEevt_passAcc && GEevt_fourMuMass > 7.5 && "
fidJpsi1Cut = "GEJpsi1_pt > {} && GEJpsi1_pt < {} && GEJpsi1_y > {} && GEJpsi1_y < {} && ".format(ptBin[0], ptBin[-1], yBin[0], yBin[-1])
fidJpsi2Cut = "GEJpsi2_pt > {} && GEJpsi2_pt < {} && GEJpsi2_y > {} && GEJpsi2_y < {} && ".format(ptBin[0], ptBin[-1], yBin[0], yBin[-1])
fidYCut = "GEJpsi1_y > {} && GEJpsi1_y < {} && GEJpsi2_y > {} && GEJpsi2_y < {} && ".format(yBin[0], yBin[-1], yBin[0], yBin[-1])
nBin_Jpsi_max, nBin_evt_max = 0, 0
nAccEvt = 0
for inTree in inTrees:
    nAccEvt += inTree.GetEntries(entryCut + fidJpsi1Cut + fidJpsi2Cut + "1")
for ptIdx in range(len(ptBin) - 1):
    for yIdy in range(len(yBin) - 1):
        binCut_Jpsi1 = entryCut + fidJpsi2Cut + "GEJpsi1_pt > {} && GEJpsi1_pt < {} && GEJpsi1_y > {} && GEJpsi1_y < {}".format(ptBin[ptIdx], ptBin[ptIdx + 1], yBin[yIdy], yBin[yIdy + 1])
        binCut_Jpsi2 = entryCut + fidJpsi1Cut + "GEJpsi2_pt > {} && GEJpsi2_pt < {} && GEJpsi2_y > {} && GEJpsi2_y < {}".format(ptBin[ptIdx], ptBin[ptIdx + 1], yBin[yIdy], yBin[yIdy + 1])
        nBin_Jpsi, nReco_Jpsi, nId_Jpsi, nVtx_Jpsi = (0.,) * 4
        for inTree in inTrees:
            nBin_Jpsi += inTree.GetEntries(binCut_Jpsi1) + inTree.GetEntries(binCut_Jpsi2)
            nReco_Jpsi += inTree.GetEntries(binCut_Jpsi1 + " && GEJpsi1_matchGEN && REmu_pt[3]") + inTree.GetEntries(binCut_Jpsi2 + " && GEJpsi2_matchGEN && REmu_pt[3]")
            nId_Jpsi += inTree.GetEntries(binCut_Jpsi1 + " && GEJpsi1_passID && REmu_pt[3]") + inTree.GetEntries(binCut_Jpsi2 + " && GEJpsi2_passID && REmu_pt[3]")
            nVtx_Jpsi += inTree.GetEntries(binCut_Jpsi1 + " && GEJpsi1_passVtx") + inTree.GetEntries(binCut_Jpsi2 + " && GEJpsi2_passVtx")
        print("\033[7mProcessing pt-y bin: {}, {}\033[0m".format(ptBin[ptIdx], yBin[yIdy]))
        print("[\033[32mJ/psi\033[0m] Vtx/Id/Reco/Total = %d/%d/%d/%d" % (nVtx_Jpsi, nId_Jpsi, nReco_Jpsi, nBin_Jpsi))
        bin = hReco_Jpsi.GetBin(ptIdx + 1, yIdy + 1)
        hBin_Jpsi.SetBinContent(bin, nBin_Jpsi)
        if nBin_Jpsi > nBin_Jpsi_max:
            nBin_Jpsi_max = nBin_Jpsi
        if nBin_Jpsi:
            hReco_Jpsi.SetBinContent(bin, nReco_Jpsi / nBin_Jpsi)
            dReco_Jpsi.SetBinContent(bin, sqrt((nBin_Jpsi - nReco_Jpsi) * nReco_Jpsi / nBin_Jpsi) / nBin_Jpsi)
        if nReco_Jpsi:
            hId_Jpsi.SetBinContent(bin, nId_Jpsi / nReco_Jpsi)
            dId_Jpsi.SetBinContent(bin, sqrt((nReco_Jpsi - nId_Jpsi) * nId_Jpsi / nReco_Jpsi) / nReco_Jpsi)
        if nId_Jpsi:
            hVtx_Jpsi.SetBinContent(bin, nVtx_Jpsi / nId_Jpsi)
            dVtx_Jpsi.SetBinContent(bin, sqrt((nId_Jpsi - nVtx_Jpsi) * nVtx_Jpsi / nId_Jpsi) / nId_Jpsi)
        mBin_Jpsi[yIdy][ptIdx] = int(nBin_Jpsi)
        mRec_Jpsi[yIdy][ptIdx] = int(nReco_Jpsi)
        mIdt_Jpsi[yIdy][ptIdx] = int(nId_Jpsi)
        mVtx_Jpsi[yIdy][ptIdx] = int(nVtx_Jpsi)
    for ptIdy in range(len(ptBin) - 1):
        binCut = entryCut + fidYCut + "GEJpsi1_pt > {} && GEJpsi1_pt < {} && GEJpsi2_pt > {} && GEJpsi2_pt < {}".format(ptBin[ptIdx], ptBin[ptIdx + 1], ptBin[ptIdy], ptBin[ptIdy + 1])
        nVtx, nHlt, nTrg = (0.,) * 3
        for inTree in inTrees:
            nVtx += inTree.GetEntries(binCut + " && GEJpsi1_passVtx && GEJpsi2_passVtx")
            nHlt += inTree.GetEntries(binCut + " && GEevt_passHLT")
            nTrg += inTree.GetEntries(binCut + " && GEevt_matchTrg")
        print("\033[7mProcessing pt-pt bin: {}, {}\033[0m".format(ptBin[ptIdx], ptBin[ptIdy]))
        print("[\033[36mevent\033[0m] Trg/Hlt/Vtx = %d/%d/%d" % (nTrg, nHlt, nVtx))
        bin = hHlt.GetBin(ptIdx + 1, ptIdy + 1)
        hBin_evt.SetBinContent(bin, nVtx)
        if nVtx > nBin_evt_max:
            nBin_evt_max = nVtx
        if nVtx:
            hHlt.SetBinContent(bin, nHlt / nVtx)
            dHlt.SetBinContent(bin, sqrt((nVtx - nHlt) * nHlt / nVtx) / nVtx)
        if nHlt:
            hTrg.SetBinContent(bin, nTrg / nHlt)
            dTrg.SetBinContent(bin, sqrt((nHlt - nTrg) * nTrg / nHlt) / nHlt)
        mVtx_evt[ptIdy][ptIdx] = int(nVtx)
        mHlt_evt[ptIdy][ptIdx] = int(nHlt)
        mTrg_evt[ptIdy][ptIdx] = int(nTrg)

# Plot [histo] with [titles(x, y, z)] and save as [name].
width, height = 3000, 1200
squarish = False
zRange = 1
def plotSave(histo, titles, name):
    c = TCanvas(name, name, width, height)
    c.cd()
    c.SetTopMargin(0.10)
    c.SetRightMargin(0.20)
    histo.GetXaxis().SetTitle(titles[0])
    histo.GetYaxis().SetTitle(titles[1])
    histo.GetZaxis().SetTitle(titles[2])
    histo.GetXaxis().SetTitleOffset(1.0)
    histo.GetYaxis().SetTitleOffset(0.8)
    histo.GetZaxis().SetTitleOffset(0.6)
    histo.GetZaxis().SetRangeUser(0., zRange) 
    histo.Draw("colzTEXT" + ("45" if squarish else ""))
    histo.GetXaxis().SetNdivisions(505)
    if squarish:
        histo.GetYaxis().SetNdivisions(505)
    if width == height:
        histo.GetYaxis().SetNdivisions(505)
    c.SaveAs(outPlotDir + name + ".png")
    del c

# Draw and save every plot.
setTDRStyle()
plotSave(hReco_Jpsi, ["pT_Jpsi", "rapidity_Jpsi", "eff_RECO_Jpsi"], "Eff_RECO_Jpsi")
plotSave(hId_Jpsi, ["pT_Jpsi", "rapidity_Jpsi", "eff_ID_Jpsi"], "Eff_ID_Jpsi")
plotSave(hVtx_Jpsi, ["pT_Jpsi", "rapidity_Jpsi", "eff_Dimuon_Jpsi"], "Eff_Dimuon_Jpsi")
plotSave(dReco_Jpsi, ["pT_Jpsi", "rapidity_Jpsi", "effErr_RECO_Jpsi"], "EffErr_RECO_Jpsi")
plotSave(dId_Jpsi, ["pT_Jpsi", "rapidity_Jpsi", "effErr_ID_Jpsi"], "EffErr_ID_Jpsi")
plotSave(dVtx_Jpsi, ["pT_Jpsi", "rapidity_Jpsi", "effErr_Dimuon_Jpsi"], "EffErr_Dimuon_Jpsi")
ROOT.gStyle.SetPaintTextFormat("1.0f");
zRange = nBin_Jpsi_max
plotSave(hBin_Jpsi, ["pT_Jpsi", "rapidity_Jpsi", "N_Jpsi"], "N_distri_Jpsi")
squarish, height = True, 2500
zRange = nBin_evt_max
plotSave(hBin_evt, ["pT_Jpsi", "pT_psi2S", "N_event"], "N_distri_dimuon_pT")
ROOT.gStyle.SetPaintTextFormat("1.2f");
zRange = 1.
plotSave(hHlt, ["pT_Jpsi", "pT_psi2S", "eff_HLT"], "Eff_HLT")
plotSave(dHlt, ["pT_Jpsi", "pT_psi2S", "effErr_HLT"], "EffErr_HLT")
plotSave(hTrg, ["pT_Jpsi", "pT_psi2S", "eff_FourMu"], "Eff_FourMu")
plotSave(dTrg, ["pT_Jpsi", "pT_psi2S", "effErr_FourMu"], "EffErr_FourMu")

# Output efficiency file format:
# n(number of pT bins), m(number of y bins), nAccEvt(number of events), n_Jpsi(number of J/psi = 2 * number of events)
# (n+1) double matrix: pT bin boundary
# (m+1) double matrix: y bin boundary
# m*(2n) integer matrix: N1(number of J/psi with vertex), N0(number of accepted J/psi)
# n*(2n) integer matrix: N1(number of trigger-matched events), N0(number of events with vertices)
with open(outPlotDir + "raw_efficiency_NLO.txt", 'w') as f:
    n_Jpsi = 0
    n, m = len(ptBin) - 1, len(yBin) - 1
    for i in range(m):
        for j in range(n):
            n_Jpsi += mBin_Jpsi[i][j]
    f.write("{} {} {} {}\n".format(n, m, nAccEvt, n_Jpsi))
    l = ""
    for pt in ptBin:
        l += str(pt) + ' '
    l += '\n'
    for y in yBin:
        l += str(y) + ' '
    l += '\n'
    for i in range(m):
        for j in range(n):
            l += str(mVtx_Jpsi[i][j]) + ' ' + str(mIdt_Jpsi[i][j]) + ' ' + str(mRec_Jpsi[i][j]) + ' ' + str(mBin_Jpsi[i][j]) + ' '
        l += '\n'
    for i in range(n):
        for j in range(n):
            l += str(mTrg_evt[i][j]) + ' ' + str(mHlt_evt[i][j]) + ' ' + str(mVtx_evt[i][j]) + ' '
        l += '\n'
    f.write(l)
