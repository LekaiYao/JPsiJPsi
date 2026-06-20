from ROOT import *
from array import array
import os
from tdrStyle import *
from math import sqrt
import numpy as np

xLO, xSPS = 0, 1# The fraction of LO xsec and SPS xsec
aLO, aNLO, aDPS = 0.23713, 0.23713, 0.18933# Gross acceptances
# Acceptance of Pythia8 generated SPS events: 0.33072
rawLO, rawNLO, rawDPS = ("./plot/raw_efficiency_NLO.txt",
    "./plot/raw_efficiency_NLO.txt",
    # "./raw_efficiency_backup.txt",
    "/eos/home-c/chensh/JPsiJPsi/SKIM_tightfilter/DPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc_eff/plot/raw_efficiency.txt")# Paths to raw efficiency maps
ptBin = array('d', [i for i in range(10, 24)] + [24, 26, 28, 30, 35, 40])
yBin = array('d', [-2, -1.75, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 1.75, 2])

mVtx_Bin_LO, mTrg_Vtx_LO, mVtx_Bin_NLO, mTrg_Vtx_NLO, mVtx_Bin_DPS, mTrg_Vtx_DPS = ([] for i in range(6))
# Read LO matrix
with open(rawLO, 'r') as f:
    l = f.readline()
    nLO = float(l.split(' ')[2])
    for i in range(2):
        l = f.readline()
    for i in range(len(yBin) - 1):
        l = f.readline()
        mVtx_Bin_LO.append(np.array(l.strip().split(), dtype=float) / nLO * xSPS * xLO * aLO)
    for i in range(len(ptBin) - 1):
        l = f.readline()
        mTrg_Vtx_LO.append(np.array(l.strip().split(), dtype=float) / nLO * xSPS * xLO * aLO)
mVtx_Bin_LO, mTrg_Vtx_LO = np.array(mVtx_Bin_LO), np.array(mTrg_Vtx_LO)
# Read NLO matrix
with open(rawNLO, 'r') as f:
    l = f.readline()
    nNLO = float(l.split(' ')[2])
    for i in range(2):
        l = f.readline()
    for i in range(len(yBin) - 1):
        l = f.readline()
        mVtx_Bin_NLO.append(np.array(l.strip().split(), dtype=float) / nNLO * xSPS * (1 - xLO) * aNLO)
    for i in range(len(ptBin) - 1):
        l = f.readline()
        mTrg_Vtx_NLO.append(np.array(l.strip().split(), dtype=float) / nNLO * xSPS * (1 - xLO) * aNLO)
mVtx_Bin_NLO, mTrg_Vtx_NLO = np.array(mVtx_Bin_NLO), np.array(mTrg_Vtx_NLO)
# Read DPS matrix
with open(rawDPS, 'r') as f:
    l = f.readline()
    nDPS = float(l.split(' ')[2])
    for i in range(2):
        l = f.readline()
    for i in range(len(yBin) - 1):
        l = f.readline()
        mVtx_Bin_DPS.append(np.array(l.strip().split(), dtype=float) / nDPS * (1 - xSPS) * aDPS)
    for i in range(len(ptBin) - 1):
        l = f.readline()
        mTrg_Vtx_DPS.append(np.array(l.strip().split(), dtype=float) / nDPS * (1 - xSPS) * aDPS)
mVtx_Bin_DPS, mTrg_Vtx_DPS = np.array(mVtx_Bin_DPS), np.array(mTrg_Vtx_DPS)

mVtx_Bin, mTrg_Vtx = mVtx_Bin_LO + mVtx_Bin_NLO + mVtx_Bin_DPS, mTrg_Vtx_LO + mTrg_Vtx_NLO + mTrg_Vtx_DPS
mVtx_Bin = mVtx_Bin.reshape(mVtx_Bin.shape[0], mVtx_Bin.shape[1] // 4, 4)
mTrg_Vtx = mTrg_Vtx.reshape(mTrg_Vtx.shape[0], mTrg_Vtx.shape[1] // 3, 3)

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

nBin_Jpsi_max, nBin_evt_max = 0, 0
for ptIdx in range(len(ptBin) - 1):
    for yIdy in range(len(yBin) - 1):
        bin = hReco_Jpsi.GetBin(ptIdx + 1, yIdy + 1)
        nVtx_Jpsi, nId_Jpsi, nReco_Jpsi, nBin_Jpsi = tuple(mVtx_Bin[yIdy, ptIdx])
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
    for ptIdy in range(len(ptBin) - 1):
        bin = hHlt.GetBin(ptIdx + 1, ptIdy + 1)
        nTrg, nHlt, nVtx = tuple(mTrg_Vtx[ptIdy, ptIdx])
        hBin_evt.SetBinContent(bin, nVtx)
        if nVtx > nBin_evt_max:
            nBin_evt_max = nVtx
        if nVtx:
            hHlt.SetBinContent(bin, nHlt / nVtx)
            dHlt.SetBinContent(bin, sqrt((nVtx - nHlt) * nHlt / nVtx) / nVtx)
        if nHlt:
            hTrg.SetBinContent(bin, nTrg / nHlt)
            dTrg.SetBinContent(bin, sqrt((nHlt - nTrg) * nTrg / nHlt) / nHlt)

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
    c.SaveAs("plot/" + name + ".png")
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
zRange = 1
plotSave(hHlt, ["pT_Jpsi", "pT_psi2S", "eff_HLT"], "Eff_HLT")
plotSave(dHlt, ["pT_Jpsi", "pT_psi2S", "effErr_HLT"], "EffErr_HLT")
plotSave(hTrg, ["pT_Jpsi", "pT_psi2S", "eff_FourMu"], "Eff_FourMu")
plotSave(dTrg, ["pT_Jpsi", "pT_psi2S", "effErr_FourMu"], "EffErr_FourMu")

with open("txt/efficiency_{}_{}.txt".format(xLO, xSPS), 'w') as f:
    n, m = len(ptBin) - 1, len(yBin) - 1
    f.write("{} {}\n".format(n, m))
    l = ""
    for pt in ptBin:
        l += str(pt) + ' '
    l += '\n'
    for y in yBin:
        l += str(y) + ' '
    l += '\n'
    for i in range(m):
        for j in range(n):
            l += str(mVtx_Bin[i][j][0]) + ' ' + str(mVtx_Bin[i][j][3]) + ' '
        l += '\n'
    for i in range(n):
        for j in range(n):
            l += str(mTrg_Vtx[i][j][0]) + ' ' + str(mTrg_Vtx[i][j][2]) + ' '
        l += '\n'
    f.write(l)