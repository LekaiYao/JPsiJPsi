import ROOT
from ROOT import TFile, TH2D, TCanvas, gStyle
from array import array
import os
from tdrStyle import setTDRStyle
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
inFileDir = "/eos/user/c/chensh/JPsiJPsi/GEN_nofilter/SPS/CMSSW_10_2_5/src/4mu_acc/direct/"  # chensh Pythia8 GEN samples (local direct/ is empty)
# inFileDir = "direct/"
inFileList = ["SPS_2016_JJ_{}.root".format(i) for i in range(1, 11)]
inFiles = {}
inTrees = []
for inFileName in inFileList:
    inFiles[inFileName] = TFile(inFileDir + inFileName, "READ")
    inTrees.append(inFiles[inFileName].Get("GenAnalyzer/gen_tree"))
for inTree in inTrees:
    print(inTree)

# Define histograms, including efficiency(h) and error(d).
# x-axis: pT, y-axis: y
hEta_Jpsi = TH2D("hEta_Jpsi", "hEta_Jpsi", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
hPt_Jpsi = TH2D("hPt_Jpsi", "hPt_Jpsi", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
dEta_Jpsi = TH2D("dEta_Jpsi", "dEta_Jpsi", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
dPt_Jpsi = TH2D("dPt_Jpsi", "dPt_Jpsi", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
hBin_Jpsi = TH2D("hBin_Jpsi", "hBin_Jpsi", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)# data point distribution(J/psi)

# Define matrices
mBin_Jpsi, mEta_Jpsi, mPt_Jpsi = tuple([np.zeros((len(yBin) - 1, len(ptBin) - 1), dtype=int) for i in range(3)])

# Loop on all rectangular bins, calculate efficiency and fill histograms.
rJpsi1, rJpsi2 = "GEN_pair_id[0].first", "GEN_pair_id[0].second"
rJpsi1_pt, rJpsi1_y, rJpsi2_pt, rJpsi2_y = "GENjpsi_pt[{}]".format(rJpsi1), "GENjpsi_y[{}]".format(rJpsi1), "GENjpsi_pt[{}]".format(rJpsi2), "GENjpsi_y[{}]".format(rJpsi2)
rJpsi1_mupt, rJpsi1_mueta, rJpsi2_mupt, rJpsi2_mueta = "GENjpsi_mu_pt[{}]".format(rJpsi1), "GENjpsi_mu_eta[{}]".format(rJpsi1), "GENjpsi_mu_pt[{}]".format(rJpsi2), "GENjpsi_mu_eta[{}]".format(rJpsi2)
evtMassCut = "GENevt_mass[0] > 7.5 && "
fidCut_Jpsi1 = " && {} > {} && {} < {} && {} > {} && {} < {}".format(rJpsi1_pt, ptBin[0], rJpsi1_pt, ptBin[-1], rJpsi1_y, yBin[0], rJpsi1_y, yBin[-1])
fidCut_Jpsi2 = " && {} > {} && {} < {} && {} > {} && {} < {}".format(rJpsi2_pt, ptBin[0], rJpsi2_pt, ptBin[-1], rJpsi2_y, yBin[0], rJpsi2_y, yBin[-1])
etaCut_Jpsi1, etaCut_Jpsi2 = " && abs({}[0]) < 2.4 && abs({}[1]) < 2.4".format(rJpsi1_mueta, rJpsi1_mueta), " && abs({}[0]) < 2.4 && abs({}[1]) < 2.4".format(rJpsi2_mueta, rJpsi2_mueta)
ptCut_Jpsi1, ptCut_Jpsi2 = " && {}[0] > 3.5 && {}[1] > 3.5".format(rJpsi1_mupt, rJpsi1_mupt), " && {}[0] > 3.5 && {}[1] > 3.5".format(rJpsi2_mupt, rJpsi2_mupt)

nBin_Jpsi_max  = 0
nEvt, nAccEvt = 0, 0
for inTree in inTrees:
    nEvt += inTree.GetEntries(evtMassCut + "1" + fidCut_Jpsi1 + fidCut_Jpsi2)
    nAccEvt += inTree.GetEntries(evtMassCut + "1" + fidCut_Jpsi1 + fidCut_Jpsi2 + etaCut_Jpsi1 + etaCut_Jpsi2 + ptCut_Jpsi1 + ptCut_Jpsi2)
for ptIdx in range(len(ptBin) - 1):
    for yIdy in range(len(yBin) - 1):
        binCut_Jpsi1 = evtMassCut + "{} > {} && {} < {} && {} > {} && {} < {}".format(rJpsi1_pt, ptBin[ptIdx], rJpsi1_pt, ptBin[ptIdx + 1], rJpsi1_y, yBin[yIdy], rJpsi1_y, yBin[yIdy + 1])
        binCut_Jpsi2 = evtMassCut + "{} > {} && {} < {} && {} > {} && {} < {}".format(rJpsi2_pt, ptBin[ptIdx], rJpsi2_pt, ptBin[ptIdx + 1], rJpsi2_y, yBin[yIdy], rJpsi2_y, yBin[yIdy + 1])
        nBin_Jpsi, nEta_Jpsi, nPt_Jpsi = (0.,) * 3
        for inTree in inTrees:
            nBin_Jpsi += inTree.GetEntries(binCut_Jpsi1 + fidCut_Jpsi2) + inTree.GetEntries(binCut_Jpsi2 + fidCut_Jpsi1)
            nEta_Jpsi += inTree.GetEntries(binCut_Jpsi1 + fidCut_Jpsi2 + etaCut_Jpsi1) + inTree.GetEntries(binCut_Jpsi2 + fidCut_Jpsi1 + etaCut_Jpsi2)
            nPt_Jpsi += inTree.GetEntries(binCut_Jpsi1 + fidCut_Jpsi2 + etaCut_Jpsi1 + ptCut_Jpsi1) + inTree.GetEntries(binCut_Jpsi2 + fidCut_Jpsi1 + etaCut_Jpsi2 + ptCut_Jpsi2)
        print("\033[7mProcessing pt-y bin: {}, {}\033[0m".format(ptBin[ptIdx], yBin[yIdy]))
        print("[\033[32mJ/psi\033[0m] PtEta/Eta/Total = %d/%d/%d" % (nPt_Jpsi, nEta_Jpsi, nBin_Jpsi))
        bin = hEta_Jpsi.GetBin(ptIdx + 1, yIdy + 1)
        hBin_Jpsi.SetBinContent(bin, nBin_Jpsi)
        if nBin_Jpsi > nBin_Jpsi_max:
            nBin_Jpsi_max = nBin_Jpsi
        if nBin_Jpsi:
            hEta_Jpsi.SetBinContent(bin, nEta_Jpsi / nBin_Jpsi)
            dEta_Jpsi.SetBinContent(bin, sqrt((nBin_Jpsi - nEta_Jpsi) * nEta_Jpsi / nBin_Jpsi) / nBin_Jpsi)
        if nEta_Jpsi:
            hPt_Jpsi.SetBinContent(bin, nPt_Jpsi / nEta_Jpsi)
            dPt_Jpsi.SetBinContent(bin, sqrt((nEta_Jpsi - nPt_Jpsi) * nPt_Jpsi / nEta_Jpsi) / nEta_Jpsi)
        mBin_Jpsi[yIdy][ptIdx] = int(nBin_Jpsi)
        mEta_Jpsi[yIdy][ptIdx] = int(nEta_Jpsi)
        mPt_Jpsi[yIdy][ptIdx] = int(nPt_Jpsi)

# Plot [histo] with [titles(x, y, z)] and save as [name].
width, height = 3000, 1200
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
    histo.Draw("colzTEXT")
    histo.GetXaxis().SetNdivisions(505)
    c.SaveAs(outPlotDir + name + ".png")
    del c

# Draw and save every plot.
setTDRStyle()
plotSave(hEta_Jpsi, ["pT_Jpsi", "rapidity_Jpsi", "acc_Eta_Jpsi"], "Acc_Eta_Jpsi")
plotSave(hPt_Jpsi, ["pT_Jpsi", "rapidity_Jpsi", "acc_Pt_Jpsi"], "Acc_Pt_Jpsi")
plotSave(dEta_Jpsi, ["pT_Jpsi", "rapidity_Jpsi", "accErr_Eta_Jpsi"], "AccErr_Eta_Jpsi")
plotSave(dPt_Jpsi, ["pT_Jpsi", "rapidity_Jpsi", "accErr_Pt_Jpsi"], "AccErr_Pt_Jpsi")
ROOT.gStyle.SetPaintTextFormat("1.0f");
zRange = nBin_Jpsi_max
plotSave(hBin_Jpsi, ["pT_Jpsi", "rapidity_Jpsi", "N_Jpsi"], "N_distri_Jpsi")

# Output efficiency file format:
# n(number of pT bins), m(number of y bins), nEvt(number of events), nAccEvt(number of accepted events)
# (n+1) double matrix: pT bin boundary
# (m+1) double matrix: y bin boundary
# m*n integer matrix: N1(number of accepted J/psi), N0(number of initial J/psi)
with open(outPlotDir + "acceptance.txt", 'w') as f:
    n, m = len(ptBin) - 1, len(yBin) - 1
    f.write("{} {} {} {}\n".format(n, m, nEvt, nAccEvt))
    l = ""
    for pt in ptBin:
        l += str(pt) + ' '
    l += '\n'
    for y in yBin:
        l += str(y) + ' '
    l += '\n'
    for i in range(m):
        for j in range(n):
            l += str(mPt_Jpsi[i][j]) + ' ' + str(mBin_Jpsi[i][j]) + ' '
        l += '\n'
    f.write(l)

# Extra map file for AN ApxB acceptance figures (keeps acceptance.txt above UNCHANGED for rephrase.cpp).
# Two-step acceptance: A_eta = mEta/mBin (mu |eta|<2.4), A_etapt = mPt/mEta (+ mu pt>3.5).
# Per-bin columns: mEta mPt mBin.
with open(outPlotDir + "acceptance_maps.txt", 'w') as f:
    n, m = len(ptBin) - 1, len(yBin) - 1
    f.write("{} {} {} {}\n".format(n, m, nEvt, nAccEvt))
    l = ""
    for pt in ptBin:
        l += str(pt) + ' '
    l += '\n'
    for y in yBin:
        l += str(y) + ' '
    l += '\n'
    for i in range(m):
        for j in range(n):
            l += str(mEta_Jpsi[i][j]) + ' ' + str(mPt_Jpsi[i][j]) + ' ' + str(mBin_Jpsi[i][j]) + ' '
        l += '\n'
    f.write(l)
