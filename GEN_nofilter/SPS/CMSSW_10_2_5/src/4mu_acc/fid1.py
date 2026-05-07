from ROOT import TFile
import os

# Keep the same input settings as efficiencyPlot.py
inFileDir = "/eos/home-c/chensh/JPsiJPsi/GEN_nofilter/SPS/CMSSW_10_2_5/src/4mu_acc/direct/"
inFileList = ["SPS_2016_JJ_{}.root".format(i) for i in range(1, 2)]

inFiles = {}
inTrees = []
for inFileName in inFileList:
    inFiles[inFileName] = TFile(inFileDir + inFileName, "READ")
    inTrees.append(inFiles[inFileName].Get("GenAnalyzer/gen_tree"))

# Build branch expressions in the same style as efficiencyPlot.py
rJpsi1, rJpsi2 = "GEN_pair_id[0].first", "GEN_pair_id[0].second"
rJpsi1_pt = "GENjpsi_pt[{}]".format(rJpsi1)
rJpsi2_pt = "GENjpsi_pt[{}]".format(rJpsi2)
rJpsi1_y = "GENjpsi_y[{}]".format(rJpsi1)
rJpsi2_y = "GENjpsi_y[{}]".format(rJpsi2)

# Event-level cuts:
# N1: both J/psi satisfy |y| < 1.2 and 20 < pT < 40
# N2: both J/psi satisfy |y| < 2.0 and 10 < pT < 40
cut_n1 = (
    "abs({y1}) < 1.2 && abs({y2}) < 1.2 && "
    "{pt1} > 20 && {pt1} < 40 && {pt2} > 20 && {pt2} < 40"
).format(y1=rJpsi1_y, y2=rJpsi2_y, pt1=rJpsi1_pt, pt2=rJpsi2_pt)

cut_n2 = (
    "abs({y1}) < 2.0 && abs({y2}) < 2.0 && "
    "{pt1} > 10 && {pt1} < 40 && {pt2} > 10 && {pt2} < 40"
).format(y1=rJpsi1_y, y2=rJpsi2_y, pt1=rJpsi1_pt, pt2=rJpsi2_pt)

n1 = 0
n2 = 0
for inTree in inTrees:
    n1 += inTree.GetEntries(cut_n1)
    n2 += inTree.GetEntries(cut_n2)

r = float(n1) / float(n2) if n2 else 0.0

outFile = "fid1.txt"
with open(outFile, "w") as f:
    f.write("N1 {}\n".format(n1))
    f.write("N2 {}\n".format(n2))
    f.write("r {:.12g}\n".format(r))

print("N1 =", n1)
print("N2 =", n2)
print("r =", r)
print("Saved to", outFile)
