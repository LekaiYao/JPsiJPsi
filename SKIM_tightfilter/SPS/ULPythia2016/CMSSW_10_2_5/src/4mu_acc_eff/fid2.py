from ROOT import TFile

# Keep the same input settings as efficiencyPlot.py
inFileDir = "/eos/home-c/chensh/JPsiJPsi/SKIM_tightfilter/SPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc_eff/"
inFileList = ["Ntuple_2016_SPS_{}.root".format(i) for i in range(1, 11)]

inFiles = {}
inTrees = []
for inFileName in inFileList:
    inFiles[inFileName] = TFile(inFileDir + inFileName, "READ")
    inTrees.append(inFiles[inFileName].Get("rootuple/oniaTree"))

# entryCut in efficiencyPlot.py:
# GEevt_valid && GEevt_passAcc && GEevt_fourMuMass > 7.5 && 
baseCut = "GEevt_valid && GEevt_passAcc"
entryCut = baseCut + " && GEevt_fourMuMass > 7.5"

n1 = 0  # pass full entryCut
n2 = 0  # pass entryCut without GEevt_fourMuMass > 7.5
for inTree in inTrees:
    n1 += inTree.GetEntries(entryCut)
    n2 += inTree.GetEntries(baseCut)

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
