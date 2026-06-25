# Fast single-pass acceptance recompute (TTree::Draw instead of GetEntries-per-bin).
# Uses the EXACT same cut expressions as the two 4mu_acc/efficiencyPlot.py, so it
# reproduces mPt/mBin in the existing acceptance.txt and adds the missing mEta needed
# for the two-step AN acc maps. Writes <dir>/plot/acceptance_maps.txt (mEta mPt mBin),
# leaving acceptance.txt untouched. Validates against acceptance.txt before writing.
# Usage: python3 acc_fast.py SPS   |   python3 acc_fast.py DPS
import sys, os
import ROOT
from ROOT import TChain, TH2D
from array import array

mode = sys.argv[1].upper() if len(sys.argv) > 1 else "SPS"
CFG = {
    "SPS": dict(
        dir="/eos/home-l/leyao/26JJ/JPsiJPsi/GEN_nofilter/SPS/CMSSW_10_2_5/src/4mu_acc",
        base="/eos/user/c/chensh/JPsiJPsi/GEN_nofilter/SPS/CMSSW_10_2_5/src/4mu_acc/direct/SPS_2016_JJ_%d.root",
        require_other_fid=True),
    "DPS": dict(
        dir="/eos/home-l/leyao/26JJ/JPsiJPsi/GEN_nofilter/DPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc",
        base="/eos/user/c/chensh/JPsiJPsi/GEN_nofilter/DPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc/direct/DPS_2016_JJ_%d.root",
        require_other_fid=False),
}[mode]

ptBin = array('d', [i for i in range(10, 24)] + [24, 26, 28, 30, 35, 40])
yBin = array('d', [-2, -1.75, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 1.75, 2])
n, m = len(ptBin) - 1, len(yBin) - 1

ch = TChain("GenAnalyzer/gen_tree")
for i in range(1, 11):
    ch.Add(CFG["base"] % i)
print("[%s] entries=%d" % (mode, ch.GetEntries()))

# Exact expressions/cuts copied from efficiencyPlot.py
r1, r2 = "GEN_pair_id[0].first", "GEN_pair_id[0].second"
p1, y1 = "GENjpsi_pt[%s]" % r1, "GENjpsi_y[%s]" % r1
p2, y2 = "GENjpsi_pt[%s]" % r2, "GENjpsi_y[%s]" % r2
met1, mpt1 = "GENjpsi_mu_eta[%s]" % r1, "GENjpsi_mu_pt[%s]" % r1
met2, mpt2 = "GENjpsi_mu_eta[%s]" % r2, "GENjpsi_mu_pt[%s]" % r2
evt = "GENevt_mass[0] > 7.5"
fid1 = " && {p} > {a} && {p} < {b} && {y} > {c} && {y} < {d}".format(p=p1, a=ptBin[0], b=ptBin[-1], y=y1, c=yBin[0], d=yBin[-1])
fid2 = " && {p} > {a} && {p} < {b} && {y} > {c} && {y} < {d}".format(p=p2, a=ptBin[0], b=ptBin[-1], y=y2, c=yBin[0], d=yBin[-1])
eta1 = " && abs(%s[0]) < 2.4 && abs(%s[1]) < 2.4" % (met1, met1)
eta2 = " && abs(%s[0]) < 2.4 && abs(%s[1]) < 2.4" % (met2, met2)
pc1 = " && %s[0] > 3.5 && %s[1] > 3.5" % (mpt1, mpt1)
pc2 = " && %s[0] > 3.5 && %s[1] > 3.5" % (mpt2, mpt2)
oth1 = fid2 if CFG["require_other_fid"] else ""   # SPS requires the partner J/psi in fiducial
oth2 = fid1 if CFG["require_other_fid"] else ""


def fill(name, xexpr, yexpr, cut):
    # Iteration$==0 collapses the per-muon array instances so each EVENT is counted
    # once (TTree::Draw otherwise fills once per muon when the cut references the
    # 2-muon arrays, double-counting eta/pt-cut numerators). Matches GetEntries semantics.
    h = TH2D(name, "", n, ptBin, m, yBin)
    ch.Draw("%s:%s>>+%s" % (yexpr, xexpr, name), "(%s) && Iteration$==0" % cut, "goff")
    return h


# role1 (binned J/psi = first) and role2 (= second); sum the two roles per bin
hBin = fill("hBin1", p1, y1, evt + oth1)
hBin.Add(fill("hBin2", p2, y2, evt + oth2))
hEta = fill("hEta1", p1, y1, evt + oth1 + eta1)
hEta.Add(fill("hEta2", p2, y2, evt + oth2 + eta2))
hPt = fill("hPt1", p1, y1, evt + oth1 + eta1 + pc1)
hPt.Add(fill("hPt2", p2, y2, evt + oth2 + eta2 + pc2))

mBin = [[int(hBin.GetBinContent(j + 1, i + 1)) for j in range(n)] for i in range(m)]
mEta = [[int(hEta.GetBinContent(j + 1, i + 1)) for j in range(n)] for i in range(m)]
mPt = [[int(hPt.GetBinContent(j + 1, i + 1)) for j in range(n)] for i in range(m)]

# event-level header counts (2 scans)
nEvt = int(ch.GetEntries(evt + fid1 + fid2))
nAccEvt = int(ch.GetEntries(evt + fid1 + fid2 + eta1 + eta2 + pc1 + pc2))

# ---- validate mPt/mBin against existing acceptance.txt (per-bin: mPt mBin) ----
acc_txt = os.path.join(CFG["dir"], "plot", "acceptance.txt")
maxdiff = -1
if os.path.exists(acc_txt):
    L = open(acc_txt).read().split("\n")
    for i in range(m):
        vals = [int(x) for x in L[3 + i].split()]
        for j in range(n):
            ref_pt, ref_bin = vals[2 * j], vals[2 * j + 1]
            maxdiff = max(maxdiff, abs(ref_pt - mPt[i][j]), abs(ref_bin - mBin[i][j]))
    print("[%s] validation vs acceptance.txt: max|diff| in (mPt,mBin) = %d" % (mode, maxdiff))
else:
    print("[%s] no acceptance.txt to validate against" % mode)

# ---- write acceptance_maps.txt (mEta mPt mBin) ----
out = os.path.join(CFG["dir"], "plot", "acceptance_maps.txt")
with open(out, "w") as f:
    f.write("{} {} {} {}\n".format(n, m, nEvt, nAccEvt))
    f.write(" ".join(str(x) for x in ptBin) + " \n")
    f.write(" ".join(str(x) for x in yBin) + " \n")
    for i in range(m):
        row = ""
        for j in range(n):
            row += "%d %d %d " % (mEta[i][j], mPt[i][j], mBin[i][j])
        f.write(row + "\n")
print("[%s] wrote %s  (nEvt=%d nAccEvt=%d)" % (mode, out, nEvt, nAccEvt))
