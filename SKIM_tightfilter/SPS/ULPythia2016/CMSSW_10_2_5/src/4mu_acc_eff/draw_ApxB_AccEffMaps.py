# Redraw all AN Appendix B acc/eff 2D maps (and the 3 main-text AccEff maps) from the
# code-side count txt tables, in CMS tdrStyle. NO official-style generator exists in the
# accessible code area, so this reproduces the AN figures from the reliable local txt:
#   - acc:  GEN_*/4mu_acc/plot/acceptance_maps.txt  (per bin: mEta mPt mBin)   Pythia8 GEN
#   - eff:  SKIM_*/4mu_acc_eff/plot/raw_efficiency*.txt (Jpsi bin: mVtx mIdt mRec mBin;
#                                                        evt  bin: mTrg mHlt mVtx)
#           SPS eff = NLO* (raw_efficiency_NLO.txt, 96 files); DPS eff = raw_efficiency.txt (65)
# Two-step acceptance (matches AN): A_eta = mEta/mBin (mu |eta|<2.4); A_etapt = mPt/mEta (+ mu pt>3.5).
# Eff stages: eRECO=mRec/mBin, eID=mIdt/mRec, eMuMu=mVtx/mIdt, eHLT=mHlt/mVtx, e4mu=mTrg/mHlt.
# Uncertainties: binomial sqrt((N0-N1)*N1/N0)/N0.
# Output: AccEffMaps_out/{Acc_2016,Eff_2016,Effevt_2016}/<AN name>.png  (+ .pdf for the 3 main-text figs)
# Run: python3 draw_ApxB_AccEffMaps.py   (PyROOT)
import os, sys, ROOT
from ROOT import TH2D, TCanvas, TLatex, gStyle, gROOT
from array import array

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.append(HERE)
from tdrStyle import setTDRStyle

gROOT.SetBatch(True)
setTDRStyle()

REPO = "/eos/home-l/leyao/26JJ/JPsiJPsi"
ACC = {
    "SPS": REPO + "/GEN_nofilter/SPS/CMSSW_10_2_5/src/4mu_acc/plot/acceptance_maps.txt",
    "DPS": REPO + "/GEN_nofilter/DPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc/plot/acceptance_maps.txt",
}
EFF = {
    "SPS": REPO + "/SKIM_tightfilter/SPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc_eff/plot/raw_efficiency_NLO.txt",
    "DPS": REPO + "/SKIM_tightfilter/DPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc_eff/plot/raw_efficiency.txt",
}
OUT = os.path.join(HERE, "AccEffMaps_out")
for sub in ("Acc_2016", "Eff_2016", "Effevt_2016"):
    os.makedirs(os.path.join(OUT, sub), exist_ok=True)

# main-text Fig/AccEff figures that also need a .pdf copy
PDF_ALSO = {"acc2d_a_eta_DPS", "recoeff2d_a_DPS", "recoeff2d_evt_SPS"}


def binom_err(n1, n0):
    if n0 <= 0:
        return 0.0
    p = n1 / n0
    return (max(0.0, p * (1.0 - p)) / n0) ** 0.5   # clamp guards any n1>n0 rounding


def read_table(path, per_bin_jpsi, per_bin_evt=0):
    """Return (ptB, yB, n, m, jpsi_rows, evt_rows). jpsi_rows[i] is a flat list."""
    lines = open(path).read().split("\n")
    hdr = lines[0].split()
    n, m = int(hdr[0]), int(hdr[1])
    ptB = array('d', [float(x) for x in lines[1].split()])
    yB = array('d', [float(x) for x in lines[2].split()])
    jpsi = [[float(x) for x in lines[3 + i].split()] for i in range(m)]
    evt = []
    if per_bin_evt:
        evt = [[float(x) for x in lines[3 + m + i].split()] for i in range(n)]
    return ptB, yB, n, m, jpsi, evt


def make_map(n, m, ptB, yB, rows, ncol, inum, iden, square=False):
    """Build (eff, err) TH2D from rows; ncol columns/bin, numerator idx inum, denom idx iden."""
    yax = ptB if square else yB
    nrow = n if square else m
    eff = TH2D("e", "", n, ptB, nrow, yax)
    err = TH2D("d", "", n, ptB, nrow, yax)
    for i in range(nrow):
        vals = rows[i]
        for j in range(n):
            cell = vals[ncol * j: ncol * j + ncol]
            num, den = cell[inum], cell[iden]
            if den > 0:
                eff.SetBinContent(j + 1, i + 1, num / den)
                err.SetBinContent(j + 1, i + 1, binom_err(num, den))
    return eff, err


def draw(h, xtitle, ytitle, ztitle, fmt, name, sublabel, sub, zlo=None, zhi=None, square=False):
    w, ht = (1500, 1300) if square else (2200, 1000)
    c = TCanvas("c", "c", w, ht)
    c.SetTopMargin(0.10); c.SetRightMargin(0.17); c.SetLeftMargin(0.10); c.SetBottomMargin(0.12)
    h.GetXaxis().SetTitle(xtitle)
    h.GetYaxis().SetTitle(ytitle)
    h.GetZaxis().SetTitle(ztitle)
    h.GetXaxis().SetTitleOffset(1.1); h.GetYaxis().SetTitleOffset(0.9); h.GetZaxis().SetTitleOffset(0.6)
    h.GetXaxis().SetNdivisions(505); h.GetYaxis().SetNdivisions(505)
    if zlo is not None:
        h.GetZaxis().SetRangeUser(zlo, zhi)
    gStyle.SetPaintTextFormat(fmt)
    h.Draw("colzTEXT")
    lat = TLatex(); lat.SetNDC(); lat.SetTextFont(42); lat.SetTextSize(0.040)
    lat.DrawLatex(0.10, 0.925, "#bf{CMS} Preliminary")
    lat.DrawLatex(0.50, 0.925, sublabel)
    lat.DrawLatex(0.80, 0.925, "#sqrt{s} = 13 TeV")
    png = os.path.join(OUT, sub, name + "_JJto4mu.png")
    c.SaveAs(png)
    if name in PDF_ALSO:
        c.SaveAs(os.path.join(OUT, sub, name + "_JJto4mu.pdf"))
    del c


XJ, YJ = "p_{T}(J/#psi) [GeV]", "y(J/#psi)"
XE, YE = "p_{T}(J/#psi_{1}) [GeV]", "p_{T}(J/#psi_{2}) [GeV]"

for S in ("SPS", "DPS"):
    lbl = "%s_JJto4mu" % S
    # ---------- acceptance (Pythia8 GEN) : cols mEta mPt mBin ----------
    if os.path.exists(ACC[S]):
        ptB, yB, n, m, arows, _ = read_table(ACC[S], 3)
        eEta, dEta = make_map(n, m, ptB, yB, arows, 3, 0, 2)            # A_eta = mEta/mBin
        ePt, dPt = make_map(n, m, ptB, yB, arows, 3, 1, 0)             # A_etapt = mPt/mEta
        draw(eEta, XJ, YJ, "A_{#eta(#mu)}(J/#psi)", "1.2f", "acc2d_a_eta_%s" % S, lbl, "Acc_2016", 0., 1.)
        draw(dEta, XJ, YJ, "uncertainty", "1.3f", "dacc2d_a_eta_%s" % S, lbl, "Acc_2016")
        draw(ePt, XJ, YJ, "A_{p_{T}(#mu)}(J/#psi)", "1.2f", "acc2d_a_etapt_%s" % S, lbl, "Acc_2016", 0., 1.)
        draw(dPt, XJ, YJ, "uncertainty", "1.3f", "dacc2d_a_etapt_%s" % S, lbl, "Acc_2016")
    else:
        print("[skip acc %s] missing %s (acc job still running)" % (S, ACC[S]))
    # ---------- efficiency : Jpsi cols mVtx mIdt mRec mBin ; evt cols mTrg mHlt mVtx ----------
    ptB, yB, n, m, jrows, erows = read_table(EFF[S], 4, per_bin_evt=1)
    eRec, dRec = make_map(n, m, ptB, yB, jrows, 4, 2, 3)            # eRECO = mRec/mBin
    eId, dId = make_map(n, m, ptB, yB, jrows, 4, 1, 2)             # eID  = mIdt/mRec
    eVtx, dVtx = make_map(n, m, ptB, yB, jrows, 4, 0, 1)           # eMuMu= mVtx/mIdt
    draw(eRec, XJ, YJ, "#epsilon_{RECO(#mu)}(J/#psi)", "1.2f", "recoeff2d_a_%s" % S, lbl, "Eff_2016", 0., 1.)
    draw(dRec, XJ, YJ, "uncertainty", "1.3f", "drecoeff2d_a_%s" % S, lbl, "Eff_2016")
    draw(eId, XJ, YJ, "#epsilon_{ID(#mu)}(J/#psi)", "1.2f", "recoeff2d_id_a_%s" % S, lbl, "Eff_2016", 0., 1.)
    draw(dId, XJ, YJ, "uncertainty", "1.3f", "drecoeff2d_id_a_%s" % S, lbl, "Eff_2016")
    draw(eVtx, XJ, YJ, "#epsilon_{#mu#mu}(J/#psi)", "1.2f", "recoeff2d_id_vtx_a_%s" % S, lbl, "Eff_2016", 0., 1.)
    draw(dVtx, XJ, YJ, "uncertainty", "1.3f", "drecoeff2d_id_vtx_a_%s" % S, lbl, "Eff_2016")
    # evt maps (pt1 vs pt2)
    eHlt, dHlt = make_map(n, n, ptB, yB, erows, 3, 1, 2, square=True)   # eHLT = mHlt/mVtx
    eEvt, dEvt = make_map(n, n, ptB, yB, erows, 3, 0, 1, square=True)   # e4mu = mTrg/mHlt
    draw(eHlt, XE, YE, "#epsilon_{HLT}", "1.2f", "recoeff2d_trg_%s" % S, lbl, "Effevt_2016", 0., 1., square=True)
    draw(dHlt, XE, YE, "uncertainty", "1.3f", "drecoeff2d_trg_%s" % S, lbl, "Effevt_2016", square=True)
    draw(eEvt, XE, YE, "#epsilon_{4#mu}", "1.2f", "recoeff2d_evt_%s" % S, lbl, "Effevt_2016", 0., 1., square=True)
    draw(dEvt, XE, YE, "uncertainty", "1.3f", "drecoeff2d_evt_%s" % S, lbl, "Effevt_2016", square=True)

print("AccEffMaps done ->", OUT)
