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
# Output: AccEffMaps_out/{Acc_2016,Eff_2016,Effevt_2016}/<AN name>.pdf
# Run: python3 draw_ApxB_AccEffMaps.py   (PyROOT)
import argparse, os, subprocess, sys, ROOT
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
parser = argparse.ArgumentParser(description="Draw the canonical acceptance and efficiency maps")
parser.add_argument("--output-dir", default=os.path.join(HERE, "AccEffMaps_out"))
args = parser.parse_args()
OUT = os.path.abspath(args.output_dir)
for sub in ("Acc_2016", "Eff_2016", "Effevt_2016"):
    os.makedirs(os.path.join(OUT, sub), exist_ok=True)

LAYOUTS = {
    "pt_y": {
        "canvas": (2400, 1200), "top": 0.11, "right": 0.18, "left": 0.11, "bottom": 0.15,
        "x_offset": 1.15, "y_offset": 1.00, "z_offset": 1.05, "bin_text_size": 0.018,
    },
    "pt_pt": {
        "canvas": (2400, 1800), "top": 0.11, "right": 0.18, "left": 0.12, "bottom": 0.16,
        "x_offset": 1.20, "y_offset": 1.10, "z_offset": 1.05, "bin_text_size": 0.012,
    },
}

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


def draw_bin_text(h, layout, value_kind):
    precision = 2 if value_kind == "efficiency" else 3
    labels = []
    for ix in range(1, h.GetNbinsX() + 1):
        for iy in range(1, h.GetNbinsY() + 1):
            value = h.GetBinContent(ix, iy)
            if value == 0.0:
                continue
            label = TLatex(h.GetXaxis().GetBinCenter(ix), h.GetYaxis().GetBinCenter(iy), ("%.*f" % (precision, value)))
            label.SetTextAlign(22); label.SetTextFont(62); label.SetTextSize(LAYOUTS[layout]["bin_text_size"])
            label.Draw(); labels.append(label)
    return labels


def draw_map(h, xtitle, ytitle, ztitle, name, sub, layout, value_kind, sample_label, zlo=None, zhi=None):
    cfg = LAYOUTS[layout]
    c = TCanvas("c_" + name, "", cfg["canvas"][0], cfg["canvas"][1])
    c.SetTopMargin(cfg["top"]); c.SetRightMargin(cfg["right"]); c.SetLeftMargin(cfg["left"]); c.SetBottomMargin(cfg["bottom"])
    h.GetXaxis().SetTitle(xtitle); h.GetYaxis().SetTitle(ytitle); h.GetZaxis().SetTitle(ztitle)
    h.GetXaxis().SetTitleOffset(cfg["x_offset"]); h.GetYaxis().SetTitleOffset(cfg["y_offset"]); h.GetZaxis().SetTitleOffset(cfg["z_offset"])
    for axis in (h.GetXaxis(), h.GetYaxis(), h.GetZaxis()):
        axis.SetTitleSize(0.050); axis.SetLabelSize(0.040)
    h.GetXaxis().SetNdivisions(505); h.GetYaxis().SetNdivisions(505)
    if zlo is not None:
        h.GetZaxis().SetRangeUser(zlo, zhi)
    h.Draw("COLZ")
    bin_labels = draw_bin_text(h, layout, value_kind)
    cms = TLatex(); cms.SetNDC(); cms.SetTextFont(61); cms.SetTextSize(0.040)
    cms.DrawLatex(cfg["left"], 0.925, "CMS")
    extra = TLatex(); extra.SetNDC(); extra.SetTextFont(52); extra.SetTextSize(0.034)
    extra.DrawLatex(cfg["left"] + 0.075, 0.925, "Simulation Preliminary")
    sample = TLatex(); sample.SetNDC(); sample.SetTextFont(42); sample.SetTextAlign(22); sample.SetTextSize(0.034)
    sample.DrawLatex(0.56, 0.925, sample_label)
    energy = TLatex(); energy.SetNDC(); energy.SetTextFont(42); energy.SetTextAlign(31); energy.SetTextSize(0.040)
    energy.DrawLatex(1.0 - cfg["right"], 0.925, "13 TeV")
    output_base = os.path.join(OUT, sub, name + "_JJto4mu")
    eps_path = output_base + ".eps"
    c.SaveAs(eps_path)
    subprocess.check_call(["ps2pdf", "-dEPSCrop", eps_path, output_base + ".pdf"])
    os.remove(eps_path)
    del c


XJ, YJ = "p_{T}(J/#psi) [GeV]", "y(J/#psi)"
XE, YE = "p_{T}(J/#psi_{1}) [GeV]", "p_{T}(J/#psi_{2}) [GeV]"

for S in ("SPS", "DPS"):
    acceptance_label = "Pythia8 %s" % S
    efficiency_label = "HELAC-Onia NLO* SPS" if S == "SPS" else "Pythia8 DPS"
    if os.path.exists(ACC[S]):
        ptB, yB, n, m, arows, _ = read_table(ACC[S], 3)
        eEta, dEta = make_map(n, m, ptB, yB, arows, 3, 0, 2)
        ePt, dPt = make_map(n, m, ptB, yB, arows, 3, 1, 0)
        draw_map(eEta, XJ, YJ, "A_{#eta(#mu)}(J/#psi)", "acc2d_a_eta_%s" % S, "Acc_2016", "pt_y", "efficiency", acceptance_label, 0., 1.)
        draw_map(dEta, XJ, YJ, "Uncertainty", "dacc2d_a_eta_%s" % S, "Acc_2016", "pt_y", "uncertainty", acceptance_label)
        draw_map(ePt, XJ, YJ, "A_{p_{T}(#mu)}(J/#psi)", "acc2d_a_etapt_%s" % S, "Acc_2016", "pt_y", "efficiency", acceptance_label, 0., 1.)
        draw_map(dPt, XJ, YJ, "Uncertainty", "dacc2d_a_etapt_%s" % S, "Acc_2016", "pt_y", "uncertainty", acceptance_label)
    else:
        print("[skip acc %s] missing %s" % (S, ACC[S]))

    ptB, yB, n, m, jrows, erows = read_table(EFF[S], 4, per_bin_evt=1)
    eRec, dRec = make_map(n, m, ptB, yB, jrows, 4, 2, 3)
    eId, dId = make_map(n, m, ptB, yB, jrows, 4, 1, 2)
    eVtx, dVtx = make_map(n, m, ptB, yB, jrows, 4, 0, 1)
    draw_map(eRec, XJ, YJ, "#epsilon_{RECO(#mu)}(J/#psi)", "recoeff2d_a_%s" % S, "Eff_2016", "pt_y", "efficiency", efficiency_label, 0., 1.)
    draw_map(dRec, XJ, YJ, "Uncertainty", "drecoeff2d_a_%s" % S, "Eff_2016", "pt_y", "uncertainty", efficiency_label)
    draw_map(eId, XJ, YJ, "#epsilon_{ID(#mu)}(J/#psi)", "recoeff2d_id_a_%s" % S, "Eff_2016", "pt_y", "efficiency", efficiency_label, 0., 1.)
    draw_map(dId, XJ, YJ, "Uncertainty", "drecoeff2d_id_a_%s" % S, "Eff_2016", "pt_y", "uncertainty", efficiency_label)
    draw_map(eVtx, XJ, YJ, "#epsilon_{#mu#mu}(J/#psi)", "recoeff2d_id_vtx_a_%s" % S, "Eff_2016", "pt_y", "efficiency", efficiency_label, 0., 1.)
    draw_map(dVtx, XJ, YJ, "Uncertainty", "drecoeff2d_id_vtx_a_%s" % S, "Eff_2016", "pt_y", "uncertainty", efficiency_label)

    eHlt, dHlt = make_map(n, n, ptB, yB, erows, 3, 1, 2, square=True)
    eEvt, dEvt = make_map(n, n, ptB, yB, erows, 3, 0, 1, square=True)
    draw_map(eHlt, XE, YE, "#epsilon_{HLT}", "recoeff2d_trg_%s" % S, "Effevt_2016", "pt_pt", "efficiency", efficiency_label, 0., 1.)
    draw_map(dHlt, XE, YE, "Uncertainty", "drecoeff2d_trg_%s" % S, "Effevt_2016", "pt_pt", "uncertainty", efficiency_label)
    draw_map(eEvt, XE, YE, "#epsilon_{evt}", "recoeff2d_evt_%s" % S, "Effevt_2016", "pt_pt", "efficiency", efficiency_label, 0., 1.)
    draw_map(dEvt, XE, YE, "Uncertainty", "drecoeff2d_evt_%s" % S, "Effevt_2016", "pt_pt", "uncertainty", efficiency_label)

subprocess.check_call([sys.executable, os.path.join(HERE, "draw_ApxB_TotalEfficiencyMaps.py"), "--output-dir", OUT])

print("AccEffMaps done ->", OUT)
