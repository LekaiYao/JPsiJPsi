# 画 AN 附录B "private vs official MC 一致性" 中 Fig4.4 的 private 面板：
#   纯 DPS 的候选数(J/psi)分布，acceptance 筛选前 (Fig4.4leftup) 与筛选后 (Fig4.4leftdown)。
# 输入: chensh GEN-only DPS 样本 direct/DPS_2016_JJ_{1..10}.root（tree GenAnalyzer/gen_tree）。
#       cut 逻辑同本目录 efficiencyPlot.py：acc 前=evtMass>7.5 内的 (pt,y) 候选；acc 后再加 mu |eta|<2.4 & pT>3.5。
#       用 TTree::Draw 单遍填直方图（不要逐 bin GetEntries，否则极慢）。
# pt 轴: 0-40，与 official 同 region（展示低 pt 真实分布、便于事例数对比）；pt<5 不画（[0,5] 留空），与 official 一致。
# 输出: plot/Fig4.4leftup.pdf, plot/Fig4.4leftdown.pdf；拷到 AN/Fig/AccEff/ 替换 private 面板（official=right* 不动）。
# 样式与 official 对齐: CMS Preliminary + J/Psi(gen.) 轴 + colorbar "Number" 0-8000。
# 运行: python3 draw_ApxB_F4.4_DPS_Ndistri.py   (ROOT/PyROOT 环境)
import sys, os, ROOT
from ROOT import TChain, TH2D, TCanvas, TLatex, gStyle, gROOT
from array import array

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.append(HERE)
from tdrStyle import setTDRStyle

gROOT.SetBatch(True)
setTDRStyle()

# pt 轴 0-40：第一个 bin [0,5] 留空（pt>5 才填），其余沿用 efficiency binning 的细分
ptB = array('d', [0, 5] + [i for i in range(6, 24)] + [24, 26, 28, 30, 35, 40])
yB  = array('d', [-2, -1.75, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 1.75, 2])
n, m = len(ptB) - 1, len(yB) - 1
ZMAX = 8000.   # colorbar 上限，与 official 一致

ch = TChain("GenAnalyzer/gen_tree")
base = "/eos/user/c/chensh/JPsiJPsi/GEN_nofilter/DPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc/direct/DPS_2016_JJ_%d.root"
for i in range(1, 11):
    ch.Add(base % i)

rJ1, rJ2 = "GEN_pair_id[0].first", "GEN_pair_id[0].second"
rJ1_pt, rJ1_y = "GENjpsi_pt[%s]" % rJ1, "GENjpsi_y[%s]" % rJ1
rJ2_pt, rJ2_y = "GENjpsi_pt[%s]" % rJ2, "GENjpsi_y[%s]" % rJ2
rJ1_meta, rJ1_mpt = "GENjpsi_mu_eta[%s]" % rJ1, "GENjpsi_mu_pt[%s]" % rJ1
rJ2_meta, rJ2_mpt = "GENjpsi_mu_eta[%s]" % rJ2, "GENjpsi_mu_pt[%s]" % rJ2
em = "GENevt_mass[0] > 7.5"
g1 = " && %s > 5" % rJ1_pt   # pt>5，保持 [0,5] 空、与 official 同 region
g2 = " && %s > 5" % rJ2_pt
eta1 = " && abs(%s[0]) < 2.4 && abs(%s[1]) < 2.4" % (rJ1_meta, rJ1_meta)
eta2 = " && abs(%s[0]) < 2.4 && abs(%s[1]) < 2.4" % (rJ2_meta, rJ2_meta)
pt1  = " && %s[0] > 3.5 && %s[1] > 3.5" % (rJ1_mpt, rJ1_mpt)
pt2  = " && %s[0] > 3.5 && %s[1] > 3.5" % (rJ2_mpt, rJ2_mpt)

hBin = TH2D("hBin", "", n, ptB, m, yB)   # before acceptance
hAcc = TH2D("hAcc", "", n, ptB, m, yB)   # after acceptance
ch.Draw("%s:%s>>+hBin" % (rJ1_y, rJ1_pt), em + g1, "goff")
ch.Draw("%s:%s>>+hBin" % (rJ2_y, rJ2_pt), em + g2, "goff")
ch.Draw("%s:%s>>+hAcc" % (rJ1_y, rJ1_pt), em + g1 + eta1 + pt1, "goff")
ch.Draw("%s:%s>>+hAcc" % (rJ2_y, rJ2_pt), em + g2 + eta2 + pt2, "goff")
print("filled: hBin max=%d ; hAcc max=%d" % (hBin.GetMaximum(), hAcc.GetMaximum()))

def draw(h, outname):
    c = TCanvas("c", "c", 2000, 900)
    c.SetTopMargin(0.10); c.SetRightMargin(0.17); c.SetLeftMargin(0.10); c.SetBottomMargin(0.12)
    h.GetXaxis().SetTitle("p_{T} (J/#Psi)(gen.)")
    h.GetYaxis().SetTitle("rapidity (J/#Psi)(gen.)")
    h.GetZaxis().SetTitle("Number")
    h.GetXaxis().SetTitleOffset(1.1); h.GetYaxis().SetTitleOffset(0.8); h.GetZaxis().SetTitleOffset(0.6)
    h.GetXaxis().SetNdivisions(505); h.GetYaxis().SetNdivisions(505)
    h.GetZaxis().SetRangeUser(0., ZMAX)
    gStyle.SetPaintTextFormat("1.0f")
    h.Draw("colzTEXT")
    lat = TLatex(); lat.SetNDC(); lat.SetTextFont(42); lat.SetTextSize(0.040)
    lat.DrawLatex(0.10, 0.925, "#bf{CMS} Preliminary")
    lat.DrawLatex(0.52, 0.925, "DPS_JJto4mu")
    lat.DrawLatex(0.80, 0.925, "#sqrt{s} = 13 TeV")
    c.SaveAs(outname)
    del c

outdir = os.path.join(HERE, "plot")
draw(hBin, os.path.join(outdir, "Fig4.4leftup.pdf"))    # private, before acceptance
draw(hAcc, os.path.join(outdir, "Fig4.4leftdown.pdf"))  # private, after acceptance
print("F4.4 done ->", outdir)
