# 画 AN 附录B "private vs official MC 一致性" 中 Fig4.5 的 private 面板：
#   纯 SPS 的 mu soft-ID 效率 map (Fig4.5leftup) 与其不确定度 (Fig4.5leftdown)。
# 输入: 本目录 plot/raw_efficiency_NLO.txt（96 个合并 SPS NLO* 样本，每 bin 存 mVtx mIdt mRec mBin）。
#       ε_ID = mIdt/mRec，不确定度按二项误差。
# 输出: plot/Fig4.5leftup.pdf, plot/Fig4.5leftdown.pdf；拷到 AN/Fig/AccEff/ 替换 private 面板（official=right* 不动）。
# 样式与 official 对齐: CMS Preliminary + J/Psi 轴 + colorbar(map) 0.45-1。
# 运行: python3 draw_ApxB_F4.5_SPS_effID.py   (ROOT/PyROOT 环境)
import sys, os, ROOT
from ROOT import TH2D, TCanvas, TLatex, gStyle, gROOT
from array import array

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.append(HERE)
from tdrStyle import setTDRStyle

gROOT.SetBatch(True)
setTDRStyle()

raw = os.path.join(HERE, "plot", "raw_efficiency_NLO.txt")
lines = open(raw).read().split("\n")
hdr = lines[0].split()
n, m = int(hdr[0]), int(hdr[1])
ptB = array('d', [float(x) for x in lines[1].split()])
yB  = array('d', [float(x) for x in lines[2].split()])

eff = TH2D("eff", "", n, ptB, m, yB)
err = TH2D("err", "", n, ptB, m, yB)
for i in range(m):
    vals = [float(x) for x in lines[3 + i].split()]
    for j in range(n):
        mVtx, mIdt, mRec, mBin = vals[4*j:4*j+4]
        if mRec > 0:
            e = mIdt / mRec
            eff.SetBinContent(j+1, i+1, e)
            err.SetBinContent(j+1, i+1, ((mRec - mIdt) * mIdt / mRec) ** 0.5 / mRec)

def draw(h, ztitle, fmt, outname, zlo=None, zhi=None):
    c = TCanvas("c", "c", 2000, 900)
    c.SetTopMargin(0.10); c.SetRightMargin(0.17); c.SetLeftMargin(0.10); c.SetBottomMargin(0.12)
    h.GetXaxis().SetTitle("p_{T} (J/#Psi)(reco.)")
    h.GetYaxis().SetTitle("rapidity (J/#Psi)(reco.)")
    h.GetZaxis().SetTitle(ztitle)
    h.GetXaxis().SetTitleOffset(1.1); h.GetYaxis().SetTitleOffset(0.8); h.GetZaxis().SetTitleOffset(0.6)
    h.GetXaxis().SetNdivisions(505); h.GetYaxis().SetNdivisions(505)
    if zlo is not None: h.GetZaxis().SetRangeUser(zlo, zhi)
    gStyle.SetPaintTextFormat(fmt)
    h.Draw("colzTEXT")
    lat = TLatex(); lat.SetNDC(); lat.SetTextFont(42); lat.SetTextSize(0.040)
    lat.DrawLatex(0.10, 0.925, "#bf{CMS} Preliminary")
    lat.DrawLatex(0.52, 0.925, "SPS_JJto4mu")
    lat.DrawLatex(0.80, 0.925, "#sqrt{s} = 13 TeV")
    c.SaveAs(outname)
    del c

outdir = os.path.join(HERE, "plot")
draw(eff, "#mu soft id. eff.", "1.2f", os.path.join(outdir, "Fig4.5leftup.pdf"), 0.45, 1.0)
draw(err, "uncertainty",        "1.3f", os.path.join(outdir, "Fig4.5leftdown.pdf"))
print("F4.5 done ->", outdir)
