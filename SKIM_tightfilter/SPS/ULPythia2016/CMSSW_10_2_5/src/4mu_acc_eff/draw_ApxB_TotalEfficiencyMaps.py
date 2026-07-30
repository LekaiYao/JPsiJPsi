#!/usr/bin/env python3
"""Build standalone factorized event total-efficiency preview maps."""
import argparse, csv, math, os, subprocess, sys
from array import array
from collections import defaultdict
from pathlib import Path

import ROOT

ROOT.gROOT.SetBatch(True)
REPO = Path("/eos/home-l/leyao/26JJ/JPsiJPsi")
HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description="Draw factorized event total-efficiency maps")
parser.add_argument("--output-dir", default=str(HERE/"AccEffMaps_out"))
args = parser.parse_args()
OUT = Path(args.output_dir)/"Effevt_2016"
PT = [float(i) for i in range(10, 24)] + [24., 26., 28., 30., 35., 40.]
Y = [-2., -1.75, -1.5, -1., -0.5, 0., 0.5, 1., 1.5, 1.75, 2.]

CONFIG = {
    "SPS": {
        "acc": REPO/"GEN_nofilter/SPS/CMSSW_10_2_5/src/4mu_acc/plot/acceptance_maps.txt",
        "eff": REPO/"SKIM_tightfilter/SPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc_eff/plot/raw_efficiency_NLO.txt",
        "files": [f"/eos/home-l/leyao/26JJ/MC_Maker/HelacOnia2016/CMSSW_10_6_20/src/NTUPLE/NLO_gpt0p8/Ntuple_2016_SPSstar_{i}.root" for i in range(1,97)],
        "label": "SPS",
    },
    "DPS": {
        "acc": REPO/"GEN_nofilter/DPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc/plot/acceptance_maps.txt",
        "eff": REPO/"SKIM_tightfilter/DPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc_eff/plot/raw_efficiency.txt",
        "files": [str(REPO/f"Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/DPS_ntuple/Ntuple_2016_DPS_{i}.root") for i in range(1,66)],
        "label": "DPS",
    },
}

def berr(num, den):
    if den <= 0: return 0.
    p=num/den
    return math.sqrt(max(0.,p*(1-p)/den))

def read(sample):
    a=open(CONFIG[sample]["acc"]).read().splitlines()
    e=open(CONFIG[sample]["eff"]).read().splitlines()
    n, m=map(int,a[0].split()[:2])
    assert n==len(PT)-1 and m==len(Y)-1
    acc=[]; eff=[]
    for iy in range(m):
        ar=list(map(float,a[3+iy].split())); er=list(map(float,e[3+iy].split()))
        acc.append([]); eff.append([])
        for ip in range(n):
            meta,mpt,mbin=ar[3*ip:3*ip+3]
            vtx,ident,reco,den=er[4*ip:4*ip+4]
            acc[-1].append((mpt/mbin if mbin else 0.,berr(mpt,mbin)))
            eff[-1].append((vtx/den if den else 0.,berr(vtx,den)))
    evt=[]
    for iy in range(n):
        row=list(map(float,e[3+m+iy].split())); evt.append([])
        for ix in range(n):
            trg,hlt,vtx=row[3*ix:3*ix+3]
            evt[-1].append((trg/vtx if vtx else 0.,berr(trg,vtx)))
    return acc,eff,evt

def ibin(x, edges):
    for i in range(len(edges)-1):
        if x>edges[i] and x<edges[i+1]: return i
    return None

def make(sample):
    acc,eff,eeff=read(sample)
    chain=ROOT.TChain("rootuple/oniaTree")
    for f in CONFIG[sample]["files"]: chain.Add(f)
    sums=defaultdict(float); counts=defaultdict(int)
    deriv=defaultdict(lambda:defaultdict(float))
    selected=0
    for t in chain:
        if not bool(t.GEevt_valid) or t.GEevt_fourMuMass<=7.5: continue
        p1,p2=float(t.GEJpsi1_pt),float(t.GEJpsi2_pt)
        y1,y2=float(t.GEJpsi1_y),float(t.GEJpsi2_y)
        i1,i2,j1,j2=ibin(p1,PT),ibin(p2,PT),ibin(y1,Y),ibin(y2,Y)
        if None in (i1,i2,j1,j2): continue
        A1,A2=acc[j1][i1][0],acc[j2][i2][0]
        E1,E2=eff[j1][i1][0],eff[j2][i2][0]
        H=eeff[i2][i1][0]
        if min(A1,A2,E1,E2,H)<=0: continue
        f=A1*A2*E1*E2*H; out=(i1,i2)
        sums[out]+=f; counts[out]+=1; selected+=1
        for key,val in ((("A",j1,i1),A1),(("A",j2,i2),A2),
                        (("E",j1,i1),E1),(("E",j2,i2),E2),(("H",i2,i1),H)):
            deriv[out][key]+=f/val
    vals={}; errs={}
    for out,nc in counts.items():
        vals[out]=sums[out]/nc; var=0.
        for key,d in deriv[out].items():
            typ=key[0]
            sig=acc[key[1]][key[2]][1] if typ=="A" else eff[key[1]][key[2]][1] if typ=="E" else eeff[key[1]][key[2]][1]
            var+=(d/nc*sig)**2
        errs[out]=math.sqrt(var)
    return vals,errs,counts,chain.GetEntries(),selected

def draw(sample, values, errors):
    ROOT.gStyle.SetOptStat(0); ROOT.gStyle.SetPaintTextFormat("1.3f")
    for kind,data,ztitle,name in [
        ("efficiency",values,"Total efficiency",f"recoeff2d_total_{sample}_JJto4mu"),
        ("uncertainty",errors,"Absolute statistical uncertainty",f"drecoeff2d_total_{sample}_JJto4mu")]:
        h=ROOT.TH2D("h_"+name,"",len(PT)-1,array('d',PT),len(PT)-1,array('d',PT))
        for (ix,iy),v in data.items(): h.SetBinContent(ix+1,iy+1,v)
        c=ROOT.TCanvas("c_"+name,"",2400,1800)
        c.SetTopMargin(.11); c.SetRightMargin(.18); c.SetLeftMargin(.12); c.SetBottomMargin(.16)
        h.GetXaxis().SetTitle("p_{T}(J/#psi_{1}) [GeV]")
        h.GetYaxis().SetTitle("p_{T}(J/#psi_{2}) [GeV]"); h.GetZaxis().SetTitle(ztitle)
        h.GetXaxis().SetTitleOffset(1.2); h.GetYaxis().SetTitleOffset(1.1); h.GetZaxis().SetTitleOffset(1.05)
        for ax in (h.GetXaxis(),h.GetYaxis(),h.GetZaxis()): ax.SetTitleSize(.05); ax.SetLabelSize(.04)
        h.GetXaxis().SetNdivisions(505); h.GetYaxis().SetNdivisions(505)
        if kind=="efficiency": h.GetZaxis().SetRangeUser(0.,1.)
        h.Draw("COLZ")
        labels=[]
        for ix in range(1,h.GetNbinsX()+1):
            for iy in range(1,h.GetNbinsY()+1):
                v=h.GetBinContent(ix,iy)
                if not v: continue
                s=f"{v:.2f}" if kind=="efficiency" else f"{v:.3f}"
                q=ROOT.TLatex(h.GetXaxis().GetBinCenter(ix),h.GetYaxis().GetBinCenter(iy),s)
                q.SetTextAlign(22); q.SetTextFont(62); q.SetTextSize(.012); q.Draw(); labels.append(q)
        cms=ROOT.TLatex(); cms.SetNDC(); cms.SetTextFont(61); cms.SetTextSize(.04); cms.DrawLatex(.12,.925,"CMS")
        ext=ROOT.TLatex(); ext.SetNDC(); ext.SetTextFont(52); ext.SetTextSize(.034); ext.DrawLatex(.195,.925,"Simulation Preliminary")
        lab=ROOT.TLatex(); lab.SetNDC(); lab.SetTextFont(42); lab.SetTextAlign(22); lab.SetTextSize(.030); lab.DrawLatex(.57,.925,CONFIG[sample]["label"])
        ene=ROOT.TLatex(); ene.SetNDC(); ene.SetTextFont(42); ene.SetTextAlign(31); ene.SetTextSize(.04); ene.DrawLatex(.82,.925,"13 TeV")
        eps=OUT/f"{name}.eps"
        c.SaveAs(str(eps))
        subprocess.check_call(["ps2pdf","-dEPSCrop",str(eps),str(OUT/f"{name}.pdf")])
        eps.unlink()

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    summary=[]
    for s in ("SPS","DPS"):
        v,e,n,alln,sel=make(s); draw(s,v,e)
        non=list(v.values())
        summary.append(f"{s}: chain entries={alln}, projected fiducial events={sel}, populated bins={len(non)}, min={min(non):.6g}, max={max(non):.6g}")
    print("\n".join(summary))

if __name__=="__main__": main()
