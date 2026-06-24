# Validation — Appendix C MC-vs-data kinematics comparison

复现并更新 AN Appendix C §"Comparison between MC and data" 的 7 张 MC-vs-data 对比图
（sideband-subtracted data vs SPS MC vs DPS MC）。本子功能从
`work/JJ/CMSSW_12_4_0/src/Differential_cross_section/Validation` 移植而来，
**唯一改动是把 SPS MC 从旧的 Pythia8(LO) 换成 nominal 的 NLO\*（HELAC-Onia）**。
DPS MC 和 sWeighted data 沿用原 pipeline 不变（DPS generator 未变）。

## 文件

- `make_validation_mc.cpp` — selector：读原始 NLO\* ntuple（`rootuple/oniaTree`），施加与
  `../rephrase.cpp` 完全相同的选择（passHLT、matchTrg、两个 J/psi pT∈[10,40]、samePV、
  M(JJ)≥7.5，倒序取第一个通过的 candidate），查 `../acceptance.txt` 与
  `../efficiency_0_0.6.txt` 算 `Weight_sum`，输出 `WeightedTree` rich 格式
  （J/psi-pair 级 + 单 J reco + 4 个 muon reco 分支）。输出 `MC_SPS_NLOstar.root`。
- `plot_validation.C` — SBS 对比图模板（占位符 `@DATAVAR@/@MCVAR@/@TITLE@/@NBINS@/@XMIN@/@XMAX@`）。
  MC 用 reco 级 shape（不加权，归一），data 用 sideband-subtraction（prompt ctau-cut）。
- `run_validation.sh` — 对 C3 的 7 个变量（DeltaPhi, DeltaY, FourMuonMass, J1Pt, J1y,
  Mu1Pt, Mu1Eta）替换占位符并出图到 `output_sbs/`。
- `Sample/MC_DPS.root` — 原 pipeline 的 DPS MC（未变，拷入）。
- `Sample/Data_MixWeighted_cut0_sw.root` — 原 pipeline 的 sWeighted data（未变，拷入）。
- `output_sbs/*_compare_SBS_MC.pdf` — 输出图。

## 运行

```bash
cd Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/Validation
root -l -q -b make_validation_mc.cpp   # 生成 MC_SPS_NLOstar.root
bash run_validation.sh                 # 生成 output_sbs/ 下 7 张对比图
```

## 与 AN 图的对应

| output_sbs 文件 | AN Fig/Kinematics |
|---|---|
| DeltaPhi_compare_SBS_MC.pdf | DeltaPhi_MCvsData |
| DeltaY_compare_SBS_MC.pdf | DeltaY_MCvsData |
| FourMuonMass_compare_SBS_MC.pdf | Minv_MCvsData |
| J1Pt_compare_SBS_MC.pdf | J1Pt_MCvsData |
| J1y_compare_SBS_MC.pdf | J1Y_MCvsData |
| Mu1Pt_compare_SBS_MC.pdf | Mu1Pt_MCvsData |
| Mu1Eta_compare_SBS_MC.pdf | Mu1Eta_MCvsData |

## 注意

- 切换到 NLO\* 后，SPS MC 在 `DeltaPhi` 末 bin（back-to-back）显著抬升（~18%），
  旧 Pythia8 几乎为 0。这与 Appendix C §"Back-to-back" 的更新一致。
- selector 默认只处理 NLO\* 的 `_1.root`（与 `rephrase.cpp` 当前 nominal 一致，
  ~15k 选出事例，对 shape 比较足够）；需要更多统计可在文件顶部加大 `suffix`。
- MC 对比图只比 shape，不加 `Weight_sum`（acc/eff 权重只影响截面归一，不影响 reco shape）。
