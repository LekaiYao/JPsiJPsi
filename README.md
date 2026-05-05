# JPsiJPsi

J/psi pair four-muon analysis code. The workflow covers GEN-level acceptance,
reconstruction efficiency, weighted data/MC ntuple production, 4D fitting,
differential fits, and summary plots.

## Main Workflow

1. Compute GEN acceptance in `4mu_acc/` and produce `plot/acceptance.txt`.
2. Compute or combine reconstruction efficiency in `4mu_acc_eff/`.
3. Copy `acceptance.txt` and `efficiency_2SPS+DPS.txt` into `NtupleAnalyzer/`.
4. Run `rephrase.cpp` to produce `WeightData.root`, `WeightSPS.root`, or `WeightDPS.root`.
5. Run `Fit_4D_tot.cpp` for the inclusive 4D fit.
6. Run `make all` for differential-bin fits.
7. Run `make Tpl` for differential cross-section summary plots.

## NtupleAnalyzer Notes

`rephrase.cpp` applies the reco-level event selection, computes kinematic
variables, applies acceptance/efficiency weights, and writes a flat tree named
`data` in `Weight*.root`.

`GEN_rephrase.cpp` is a preserved variant of `rephrase.cpp` that also reads
`GEJpsi1_phi` and `GEJpsi2_phi` and writes `delta_phi_GEN` to the output tree.
Use this variant when downstream plots need GEN-level delta-phi.

`mix.cpp` randomly mixes LO and NLO ntuple events while preserving the original
`rootuple/oniaTree` format. Its top-level settings define the LO/NLO input
paths, `n_NLO`, the LO:NLO ratio `r`, the random seed, and the output file.

Run the mixer from `NtupleAnalyzer/`:

```bash
root -l -b -q mix.cpp
```

Do not run it as `.L mix.cpp+` unless compiled ROOT artifacts are intended.

`plot.cpp` reads one `Weight*.root` file and draws distributions into `plots/`.
The main switches are at the top of the file:

```cpp
const string inputFileName = "WeightSPS_GEN.root";
const bool drawWeightComparison = false;
```

By default, `plot.cpp` writes only weighted distributions. Set
`drawWeightComparison = true` to draw unweighted and weighted histograms side by
side. Current outputs include:

- `plots/JpsiMass1.png`
- `plots/JpsiMass2.png`
- `plots/JpsiCtau1.png`
- `plots/JpsiCtau2.png`
- `plots/evtMass.png`
- `plots/deltaY.png`
- `plots/deltaPhi.png`
- `plots/deltaPhiGEN.png` when the input tree contains `delta_phi_GEN`

Run from `NtupleAnalyzer/`:

```bash
root -l -q -b plot.cpp
```

## Environment

The ROOT/RooFit macros are old-style analysis macros and are intended for the
lxplus EL7/CMSSW environment.
