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

`plot.cpp` reads one `Weight*.root` file and draws distributions into `plots/`.
The main switches are at the top of the file:

```cpp
const string inputFileName = "WeightSPS.root";
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

Run from `NtupleAnalyzer/`:

```bash
root -l -q -b plot.cpp
```

## Environment

The ROOT/RooFit macros are old-style analysis macros and are intended for the
lxplus EL7/CMSSW environment.
