#!/usr/bin/env python3
"""Run the accepted fixed-half method for the 0.85:0.15 correction generation."""

from pathlib import Path

import run_newdata_nominal_chain as chain


UPSTREAM_TAG = (
    "newdata167_accmix23_effmix19_sps0p85_dps0p15_lumi36p684_"
    "pairsym_cseed50_root640_strategy2_nominal_v4_20260902"
)
UPSTREAM = (
    chain.REPO
    / "Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/fit_results"
    / UPSTREAM_TAG
)

chain.UPSTREAM_TAG = UPSTREAM_TAG
chain.UPSTREAM = UPSTREAM
chain.MODEL_SOURCE = UPSTREAM / "total/Model_4D_tot.root"
chain.FIT_RESULT_SOURCE = (
    UPSTREAM / "total/Fit_4D_tot_native_asymptotic_unseeded.root"
)
chain.WEIGHT_DATA_SOURCE = (
    UPSTREAM
    / "inputs/WeightData_newdata167_accmix23_effmix19_sps0p85_dps0p15_"
    "pairsym_cseed50_v1_20260902.root"
)
chain.DATA_MANIFEST_SOURCE = (
    chain.REPO
    / "Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/"
    "data_ntuple_manifest_chensh_newdata_v1.list"
)
chain.ACCEPTANCE_SOURCE = (
    chain.REPO
    / "GEN_nofilter/DPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc/"
    "closure_results/nominal_mixed23_sps0p85_dps0p15_avgacc_v2_20260902/"
    "acceptance_sps0p85_dps0p15.txt"
)
chain.EFFICIENCY_SOURCE = (
    chain.REPO
    / "GEN_nofilter/DPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc/"
    "closure_results/inputs_mixed23_sps0p85_dps0p15_avgacc_v1_20260902/"
    "efficiency_sps0p85_dps0p15_avgacc_19x10.txt"
)
chain.SELECTION_SOURCE = UPSTREAM / "inputs/rephrase.cpp"
chain.EXPECTED = {
    chain.MODEL_SOURCE:
        "822b924c25b6f97121e6bf66adb708d41b4a2acbb08737a0718bdb03b4d00b7b",
    chain.FIT_RESULT_SOURCE:
        "27cd0530309dd8fb501586c9ecbf84319372ce5e3a4b47295cdd08635a54b0c6",
    chain.WEIGHT_DATA_SOURCE:
        "bc73bd9411c4e84ce343e183a68767377280e905cb042ee6bca61580af6ad7ff",
    chain.DATA_MANIFEST_SOURCE:
        "21c2830dc65e6be490f3ed262be30800b23cc92da61ea197305fed640fb84dee",
    chain.ACCEPTANCE_SOURCE:
        "8562012708b775d896b82f2d7d8e0bd8f424aa39d611019b98818105cfa6d52f",
    chain.EFFICIENCY_SOURCE:
        "fa21180118ad48257e9aa4be7f69b60052c431f36cad7d51696e9e57f53b9fc3",
    chain.SELECTION_SOURCE:
        "a4cd8bd802a452825025e0c70fb83ca55d8b2b13c77ffcc8887dab50e46a2541",
    UPSTREAM / "total/Fit_4D_tot.cpp":
        "9cccff5a8235a576a9fadc20dc936bdbb1ca1cae8f0dd1f60f7ae6916ab812c1",
    UPSTREAM / "total/total_fit_summary.csv":
        "37f415988f93ffcafde46553d791776d1b831456d18ee16d0da68a4315cbc2ee",
    UPSTREAM / "total/run.metadata.txt":
        "ef1520c2ecb9fc7851f277f14e907bf331f8f213a4fdfb552dcf82aea97811dc",
}
chain.WEIGHT_DATA_SNAPSHOT_NAME = (
    "WeightData_newdata167_accmix23_effmix19_sps0p85_dps0p15_"
    "pairsym_cseed50.root"
)
chain.ACCEPTANCE_SNAPSHOT_NAME = "acceptance_sps0p85_dps0p15.txt"
chain.EFFICIENCY_SNAPSHOT_NAME = (
    "efficiency_sps0p85_dps0p15_avgacc_19x10.txt"
)
chain.DPS_EXPECTED_FILES = 60
chain.DPS_EXCLUDED_INDICES = {11, 12, 13, 14, 15}
chain.DATA_LABEL_SYMMETRIZE = True
chain.DATA_LABEL_SEED = 50
chain.DPS_LABEL_SYMMETRIZE = False
chain.DPS_LABEL_SEED = 50
chain.EXPECTED_DATA_RANDOMIZATION = {
    "label_randomization_reseeded_files": 167,
    "randomization_draws": 5742,
    "randomized_candidate_swaps": 2784,
    "randomized_selected_event_swaps": 1403,
}
chain.CR_DEFINITIONS = [
    ("phi70", 1.8, chain.PHI70),
    ("phi110", 1.8, chain.PHI110),
    ("dy2p0", 2.0, chain.PHI90),
]
chain.EXPECTED_PP_FITS = 51
chain.RUNNER_RELATIVE = "Data_driven/run_sps0p85_dps0p15_nominal_chain.py"
chain.GENERATION_LABEL = "newdata167_sps0p85_dps0p15_pairsym_cseed50_dedup60"


if __name__ == "__main__":
    chain.main()
