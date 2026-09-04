#!/usr/bin/env python3
"""Run the 0.8:0.2 correction generation with the deduplicated DPS inventory."""

import run_newdata_nominal_chain as chain


UPSTREAM_TAG = "newdata167_accmix23_effmix19_lumi36p684_root640_v1_20260902"
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
    UPSTREAM / "inputs/WeightData_newdata167_accmix23_effmix19_v1_20260902.root"
)
chain.DATA_MANIFEST_SOURCE = (
    chain.REPO
    / "Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/"
    "data_ntuple_manifest_chensh_newdata_v1.list"
)
chain.ACCEPTANCE_SOURCE = (
    chain.REPO
    / "GEN_nofilter/DPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc/"
    "closure_results/nominal_mixed23_sps0p8_dps0p2_avgacc_v1_20260902/"
    "acceptance_sps0p8_dps0p2.txt"
)
chain.EFFICIENCY_SOURCE = (
    chain.REPO
    / "GEN_nofilter/DPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc/"
    "closure_results/inputs_mixed23_sps0p8_dps0p2_avgacc_v1_20260902/"
    "efficiency_sps0p8_dps0p2_avgacc_19x10.txt"
)
chain.SELECTION_SOURCE = UPSTREAM / "inputs/rephrase.cpp"
chain.EXPECTED = {
    chain.MODEL_SOURCE:
        "9987bf734f9308f4177b290828cff82c7cb470d035cca8bbdba337534af61238",
    chain.FIT_RESULT_SOURCE:
        "d0a237af9042d663cce10a217b5a869963bf6c697887a5ef6acabfb93f3cf25a",
    chain.WEIGHT_DATA_SOURCE:
        "8cd1abae2a896f5ebf82a805e62009e96a7f7e4102350c1f6e18c2568c1e8957",
    chain.DATA_MANIFEST_SOURCE:
        "21c2830dc65e6be490f3ed262be30800b23cc92da61ea197305fed640fb84dee",
    chain.ACCEPTANCE_SOURCE:
        "71cfd7a62f4bb2c1a4bff045b96bbe78498a9798d30f0f6bd38b6bc6d262c077",
    chain.EFFICIENCY_SOURCE:
        "e3ced212b79d33bbb5bc218648867c2984781efcfb4ca8ac6299baed4f92a256",
    chain.SELECTION_SOURCE:
        "464be3072dbc567ba635c7b4a16d702fd01ad3c3b0c4512cc976141cdea55fd0",
    UPSTREAM / "total/Fit_4D_tot.cpp":
        "1dbb97f7a91ebb6ed364085b20f606404d925a95dfc68154af08c79a53488738",
    UPSTREAM / "total/total_fit_summary.csv":
        "9e0161b54191c927369e989ab3e758daaffd8e544a844c051d11de4f8e3c7e7d",
    UPSTREAM / "total/run.metadata.txt":
        "d68289da64f00dd4d3934e2211c1f767e2cfadf245b33dbc5e3f38efd8ab857f",
}
chain.WEIGHT_DATA_SNAPSHOT_NAME = "WeightData_newdata167_accmix23_effmix19.root"
chain.ACCEPTANCE_SNAPSHOT_NAME = "acceptance_sps0p8_dps0p2.txt"
chain.EFFICIENCY_SNAPSHOT_NAME = "efficiency_sps0p8_dps0p2_avgacc_19x10.txt"
chain.DPS_EXPECTED_FILES = 60
chain.DPS_EXCLUDED_INDICES = {11, 12, 13, 14, 15}
chain.DATA_LABEL_SYMMETRIZE = False
chain.DATA_LABEL_SEED = 50
chain.DPS_LABEL_SYMMETRIZE = False
chain.DPS_LABEL_SEED = 50
chain.EXPECTED_DATA_RANDOMIZATION = None
chain.CR_DEFINITIONS = [
    ("phi70", 1.8, chain.PHI70),
    ("phi110", 1.8, chain.PHI110),
    ("dy2p0", 2.0, chain.PHI90),
]
chain.EXPECTED_PP_FITS = 51
chain.RUNNER_RELATIVE = "Data_driven/run_sps0p8_dps0p2_dedup_nominal_chain.py"
chain.GENERATION_LABEL = "newdata167_sps0p8_dps0p2_dedup60"


if __name__ == "__main__":
    chain.main()
