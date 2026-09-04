#!/usr/bin/env bash
set -euo pipefail

run_tag="${1:-newdata167_accmix23_effmix19_lumi36p684_pairsym_cseed50_root640_v1_20260902}"
if [[ ! "${run_tag}" =~ ^[A-Za-z0-9._-]+$ ]]; then
    echo "Invalid run tag: ${run_tag}" >&2
    exit 2
fi

script_path="$(readlink -f "$0")"
if [[ "${WEIGHTDATA_PAIRSYM_INSIDE_EL7:-0}" != "1" ]]; then
    exec /cvmfs/cms.cern.ch/common/cmssw-el7 --command-to-run \
        "export WEIGHTDATA_PAIRSYM_INSIDE_EL7=1; exec ${script_path} ${run_tag}"
fi

analysis="/eos/home-l/leyao/26JJ/JPsiJPsi/Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer"
cmssw_src="/eos/home-l/leyao/26JJ/JPsiJPsi/Data/ULntuple16/CMSSW_10_6_20/src"
result_dir="${analysis}/fit_results/${run_tag}"
input_dir="${result_dir}/inputs"
log_dir="${result_dir}/logs"
weight_name="WeightData_newdata167_accmix23_effmix19_pairsym_cseed50_v1_20260902.root"
metadata_name="WeightData_newdata167_accmix23_effmix19_pairsym_cseed50_v1_20260902.metadata.txt"
acceptance="/eos/home-l/leyao/26JJ/JPsiJPsi/GEN_nofilter/DPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc/closure_results/nominal_mixed23_sps0p8_dps0p2_avgacc_v1_20260902/acceptance_sps0p8_dps0p2.txt"
efficiency="/eos/home-l/leyao/26JJ/JPsiJPsi/GEN_nofilter/DPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc/closure_results/inputs_mixed23_sps0p8_dps0p2_avgacc_v1_20260902/efficiency_sps0p8_dps0p2_avgacc_19x10.txt"
provenance="/eos/home-l/leyao/26JJ/JPsiJPsi/GEN_nofilter/DPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc/closure_inputs_mixed23_sps0p8_dps0p2_avgacc_v1.json"
nominal_dir="/eos/home-l/leyao/26JJ/JPsiJPsi/GEN_nofilter/DPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc/closure_results/nominal_mixed23_sps0p8_dps0p2_avgacc_v1_20260902"
correction_input_dir="/eos/home-l/leyao/26JJ/JPsiJPsi/GEN_nofilter/DPS/ULPythia2016/CMSSW_10_2_5/src/4mu_acc/closure_results/inputs_mixed23_sps0p8_dps0p2_avgacc_v1_20260902"
chensh_reference="/eos/user/c/chensh/JPsiJPsi/Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/rephrase.cpp"

if [[ -e "${result_dir}" ]]; then
    echo "Refusing to reuse an existing run tag: ${result_dir}" >&2
    exit 17
fi

acceptance_sha="$(sha256sum "${acceptance}" | awk '{print $1}')"
efficiency_sha="$(sha256sum "${efficiency}" | awk '{print $1}')"
chensh_reference_sha="$(sha256sum "${chensh_reference}" | awk '{print $1}')"
[[ "${acceptance_sha}" == "71cfd7a62f4bb2c1a4bff045b96bbe78498a9798d30f0f6bd38b6bc6d262c077" ]]
[[ "${efficiency_sha}" == "e3ced212b79d33bbb5bc218648867c2984781efcfb4ca8ac6299baed4f92a256" ]]
[[ "${chensh_reference_sha}" == "b6bd71f8bc03f0a6a316b1c5e43180cf9f0ebecb624adeb326ce3bb85252140b" ]]

mkdir -p "${input_dir}/provenance" "${log_dir}"
cp "${analysis}/rephrase.cpp" "${input_dir}/rephrase.cpp"
cp "${analysis}/data_ntuple_manifest_chensh_newdata_v1.list" "${input_dir}/data_ntuple_manifest_chensh_newdata_v1.list"
cp "${chensh_reference}" "${input_dir}/provenance/rephrase_chensh_reference.cpp"
cp "${acceptance}" "${input_dir}/acceptance_sps0p8_dps0p2.txt"
cp "${efficiency}" "${input_dir}/efficiency_sps0p8_dps0p2_avgacc_19x10.txt"
cp "${provenance}" "${input_dir}/provenance/"
cp "${nominal_dir}/metadata.json" "${input_dir}/provenance/nominal_metadata.json"
cp "${nominal_dir}/checksums.sha256" "${input_dir}/provenance/nominal_checksums.sha256"
cp "${nominal_dir}/sys3.csv" "${input_dir}/provenance/"
cp "${nominal_dir}/sys3_cpp.txt" "${input_dir}/provenance/"
cp "${correction_input_dir}/metadata.json" "${input_dir}/provenance/correction_inputs_metadata.json"
cp "${correction_input_dir}/checksums.sha256" "${input_dir}/provenance/correction_inputs_checksums.sha256"

scratch_parent="${_CONDOR_SCRATCH_DIR:-/tmp}"
scratch="$(mktemp -d "${scratch_parent}/jpsijpsi_weightdata_pairsym_cseed50.XXXXXX")"
cd "${cmssw_src}"
eval "$(scramv1 runtime -sh)"
cd "${scratch}"
cp "${input_dir}/rephrase.cpp" .
cp "${input_dir}/data_ntuple_manifest_chensh_newdata_v1.list" .

start_epoch="$(date +%s)"
root -l -b -q 'rephrase.cpp(true,50)' > "${log_dir}/rephrase.log" 2>&1
end_epoch="$(date +%s)"

test -s WeightData.root
root -l -b -q -e 'TFile f("WeightData.root","READ"); TTree *t=dynamic_cast<TTree *>(f.Get("data")); if(f.IsZombie() || !t || !t->GetBranch("evt_weight")) gSystem->Exit(2); const Long64_t n=t->GetEntries(); const Long64_t nonfinite=t->GetEntries("!TMath::Finite(evt_weight)"); const Long64_t nonpositive=t->GetEntries("evt_weight<=0"); const Long64_t pt1gt=t->GetEntries("Jpsi_pt1>Jpsi_pt2"); const Long64_t pt1eq=t->GetEntries("Jpsi_pt1==Jpsi_pt2"); const double minimum=t->GetMinimum("evt_weight"); const double maximum=t->GetMaximum("evt_weight"); printf("entries=%lld\nnonfinite_weights=%lld\nnonpositive_weights=%lld\nminimum_weight=%.17g\nmaximum_weight=%.17g\npt1_gt_pt2=%lld\npt1_eq_pt2=%lld\npt1_gt_pt2_fraction=%.17g\n",n,nonfinite,nonpositive,minimum,maximum,pt1gt,pt1eq,n>0?(double)pt1gt/n:0.0); if(n!=2954 || nonfinite || nonpositive || maximum>1000.) gSystem->Exit(3);' > "${log_dir}/weightdata_qa.txt" 2>&1

randomization_draws="$(sed -n 's/^Randomization draws after HLT+trigger matching: //p' "${log_dir}/rephrase.log")"
candidate_swaps="$(sed -n 's/^Randomized candidate swaps: //p' "${log_dir}/rephrase.log")"
selected_swaps="$(sed -n 's/^Randomized selected-event swaps: //p' "${log_dir}/rephrase.log")"
test -n "${randomization_draws}"
test -n "${candidate_swaps}"
test -n "${selected_swaps}"

output_sha="$(sha256sum WeightData.root | awk '{print $1}')"
cp WeightData.root "${input_dir}/.${weight_name}.tmp"
mv "${input_dir}/.${weight_name}.tmp" "${input_dir}/${weight_name}"

{
    echo "status=complete"
    echo "task=weightdata_newdata167_accmix23_effmix19_pairsym_cseed50_v1"
    echo "run_tag=${run_tag}"
    echo "data_manifest=data_ntuple_manifest_chensh_newdata_v1.list"
    echo "data_files=167"
    echo "data_manifest_sha256=$(sha256sum "${analysis}/data_ntuple_manifest_chensh_newdata_v1.list" | awk '{print $1}')"
    echo "label_randomization=enabled"
    echo "label_randomization_reference=${chensh_reference}"
    echo "label_randomization_reference_sha256=${chensh_reference_sha}"
    echo "label_randomization_rng=C library srand/rand"
    echo "label_randomization_seed=50"
    echo "label_randomization_scope=reseed at the start of every input ROOT file"
    echo "label_randomization_draw_timing=after passHLT and matchTrg, before pT/samePV/pair-mass selection, branch assignment, and calWeight"
    echo "label_randomization_rule=keep original IDs when rand/(RAND_MAX+1.0)>0.5; otherwise swap"
    echo "randomization_draws=${randomization_draws}"
    echo "randomized_candidate_swaps=${candidate_swaps}"
    echo "randomized_selected_event_swaps=${selected_swaps}"
    echo "acceptance_file=${acceptance}"
    echo "acceptance_sha256=${acceptance_sha}"
    echo "efficiency_file=${efficiency}"
    echo "efficiency_sha256=${efficiency_sha}"
    echo "output=${input_dir}/${weight_name}"
    echo "output_sha256=${output_sha}"
    echo "root_version=$(root-config --version)"
    echo "cmssw_base=${CMSSW_BASE}"
    echo "hostname=$(hostname)"
    echo "start_epoch=${start_epoch}"
    echo "end_epoch=${end_epoch}"
    echo "wall_seconds=$((end_epoch - start_epoch))"
    sha256sum "${input_dir}/rephrase.cpp" "${input_dir}/provenance/rephrase_chensh_reference.cpp" "${input_dir}/${weight_name}" "${log_dir}/weightdata_qa.txt"
} > "${input_dir}/${metadata_name}"

echo "WEIGHTDATA_PAIRSYM_COMPLETE output=${input_dir}/${weight_name} sha256=${output_sha}"
