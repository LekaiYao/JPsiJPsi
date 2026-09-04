#!/usr/bin/env bash
set -euo pipefail

run_tag="${1:-newdata167_accmix23_effmix19_sps0p85_dps0p15_lumi36p684_pairsym_cseed50_root640_strategy2_nominal_v4_20260902}"
if [[ ! "${run_tag}" =~ ^[A-Za-z0-9._-]+$ ]]; then
    echo "Invalid run tag: ${run_tag}" >&2
    exit 2
fi

analysis="/eos/home-l/leyao/26JJ/JPsiJPsi/Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer"
source_weight_tag="newdata167_accmix23_effmix19_sps0p85_dps0p15_lumi36p684_pairsym_cseed50_root640_strategy2_nominal_v1_20260902"
source_input_dir="${analysis}/fit_results/${source_weight_tag}/inputs"
result_dir="${analysis}/fit_results/${run_tag}"
input_dir="${result_dir}/inputs"
total_dir="${result_dir}/total"
weight_name="WeightData_newdata167_accmix23_effmix19_sps0p85_dps0p15_pairsym_cseed50_v1_20260902.root"

if [[ "$(root-config --version)" != "6.40.02" ]]; then
    echo "This workflow requires ROOT 6.40.02" >&2
    exit 3
fi
if [[ ! -s "${source_input_dir}/${weight_name}" ]]; then
    echo "Missing source WeightData: ${source_input_dir}/${weight_name}" >&2
    exit 4
fi
if [[ -e "${result_dir}" ]]; then
    echo "Refusing to reuse total-fit tag: ${result_dir}" >&2
    exit 17
fi

mkdir -p "${result_dir}"
cp -a "${source_input_dir}" "${result_dir}/inputs"
mkdir -p "${total_dir}/logs"
cp "${analysis}/Fit_4D_tot.cpp" "${total_dir}/"
cp "${analysis}/Plot_4D.hpp" "${total_dir}/"
ln -s "../inputs/${weight_name}" "${total_dir}/WeightData.root"

cd "${total_dir}"
start_epoch="$(date +%s)"
root -l -b -q 'Fit_4D_tot.cpp(false,false,false,true,false,2)' > logs/nominal.log 2>&1
end_epoch="$(date +%s)"

test -s Model_4D_tot.root
test -s Fit_4D_tot_native_asymptotic_unseeded.root
printf '%s\n' 'fit,status,covQual,edm,n_P_P,n_P_P_error' > total_fit_summary.csv
root -l -b -q -e 'TFile f("Fit_4D_tot_native_asymptotic_unseeded.root","READ"); RooFitResult *r=dynamic_cast<RooFitResult *>(f.Get("fit_result_native_asymptotic")); RooRealVar *n=r ? dynamic_cast<RooRealVar *>(r->floatParsFinal().find("n_P_P")) : nullptr; if(f.IsZombie() || !r || !n || r->status()!=0 || r->covQual()!=3 || r->edm()>=0.01) gSystem->Exit(5); printf("FIT_SUMMARY,nominal,%d,%d,%.17g,%.17g,%.17g\n",r->status(),r->covQual(),r->edm(),n->getVal(),n->getError()); gSystem->Exit(0);' 2>/dev/null | sed -n 's/^FIT_SUMMARY,//p' >> total_fit_summary.csv
[[ "$(wc -l < total_fit_summary.csv)" -eq 2 ]]

printf '%s\n' 'parameter,value,error,minimum,maximum,distance_to_nearest_boundary_in_errors' > fit_parameters.csv
root -l -b -q -e 'TFile f("Fit_4D_tot_native_asymptotic_unseeded.root","READ"); RooFitResult *r=dynamic_cast<RooFitResult *>(f.Get("fit_result_native_asymptotic")); if(!r) gSystem->Exit(5); const RooArgList &p=r->floatParsFinal(); for(int i=0;i<p.getSize();++i){RooRealVar *v=dynamic_cast<RooRealVar *>(p.at(i)); double d=TMath::Min(v->getVal()-v->getMin(),v->getMax()-v->getVal())/v->getError(); printf("FIT_PARAM,%s,%.17g,%.17g,%.17g,%.17g,%.17g\n",v->GetName(),v->getVal(),v->getError(),v->getMin(),v->getMax(),d);} gSystem->Exit(0);' 2>/dev/null | sed -n 's/^FIT_PARAM,//p' >> fit_parameters.csv
[[ "$(wc -l < fit_parameters.csv)" -gt 2 ]]

plot_count="$(find fig/native_asymptotic/unseeded -maxdepth 1 -type f -name '*.pdf' | wc -l)"
[[ "${plot_count}" -eq 4 ]]
{
    echo "status=complete"
    echo "task=total_4d_newdata167_mix0p85_pairsym_cseed50_strategy2_nominal"
    echo "run_tag=${run_tag}"
    echo "source_weight_tag=${source_weight_tag}"
    echo "input_label_randomization=nominal chensh-compatible C srand/rand seed 50, reseeded per input ROOT file"
    echo "error_convention=ROOT native AsymptoticError(true), no external seed"
    echo "minimizer_strategy=2"
    echo "central_prefit=same-process weighted extended-likelihood prefit with SumW2Error(false), Strategy(2), no external seed, up to 6 attempts; discarded after establishing the endpoint"
    echo "corrected_fit_retry=up to 6 native AsymptoticError(true) attempts from the preceding in-memory endpoint, Strategy(2), Offset(false)"
    echo "fit_scope=nominal Crystal-Ball-plus-Gaussian model only; alternate Gaussian-core reference not requested or run"
    echo "fit_gate=status=0,covQual=3,EDM<0.01"
    echo "luminosity_fb=36.684 (not used by the likelihood)"
    echo "weightdata=${input_dir}/${weight_name}"
    echo "weightdata_sha256=$(sha256sum "${input_dir}/${weight_name}" | awk '{print $1}')"
    echo "acceptance_sha256=$(sha256sum "${input_dir}/acceptance_sps0p85_dps0p15.txt" | awk '{print $1}')"
    echo "efficiency_sha256=$(sha256sum "${input_dir}/efficiency_sps0p85_dps0p15_avgacc_19x10.txt" | awk '{print $1}')"
    echo "root_version=$(root-config --version)"
    echo "git_head=$(git -C /eos/home-l/leyao/26JJ/JPsiJPsi rev-parse HEAD)"
    echo "hostname=$(hostname)"
    echo "start_epoch=${start_epoch}"
    echo "end_epoch=${end_epoch}"
    echo "wall_seconds=$((end_epoch - start_epoch))"
    sha256sum Fit_4D_tot.cpp Plot_4D.hpp Model_4D_tot.root \
        Fit_4D_tot_native_asymptotic_unseeded.root total_fit_summary.csv fit_parameters.csv
} > run.metadata.txt

echo "TOTAL4D_MIX0P85_PAIRSYM_COMPLETE result=${total_dir}/total_fit_summary.csv"
