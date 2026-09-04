#!/usr/bin/env bash
set -euo pipefail

run_tag="${1:-newdata167_accmix23_effmix19_lumi36p684_pairsym_cseed50_root640_v1_20260902}"
if [[ ! "${run_tag}" =~ ^[A-Za-z0-9._-]+$ ]]; then
    echo "Invalid run tag: ${run_tag}" >&2
    exit 2
fi

analysis="/eos/home-l/leyao/26JJ/JPsiJPsi/Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer"
result_dir="${analysis}/fit_results/${run_tag}"
input_dir="${result_dir}/inputs"
total_dir="${result_dir}/total"
weight_name="WeightData_newdata167_accmix23_effmix19_pairsym_cseed50_v1_20260902.root"

if [[ "$(root-config --version)" != "6.40.02" ]]; then
    echo "This workflow requires ROOT 6.40.02" >&2
    exit 3
fi
if [[ ! -s "${input_dir}/${weight_name}" ]]; then
    echo "Missing tagged WeightData: ${input_dir}/${weight_name}" >&2
    exit 4
fi
if [[ -e "${total_dir}" ]]; then
    echo "Refusing to reuse total-fit directory: ${total_dir}" >&2
    exit 17
fi

mkdir -p "${total_dir}/logs"
cp "${analysis}/Fit_4D_tot.cpp" "${total_dir}/"
cp "${analysis}/Plot_4D.hpp" "${total_dir}/"
ln -s "../inputs/${weight_name}" "${total_dir}/WeightData.root"

cd "${total_dir}"
start_epoch="$(date +%s)"
root -l -b -q 'Fit_4D_tot.cpp(false,false,false,true,false)' > logs/nominal.log 2>&1
root -l -b -q 'Fit_4D_tot.cpp(true,false,false,true,false)' > logs/reference.log 2>&1
end_epoch="$(date +%s)"

for required in Model_4D_tot.root Model_4D_tot_ref.root \
    Fit_4D_tot_native_asymptotic_unseeded.root \
    Fit_4D_tot_ref_native_asymptotic_unseeded.root; do
    test -s "${required}"
done

printf '%s\n' 'fit,status,covQual,edm,n_P_P,n_P_P_error' > total_fit_summary.csv
for spec in \
    'nominal Fit_4D_tot_native_asymptotic_unseeded.root' \
    'reference Fit_4D_tot_ref_native_asymptotic_unseeded.root'; do
    read -r label file_name <<< "${spec}"
    root -l -b -q -e "TFile f(\"${file_name}\",\"READ\"); RooFitResult *r=dynamic_cast<RooFitResult *>(f.Get(\"fit_result_native_asymptotic\")); RooRealVar *n=r ? dynamic_cast<RooRealVar *>(r->floatParsFinal().find(\"n_P_P\")) : nullptr; if(f.IsZombie() || !r || !n || r->status()!=0 || r->covQual()!=3 || r->edm()>=0.01) gSystem->Exit(5); printf(\"FIT_SUMMARY,${label},%d,%d,%.17g,%.17g,%.17g\\n\",r->status(),r->covQual(),r->edm(),n->getVal(),n->getError()); gSystem->Exit(0);" 2>/dev/null | sed -n 's/^FIT_SUMMARY,//p' >> total_fit_summary.csv
done
[[ "$(wc -l < total_fit_summary.csv)" -eq 3 ]]

{
    echo "status=complete"
    echo "task=total_4d_newdata167_accmix23_effmix19_pairsym_cseed50_root640_v1"
    echo "run_tag=${run_tag}"
    echo "input_label_randomization=chensh-compatible C srand/rand seed 50, reseeded per input ROOT file"
    echo "error_convention=ROOT native AsymptoticError(true), no external seed"
    echo "fits=nominal,reference"
    echo "fit_gate=status=0,covQual=3,EDM<0.01"
    echo "luminosity_fb=36.684 (plot annotation and downstream cross-section input; not used by the likelihood)"
    echo "weightdata=${input_dir}/${weight_name}"
    echo "weightdata_sha256=$(sha256sum "${input_dir}/${weight_name}" | awk '{print $1}')"
    echo "root_version=$(root-config --version)"
    echo "git_head=$(git -C /eos/home-l/leyao/26JJ/JPsiJPsi rev-parse HEAD)"
    echo "hostname=$(hostname)"
    echo "start_epoch=${start_epoch}"
    echo "end_epoch=${end_epoch}"
    echo "wall_seconds=$((end_epoch - start_epoch))"
    sha256sum Fit_4D_tot.cpp Plot_4D.hpp \
        Model_4D_tot.root Model_4D_tot_ref.root \
        Fit_4D_tot_native_asymptotic_unseeded.root \
        Fit_4D_tot_ref_native_asymptotic_unseeded.root \
        total_fit_summary.csv
} > run.metadata.txt

echo "TOTAL4D_PAIRSYM_COMPLETE result=${total_dir}/total_fit_summary.csv"
