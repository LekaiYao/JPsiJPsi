#!/bin/bash
set -euo pipefail

if [ "${WEIGHTDATA_INSIDE_EL7:-0}" != "1" ]; then
  exec /cvmfs/cms.cern.ch/common/cmssw-el7 --command-to-run \
    "export WEIGHTDATA_INSIDE_EL7=1; exec $0"
fi

ANALYSIS="/eos/home-l/leyao/26JJ/JPsiJPsi/Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer"
CMSSW_SRC="/eos/home-l/leyao/26JJ/JPsiJPsi/Data/ULntuple16/CMSSW_10_6_20/src"
FINAL="${ANALYSIS}/WeightData.root"
FINAL_META="${ANALYSIS}/WeightData_acc10_eff0p8_dedup60_v1.metadata.txt"
CONDOR_SCRATCH="${_CONDOR_SCRATCH_DIR:-/tmp}"
SCRATCH="${CONDOR_SCRATCH}/jpsijpsi_weightdata_corrected_v1_$$"

if [ -e "${FINAL}" ] || [ -e "${FINAL_META}" ]; then
  echo "Refusing to overwrite existing corrected output" >&2
  exit 17
fi

mkdir -p "${SCRATCH}"
cd "${CMSSW_SRC}"
eval "$(scramv1 runtime -sh)"
cd "${SCRATCH}"

cp "${ANALYSIS}/rephrase.cpp" .
cp "${ANALYSIS}/acceptance_sps_full10_v1.txt" .
cp "${ANALYSIS}/efficiency_sps0p8_dps0p2_dedup60_v1.txt" .

START_EPOCH=$(date +%s)
root -l -b -q rephrase.cpp
END_EPOCH=$(date +%s)

test -s WeightData.root
root -l -b -q -e 'TFile f("WeightData.root","READ"); TTree *t=(TTree*)f.Get("data"); if(f.IsZombie() || !t || t->GetEntries()<=0) gSystem->Exit(2);'

OUTPUT_SHA256=$(sha256sum WeightData.root | awk '{print $1}')
{
  echo "status=complete"
  echo "task=weightdata_acc10_eff0p8_dedup60_v1"
  echo "data_files=119"
  echo "data_suffixes=B20_C9_D14_E3_F8_G29_H36"
  echo "acceptance_file=acceptance_sps_full10_v1.txt"
  echo "acceptance_sha256=0313a2374fd974db331ce55db3958f14ddcaba6279a8cc361a40f11390e8afa6"
  echo "efficiency_file=efficiency_sps0p8_dps0p2_dedup60_v1.txt"
  echo "efficiency_sha256=3d8c72cf35f6704f570c64fb2d4d2b8df6c6bcd83ed301fd9a94dd1c90a23fd6"
  echo "output_sha256=${OUTPUT_SHA256}"
  echo "root_version=$(root-config --version)"
  echo "cmssw_base=${CMSSW_BASE}"
  echo "hostname=$(hostname)"
  echo "start_epoch=${START_EPOCH}"
  echo "end_epoch=${END_EPOCH}"
  echo "wall_seconds=$((END_EPOCH - START_EPOCH))"
} > WeightData_acc10_eff0p8_dedup60_v1.metadata.txt

COPY_ROOT="${ANALYSIS}/.WeightData.root.tmp.$$"
COPY_META="${ANALYSIS}/.WeightData_acc10_eff0p8_dedup60_v1.metadata.txt.tmp.$$"
cp WeightData.root "${COPY_ROOT}"
cp WeightData_acc10_eff0p8_dedup60_v1.metadata.txt "${COPY_META}"
mv "${COPY_ROOT}" "${FINAL}"
mv "${COPY_META}" "${FINAL_META}"

echo "WEIGHTDATA_CORRECTED_V1_COMPLETE output=${FINAL} sha256=${OUTPUT_SHA256}"
