#!/usr/bin/env bash
set -euo pipefail

REPO="/eos/home-l/leyao/26JJ/JPsiJPsi"
SCRIPT="${REPO}/Data_driven/presentation/scripts/render_reference_f04.sh"
REFERENCE="${REPO}/Data_driven/reference/preupdates_20260829/h015"
ACLIC_DIR="${JJ_DD_F04_ACLIC_DIR:-/tmp/leyao/data_driven_reference_f04_aclic}"
TEMP_EPS="${ACLIC_DIR}/single_jpsi_fit_summary.eps"

if [ "${JJ_DD_F04_INSIDE_EL7:-0}" != "1" ]; then
  /cvmfs/cms.cern.ch/common/cmssw-el7 --command-to-run \
    "export JJ_DD_F04_INSIDE_EL7=1; export JJ_DD_F04_ACLIC_DIR=${ACLIC_DIR}; exec ${SCRIPT}"
  ps2pdf -dEPSCrop "${TEMP_EPS}" \
    "${REPO}/Data_driven/presentation/figures/single_jpsi_fit_summary.pdf"
  test -s "${REPO}/Data_driven/presentation/figures/single_jpsi_fit_summary.pdf"
  test -s "${REPO}/Data_driven/presentation/figures/single_jpsi_fit_summary.png"
  echo "REFERENCE_F04_COMPLETE runtime=ROOT_6.14/09"
  exit 0
fi

cd "${REPO}/Data/ULntuple16/CMSSW_10_6_20/src"
eval "$(scramv1 runtime -sh)"
test "$(root-config --version)" = "6.14/09"
cd "${REPO}"
mkdir -p "${ACLIC_DIR}" "${REPO}/Data_driven/presentation/figures"
root -l -b -q -e \
  "gSystem->SetBuildDir(\"${ACLIC_DIR}\",true); gROOT->ProcessLine(\".L Data_driven/presentation/scripts/plot_single_jpsi_fit_summary.cpp+\"); plot_single_jpsi_fit_summary(\"${TEMP_EPS}\",\"Data_driven/presentation/figures/single_jpsi_fit_summary.png\",\"${REFERENCE}/Model_4D_tot.root\",\"${REFERENCE}/jpsi1_sweights.root\",\"${REFERENCE}/jpsi2_sweights.root\");"
test -s "${TEMP_EPS}"
test -s "${REPO}/Data_driven/presentation/figures/single_jpsi_fit_summary.png"
echo "REFERENCE_F04_DRAW_COMPLETE runtime=ROOT_$(root-config --version)"
