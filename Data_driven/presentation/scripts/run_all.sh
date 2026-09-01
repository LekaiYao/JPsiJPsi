#!/usr/bin/env bash
set -euo pipefail

REPO="/eos/home-l/leyao/26JJ/JPsiJPsi"
PACKAGE="${REPO}/Data_driven/presentation"
cd "${REPO}"

python3 "${PACKAGE}/scripts/plot_figures.py"
"${PACKAGE}/scripts/render_reference_f04.sh"
python3 "${PACKAGE}/scripts/build_package.py"
echo "PRESENTATION_REFERENCE_FIGURES_COMPLETE output=${PACKAGE}/figures"
