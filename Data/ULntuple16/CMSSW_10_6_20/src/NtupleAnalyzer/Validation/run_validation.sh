#!/bin/bash
# Regenerate the 7 Appendix C MC-vs-data validation figures with the nominal NLO* SPS MC.
# Each line: DATAVAR | MCVAR | TITLE | NBINS | XMIN | XMAX
set -e
mkdir -p output_sbs

VARS=(
  "DeltaPhi|DeltaPhi_|#Delta#phi(J/#psi_{1},J/#psi_{2})|16|0|3.14159"
  "DeltaY|DeltaY_|#Deltay(J/#psi_{1},J/#psi_{2})|20|0|4.0"
  "FourMuonMass|FourMuonMass_|M(J/#psi J/#psi) [GeV]|40|7.5|107.5"
  "J1Pt|J1PtReco_|p_{T}(J/#psi_{1}) [GeV]|30|10|40"
  "J1y|J1yReco_|y(J/#psi_{1})|20|-2.0|2.0"
  "Mu1Pt|Mu1PtReco_|p_{T}(#mu_{1}) [GeV]|30|0|40"
  "Mu1Eta|Mu1EtaReco_|#eta(#mu_{1})|22|-2.2|2.2"
)

for entry in "${VARS[@]}"; do
  IFS='|' read -r DATAVAR MCVAR TITLE NBINS XMIN XMAX <<< "$entry"
  echo "Processing: $DATAVAR (MC branch $MCVAR)"
  script="plot_${DATAVAR}.C"
  sed -e "s|@DATAVAR@|${DATAVAR}|g" \
      -e "s|@MCVAR@|${MCVAR}|g" \
      -e "s|@TITLE@|${TITLE}|g" \
      -e "s|@NBINS@|${NBINS}|g" \
      -e "s|@XMIN@|${XMIN}|g" \
      -e "s|@XMAX@|${XMAX}|g" \
      plot_validation.C > "$script"
  root -l -b -q "$script"
  rm "$script"
done
echo "Done. Figures in output_sbs/"
