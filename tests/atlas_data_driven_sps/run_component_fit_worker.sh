#!/bin/bash
set -euo pipefail
var="$1"
bin="$2"
comp="$3"
base="${4:-tests/atlas_data_driven_sps/results/component_fits_adaptive14}"
input="${5:-tests/atlas_data_driven_sps/results/component_trees_adaptive14/component_data.root}"
out="${base}/jobs/${var}_bin${bin}_${comp}"
mkdir -p "$out"
start_time=$(date +%s)
root -l -b -q "tests/atlas_data_driven_sps/fit_component_variables.cpp+(\"${var}\",${bin},\"${comp}\",\"${out}\",\"${input}\")" > "$out/run.log" 2>&1
result="$out/fit_results.csv"
test -f "$result"
test "$(stat -c %Y "$result")" -ge "$start_time"
grep -q '^single_fit_result_written=' "$out/run.log"
