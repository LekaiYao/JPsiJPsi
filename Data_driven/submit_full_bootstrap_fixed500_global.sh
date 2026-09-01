#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <confirmed-central-tag> <unique-bootstrap-tag>" >&2
  exit 64
fi
CENTRAL_TAG="$1"
BOOTSTRAP_TAG="$2"
for tag in "${CENTRAL_TAG}" "${BOOTSTRAP_TAG}"; do
  [[ "${tag}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || {
    echo "Unsafe tag: ${tag}" >&2
    exit 64
  }
done

REPO="/eos/home-l/leyao/26JJ/JPsiJPsi"
source "${REPO}/Data_driven/require_root640.sh"
BASE_SEED="${JJ_DD_BOOTSTRAP_BASE_SEED:-20260806}"
CENTRAL="${REPO}/Data_driven/results/${CENTRAL_TAG}"
EOS_OUTPUT_ROOT="${REPO}/Data_driven/results/${BOOTSTRAP_TAG}"
MODEL="${REPO}/Data_driven/inputs/Model_4D_tot.root"
SOURCE_CANDIDATES="${CENTRAL}/input/jpsi12_candidates.root"
SOURCE_DATA="${CENTRAL}/input/WeightData_mJJ7p5.root"
AFS_STAGE="/afs/cern.ch/user/l/leyao/private/JJ/JPsiJPsi_${BOOTSTRAP_TAG}"
WORKER_SOURCE="${REPO}/Data_driven/run_full_bootstrap_batch_worker.sh"
WORKER_STAGE="${AFS_STAGE}/run_full_bootstrap_batch_worker.sh"
SUBMIT_FILE="${AFS_STAGE}/full_bootstrap_fixed500_global.sub"

for required in "${MODEL}" "${SOURCE_CANDIDATES}" "${SOURCE_DATA}"; do
  test -s "${required}" || { echo "Missing required input: ${required}" >&2; exit 66; }
done
test -f "${CENTRAL}/complete.marker" || {
  echo "Central tag has no completion marker: ${CENTRAL}" >&2
  exit 66
}
if [ -e "${EOS_OUTPUT_ROOT}" ]; then
  echo "Refusing duplicate fixed-500 bootstrap tag: ${EOS_OUTPUT_ROOT}" >&2
  exit 68
fi

required_macros=(
  build_full_bootstrap_replica_input
  fit_single_jpsi_category_splot
  fit_single_jpsi_global_calibrated_splot
  build_crossslot_splot_mixed
  fit_pp_2d_adaptive
  build_2d_templates_adaptive
)
cd "${REPO}"
for macro in "${required_macros[@]}"; do
  source_file="${REPO}/Data_driven/${macro}.cpp"
  library="${REPO}/Data_driven/${macro}_cpp.so"
  if [ ! -f "${library}" ] || [ "${source_file}" -nt "${library}" ]; then
    root -l -b -q -e "gROOT->ProcessLine(\".L Data_driven/${macro}.cpp+\");"
  fi
  test -s "${library}" || {
    echo "ROOT 6.40.02 ACLiC compilation failed: ${source_file}" >&2
    exit 67
  }
done

mkdir -p "${AFS_STAGE}/logs" "${EOS_OUTPUT_ROOT}"
cp "${WORKER_SOURCE}" "${WORKER_STAGE}"
chmod +x "${WORKER_STAGE}"
sed -e "s|@EXECUTABLE@|${WORKER_STAGE}|g" \
    -e "s|@BASE_SEED@|${BASE_SEED}|g" \
    -e "s|@EOS_OUTPUT_ROOT@|${EOS_OUTPUT_ROOT}|g" \
    -e "s|@LOG_DIR@|${AFS_STAGE}/logs|g" \
    -e "s|@MODEL@|${MODEL}|g" \
    -e "s|@SOURCE_CANDIDATES@|${SOURCE_CANDIDATES}|g" \
    -e "s|@SOURCE_DATA@|${SOURCE_DATA}|g" \
    "${REPO}/Data_driven/full_bootstrap_fixed500_global.sub.in" > "${SUBMIT_FILE}"

sha256_of() { sha256sum "$1" | awk '{print $1}'; }
write_submission_metadata() {
  local state="$1" cluster_id="${2:-none}" tmp
  tmp="${EOS_OUTPUT_ROOT}/.submission_metadata.tmp.$$"
  {
    echo "run_tag=${BOOTSTRAP_TAG}"
    echo "central_tag=${CENTRAL_TAG}"
    echo "base_seed=${BASE_SEED}"
    echo "root_version=$(root-config --version)"
    echo "nominal_definition=global_calibrated_allpairs_singleCR_12cell_mergeddy"
    echo "jobs=100"
    echo "replicas_per_job=5"
    echo "max_attempts_per_job=5"
    echo "attempt_stride=10"
    echo "replacement_seeds=false"
    echo "target_total_replicas=500"
    echo "model=${MODEL}"
    echo "model_sha256=$(sha256_of "${MODEL}")"
    echo "source_candidates=${SOURCE_CANDIDATES}"
    echo "source_candidates_sha256=$(sha256_of "${SOURCE_CANDIDATES}")"
    echo "source_data=${SOURCE_DATA}"
    echo "source_data_sha256=$(sha256_of "${SOURCE_DATA}")"
    echo "replica_worker_sha256=$(sha256_of "${REPO}/Data_driven/run_full_bootstrap_replica_worker.sh")"
    echo "batch_worker_sha256=$(sha256_of "${WORKER_SOURCE}")"
    echo "schedd=bigbird08.cern.ch"
    echo "submission_state=${state}"
    echo "cluster_id=${cluster_id}"
  } > "${tmp}"
  mv "${tmp}" "${EOS_OUTPUT_ROOT}/submission_metadata.txt"
}

write_submission_metadata prepared
cd "${AFS_STAGE}"
submit_output=$(condor_submit -name bigbird08.cern.ch "${SUBMIT_FILE}")
echo "${submit_output}"
cluster_id=$(printf '%s\n' "${submit_output}" | sed -n 's/.*cluster \([0-9][0-9]*\).*/\1/p' | tail -n 1)
if [ -z "${cluster_id}" ]; then
  echo "Submission returned no cluster id" >&2
  exit 69
fi
write_submission_metadata submitted "${cluster_id}"
