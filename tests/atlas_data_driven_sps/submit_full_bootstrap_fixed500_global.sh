#!/bin/bash
set -euo pipefail

REPO="/eos/home-l/leyao/26JJ/JPsiJPsi"
BASE_SEED="20260806"
RUN_NAME="bootstrap_fixed500_global_nominal_v1"
EOS_OUTPUT_ROOT="${REPO}/tests/atlas_data_driven_sps/results/route9p3_global_calibrated_allpairs/${RUN_NAME}"
AFS_STAGE="/afs/cern.ch/user/l/leyao/private/JJ/JPsiJPsi_${RUN_NAME}"
WORKER_SOURCE="${REPO}/tests/atlas_data_driven_sps/run_full_bootstrap_batch_worker.sh"
WORKER_STAGE="${AFS_STAGE}/run_full_bootstrap_batch_worker.sh"
SUBMIT_FILE="${AFS_STAGE}/full_bootstrap_fixed500_global.sub"

required_macros=(
  build_full_bootstrap_replica_input
  fit_single_jpsi_category_splot
  fit_single_jpsi_global_calibrated_splot
  build_crossslot_splot_mixed
  fit_pp_2d_adaptive
  build_2d_templates_adaptive
)
for macro in "${required_macros[@]}"; do
  source_file="${REPO}/tests/atlas_data_driven_sps/${macro}.cpp"
  library="${REPO}/tests/atlas_data_driven_sps/${macro}_cpp.so"
  if [ ! -f "${library}" ] || [ "${source_file}" -nt "${library}" ]; then
    echo "Missing or stale ROOT 6.14 ACLiC library: ${library}" >&2
    exit 67
  fi
done
if [ -e "${EOS_OUTPUT_ROOT}/submission_metadata.txt" ]; then
  echo "Refusing duplicate global fixed-500 submission: ${EOS_OUTPUT_ROOT}" >&2
  exit 68
fi

mkdir -p "${AFS_STAGE}/logs" "${EOS_OUTPUT_ROOT}"
cp "${WORKER_SOURCE}" "${WORKER_STAGE}"
chmod +x "${WORKER_STAGE}"
sed -e "s|@EXECUTABLE@|${WORKER_STAGE}|g" \
    -e "s|@BASE_SEED@|${BASE_SEED}|g" \
    -e "s|@EOS_OUTPUT_ROOT@|${EOS_OUTPUT_ROOT}|g" \
    -e "s|@LOG_DIR@|${AFS_STAGE}/logs|g" \
    "${REPO}/tests/atlas_data_driven_sps/full_bootstrap_fixed500_global.sub.in" > "${SUBMIT_FILE}"

sha256_of() { sha256sum "$1" | awk '{print $1}'; }
write_submission_metadata() {
  local state="$1" cluster_id="${2:-none}" tmp
  tmp="${EOS_OUTPUT_ROOT}/.submission_metadata.tmp.$$"
  {
    echo "run_name=${RUN_NAME}"
    echo "base_seed=${BASE_SEED}"
    echo "nominal_definition=calibrated_global_allpairs_singleCR_12cell_mergeddy"
    echo "single_j_mode=global_calibrated"
    echo "single_j_fit=two_stage_resolution_calibration_then_frozen_scale_yield_splot"
    echo "binning_mode=nominal12_mergeddy"
    echo "jobs=100"
    echo "replicas_per_job=5"
    echo "max_attempts_per_job=5"
    echo "attempt_stride=10"
    echo "fixed_attempt_ids=job_index_times_10_plus_local_0_to_4"
    echo "replacement_seeds=false"
    echo "target_total_replicas=500"
    echo "fit_single_jpsi_global_sha256=$(sha256_of "${REPO}/tests/atlas_data_driven_sps/fit_single_jpsi_global_calibrated_splot.cpp")"
    echo "fit_single_jpsi_base_sha256=$(sha256_of "${REPO}/tests/atlas_data_driven_sps/fit_single_jpsi_category_splot.cpp")"
    echo "fit_pp_sha256=$(sha256_of "${REPO}/tests/atlas_data_driven_sps/fit_pp_2d_adaptive.cpp")"
    echo "template_builder_sha256=$(sha256_of "${REPO}/tests/atlas_data_driven_sps/build_2d_templates_adaptive.cpp")"
    echo "replica_worker_sha256=$(sha256_of "${REPO}/tests/atlas_data_driven_sps/run_full_bootstrap_replica_worker.sh")"
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
