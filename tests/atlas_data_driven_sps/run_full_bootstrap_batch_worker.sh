#!/bin/bash
set -euo pipefail

if [ "$#" -lt 5 ] || [ "$#" -gt 8 ]; then
  echo "Usage: $0 <job_index> <base_seed> <eos_output_root> <target_successes> <max_attempts> [legacy14|nominal12_mergeddy] [attempt_stride] [absy1p2|global_calibrated]" >&2
  exit 64
fi
if [ "${FULL_BOOTSTRAP_BATCH_INSIDE_EL7:-0}" != "1" ]; then
  exec /cvmfs/cms.cern.ch/common/cmssw-el7 --command-to-run \
    "export FULL_BOOTSTRAP_BATCH_INSIDE_EL7=1; exec $0 $1 $2 $3 $4 $5 ${6:-legacy14} ${7:-$5} ${8:-absy1p2}"
fi

JOB_INDEX="$1"; BASE_SEED="$2"; EOS_ROOT="$3"
TARGET_SUCCESSES="$4"; MAX_ATTEMPTS="$5"
BINNING_MODE="${6:-legacy14}"
ATTEMPT_STRIDE="${7:-${MAX_ATTEMPTS}}"
SINGLE_J_MODE="${8:-absy1p2}"
case "${BINNING_MODE}" in
  legacy14|nominal12_mergeddy) ;;
  *) echo "Unsupported binning mode: ${BINNING_MODE}" >&2; exit 64 ;;
esac
case "${SINGLE_J_MODE}" in
  absy1p2|global_calibrated) ;;
  *) echo "Unsupported single-J mode: ${SINGLE_J_MODE}" >&2; exit 64 ;;
esac
for value in "${JOB_INDEX}" "${BASE_SEED}" "${TARGET_SUCCESSES}" "${MAX_ATTEMPTS}" "${ATTEMPT_STRIDE}"; do
  [[ "${value}" =~ ^[0-9]+$ ]] || { echo "integer argument required: ${value}" >&2; exit 64; }
done
if [ "${TARGET_SUCCESSES}" -lt 1 ] || [ "${MAX_ATTEMPTS}" -lt "${TARGET_SUCCESSES}" ]; then
  echo "Require 1 <= target_successes <= max_attempts" >&2
  exit 64
fi
if [ "${ATTEMPT_STRIDE}" -lt "${MAX_ATTEMPTS}" ]; then
  echo "Require attempt_stride >= max_attempts" >&2
  exit 64
fi

REPO="/eos/home-l/leyao/26JJ/JPsiJPsi"
CMSSW_SRC="${REPO}/Data/ULntuple16/CMSSW_10_6_20/src"
REPLICA_WORKER="${FULL_BOOTSTRAP_REPLICA_WORKER_OVERRIDE:-${REPO}/tests/atlas_data_driven_sps/run_full_bootstrap_replica_worker.sh}"
JOB_TAG=$(printf 'job_%03d' "${JOB_INDEX}")
JOB_EOS="${EOS_ROOT}/${JOB_TAG}"
CONDOR_SCRATCH="${_CONDOR_SCRATCH_DIR:-/tmp}"
JOB_SCRATCH="${CONDOR_SCRATCH}/jpsijpsi_bootstrap_${JOB_TAG}"

case "${JOB_SCRATCH}" in
  /tmp/jpsijpsi_bootstrap_job_*|/pool/condor/*/jpsijpsi_bootstrap_job_*) ;;
  *) echo "Unsafe scratch path: ${JOB_SCRATCH}" >&2; exit 66 ;;
esac
test -x "${REPLICA_WORKER}" || { echo "Replica worker is not executable: ${REPLICA_WORKER}" >&2; exit 66; }
mkdir -p "${JOB_EOS}" "${JOB_SCRATCH}"
cd "${CMSSW_SRC}"; eval "$(scramv1 runtime -sh)"; cd "${REPO}"

write_job_summary() {
  local state="$1" end_epoch tmp
  end_epoch=$(date +%s)
  tmp="${JOB_EOS}/.job_summary.tmp.$$"
  {
    echo "job_index=${JOB_INDEX}"
    echo "target_successes=${TARGET_SUCCESSES}"
    echo "max_attempts=${MAX_ATTEMPTS}"
    echo "attempt_stride=${ATTEMPT_STRIDE}"
    echo "binning_mode=${BINNING_MODE}"
    echo "single_j_mode=${SINGLE_J_MODE}"
    echo "attempted=${ATTEMPTED}"
    echo "successes=${SUCCESSES}"
    echo "failures=${FAILURES}"
    echo "hostname=$(hostname)"
    echo "root_version=$(root-config --version)"
    echo "state=${state}"
    echo "updated_epoch=${end_epoch}"
  } > "${tmp}"
  mv "${tmp}" "${JOB_EOS}/job_summary.txt"
}

checkpoint_attempt() {
  local scratch="$1" final="$2" status="$3" exit_code="$4" attempt_id="$5" local_attempt="$6"
  local tmp="${JOB_EOS}/.attempt_${attempt_id}.tmp.$$" rel
  mkdir -p "${tmp}"
  for rel in \
    metadata.txt timing.csv \
    input/metadata.txt input/source_event_multiplicities.csv \
    splot/category_splot_summary_jpsi1.txt splot/category_splot_summary_jpsi2.txt \
    splot/global_calibrated_splot_summary_jpsi1.txt splot/global_calibrated_splot_summary_jpsi2.txt \
    pp_data_2d/summary.txt pp_data_2d/fit_results.csv \
    templates_2d/summary.txt templates_2d/cells.csv \
    fit_jpsi1.log fit_jpsi2.log build_mixing.log fit_pp.log build_templates.log worker.log; do
    if [ -f "${scratch}/${rel}" ]; then
      mkdir -p "${tmp}/$(dirname "${rel}")"
      cp "${scratch}/${rel}" "${tmp}/${rel}"
    fi
  done
  {
    echo "status=${status}"
    echo "exit_code=${exit_code}"
    echo "job_index=${JOB_INDEX}"
    echo "local_attempt=${local_attempt}"
    echo "attempt_id=${attempt_id}"
    echo "base_seed=${BASE_SEED}"
    echo "seed=$((BASE_SEED + 1000003 * attempt_id))"
  } > "${tmp}/attempt_status.txt"
  mv "${tmp}" "${final}"
}

SUCCESSES=0; FAILURES=0; ATTEMPTED=0
for ((local_attempt=0; local_attempt<MAX_ATTEMPTS; ++local_attempt)); do
  attempt_id=$((JOB_INDEX * ATTEMPT_STRIDE + local_attempt))
  attempt_tag=$(printf 'attempt_%06d' "${attempt_id}")
  final_attempt="${JOB_EOS}/${attempt_tag}"
  if [ -f "${final_attempt}/attempt_status.txt" ]; then
    status=$(awk -F= '$1=="status"{print $2}' "${final_attempt}/attempt_status.txt")
    ATTEMPTED=$((ATTEMPTED + 1))
    if [ "${status}" = "success" ]; then SUCCESSES=$((SUCCESSES + 1)); else FAILURES=$((FAILURES + 1)); fi
  fi
done
if [ "${SUCCESSES}" -ge "${TARGET_SUCCESSES}" ]; then
  write_job_summary complete
  echo "BOOTSTRAP_BATCH_ALREADY_COMPLETE job=${JOB_INDEX} successes=${SUCCESSES} attempted=${ATTEMPTED}"
  exit 0
fi
write_job_summary running

for ((local_attempt=0; local_attempt<MAX_ATTEMPTS && SUCCESSES<TARGET_SUCCESSES; ++local_attempt)); do
  attempt_id=$((JOB_INDEX * ATTEMPT_STRIDE + local_attempt))
  attempt_tag=$(printf 'attempt_%06d' "${attempt_id}")
  final_attempt="${JOB_EOS}/${attempt_tag}"
  [ -f "${final_attempt}/attempt_status.txt" ] && continue
  scratch_attempt="${JOB_SCRATCH}/${attempt_tag}"
  mkdir -p "${scratch_attempt}"
  set +e
  FULL_BOOTSTRAP_INSIDE_EL7=1 "${REPLICA_WORKER}" \
    "${attempt_id}" "${BASE_SEED}" "${scratch_attempt}" "${BINNING_MODE}" "${SINGLE_J_MODE}" > "${scratch_attempt}/worker.log" 2>&1
  exit_code=$?
  set -e
  status=failed
  if [ "${exit_code}" -eq 0 ] && grep -q '^status=complete$' "${scratch_attempt}/metadata.txt" 2>/dev/null; then
    status=success
  fi
  checkpoint_attempt "${scratch_attempt}" "${final_attempt}" "${status}" \
    "${exit_code}" "${attempt_id}" "${local_attempt}"
  ATTEMPTED=$((ATTEMPTED + 1))
  if [ "${status}" = "success" ]; then
    SUCCESSES=$((SUCCESSES + 1))
  else
    FAILURES=$((FAILURES + 1))
  fi
  rm -rf -- "${scratch_attempt}"
  write_job_summary running
  echo "BOOTSTRAP_ATTEMPT_CHECKPOINTED job=${JOB_INDEX} attempt=${attempt_id} status=${status} successes=${SUCCESSES} failures=${FAILURES}"
done

if [ "${SUCCESSES}" -eq "${TARGET_SUCCESSES}" ]; then
  write_job_summary complete
  echo "BOOTSTRAP_BATCH_COMPLETE job=${JOB_INDEX} successes=${SUCCESSES} failures=${FAILURES} attempted=${ATTEMPTED}"
  exit 0
fi
write_job_summary exhausted
echo "BOOTSTRAP_BATCH_EXHAUSTED job=${JOB_INDEX} successes=${SUCCESSES} failures=${FAILURES} attempted=${ATTEMPTED}" >&2
exit 20
