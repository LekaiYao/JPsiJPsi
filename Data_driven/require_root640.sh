#!/usr/bin/env bash

set -euo pipefail

JJ_REPO_ROOT="${JJ_REPO_ROOT:-/eos/home-l/leyao/26JJ/JPsiJPsi}"
JJ_DATA_DRIVEN_ROOT="${JJ_DATA_DRIVEN_ROOT:-${JJ_REPO_ROOT}/Data_driven}"
JJ_REQUIRED_ROOT_VERSION="6.40.02"

if ! command -v root >/dev/null 2>&1 || ! command -v root-config >/dev/null 2>&1; then
  echo "ROOT is not available in PATH" >&2
  return 65 2>/dev/null || exit 65
fi
if [ "$(root-config --version)" != "${JJ_REQUIRED_ROOT_VERSION}" ]; then
  echo "ROOT ${JJ_REQUIRED_ROOT_VERSION} is required; found $(root-config --version)" >&2
  return 65 2>/dev/null || exit 65
fi
export JJ_REPO_ROOT JJ_DATA_DRIVEN_ROOT JJ_REQUIRED_ROOT_VERSION
