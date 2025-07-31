#!/usr/bin/env bash
set -euo pipefail

if [[ $(uname -m) != "arm64" ]]; then
  echo "Error: This script is intended to run on an Apple Silicon (arm64) Mac."
  exit 1
fi

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <git-repo-url>"
  exit 1
fi

# --- Initial Setup ---
REPO_URL="$1"
WORK_DIR="$(pwd)"
WHEEL_DIR="${WORK_DIR}/wheels_macos_aarch64_tmp"
FINAL_WHEEL_DIR="${WORK_DIR}/wheels_macos_aarch64"

echo "➤ Creating wheel output directories"
mkdir -p "${WHEEL_DIR}"
mkdir -p "${FINAL_WHEEL_DIR}"

TMPDIR="$(mktemp -d)"
echo "➤ Cloning ${REPO_URL} → ${TMPDIR}/project"
git clone --depth 1 "${REPO_URL}" "${TMPDIR}/project"
cd "${TMPDIR}/project"

# --- macOS Build Configuration ---
export MACOSX_DEPLOYMENT_TARGET=12.0

PYTHON_VERSIONS=("3.8.20" "3.9.23" "3.10.18" "3.11.13" "3.12.11" "3.13.5" "3.13.5t" "3.14.0rc1" "3.14.0rc1t" 'pypy3.10-7.3.19' "pypy3.11-7.3.20")

if ! command -v pyenv &> /dev/null; then
    echo "Error: pyenv not found. Please install it to manage Python versions."
    exit 1
fi
eval "$(pyenv init -)"

for PY_VERSION in "${PYTHON_VERSIONS[@]}"; do
  (
    set -e

    if ! pyenv versions --bare --filter="${PY_VERSION}." &>/dev/null; then
        echo "⚠️  Python version matching '${PY_VERSION}.*' not found by pyenv. Skipping."
        continue
    fi
    
    pyenv shell "${PY_VERSION}"

    PY="$(pyenv which python)"
    echo "➤ Building for $($PY --version) on macOS arm64 (target: ${MACOSX_DEPLOYMENT_TARGET})"

    "$PY" -m pip install --upgrade pip setuptools wheel meson meson-python delocate

    CC=clang CXX=clang++ "$PY" -m pip wheel . \
      --config-settings=setup-args=-Dunity=on \
      -w "${WHEEL_DIR}" -vv

    echo "➤ Repairing wheel(s) with delocate"
    delocate-wheel -w "${FINAL_WHEEL_DIR}" "${WHEEL_DIR}"/*.whl

    rm "${WHEEL_DIR}"/*.whl

  )
done

rm -rf "${TMPDIR}"
rm -rf "${WHEEL_DIR}"

