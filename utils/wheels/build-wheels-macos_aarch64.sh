#!/usr/bin/env bash
set -euo pipefail

# 1. Validation
if [[ $(uname -m) != "arm64" ]]; then
  echo "Error: This script is intended to run on an Apple Silicon (arm64) Mac."
  exit 1
fi

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <git-repo-url>"
  exit 1
fi

# 2. Setup Directories
REPO_URL="$1"
WORK_DIR="$(pwd)"
WHEEL_DIR="${WORK_DIR}/wheels_macos_aarch64_tmp"
FINAL_WHEEL_DIR="${WORK_DIR}/wheels_macos_aarch64"
RPATH_SCRIPT="${WORK_DIR}/../../build-python/fix_rpaths.py" # Assumes script is in this location relative to execution

echo "➤ Creating wheel output directories"
mkdir -p "${WHEEL_DIR}"
mkdir -p "${FINAL_WHEEL_DIR}"

TMPDIR="$(mktemp -d)"
echo "➤ Cloning ${REPO_URL} → ${TMPDIR}/project"
git clone --depth 1 "${REPO_URL}" "${TMPDIR}/project"
cd "${TMPDIR}/project"

# 3. Build Configuration
export MACOSX_DEPLOYMENT_TARGET=15.0

PYTHON_VERSIONS=("3.9.23" "3.10.18" "3.11.13" "3.12.11" "3.13.5" "3.13.5t" "3.14.0rc1" "3.14.0rc1t" 'pypy3.10-7.3.19' "pypy3.11-7.3.20")

if ! command -v pyenv &> /dev/null; then
    echo "Error: pyenv not found. Please install it to manage Python versions."
    exit 1
fi
eval "$(pyenv init -)"

# 4. Build Loop
for PY_VERSION in "${PYTHON_VERSIONS[@]}"; do
  (
    set -e

    pyenv shell "${PY_VERSION}"
    PY="$(pyenv which python)"
    
    echo "----------------------------------------------------------------"
    echo "➤ Building for $($PY --version) on macOS arm64"
    echo "----------------------------------------------------------------"

    # Install build deps explicitly so we can skip build isolation
    "$PY" -m pip install --upgrade pip setuptools wheel meson-python delocate
    "$PY" -m pip install meson==1.9.1

    echo "➤ Building wheel with ccache enabled"
    echo "➤ Found meson version $(meson --version)"

    # for every single build, saving significant I/O and network time.
    CC="ccache clang" CXX="ccache clang++" "$PY" -m pip wheel . --no-build-isolation -w "${WHEEL_DIR}" -v

    # We expect exactly one new wheel in the tmp dir per iteration
    CURRENT_WHEEL=$(find "${WHEEL_DIR}" -name "*.whl" | head -n 1)

    echo "➤ Repairing wheel with delocate"
    # Delocate moves the repaired wheel to FINAL_WHEEL_DIR
    delocate-wheel -w "${FINAL_WHEEL_DIR}" "$CURRENT_WHEEL"

    # Clean up the intermediate wheel from this iteration so it doesn't confuse the next
    rm "$CURRENT_WHEEL"
  )
done

# Cleanup
rm -rf "${TMPDIR}"
rm -rf "${WHEEL_DIR}"

echo "✅ All builds complete. Artifacts in ${FINAL_WHEEL_DIR}"
