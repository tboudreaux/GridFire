#!/usr/bin/env bash
set -euo pipefail

# 1. Validation
if [[ $(uname -m) != "arm64" ]]; then
  echo "Error: This script is intended to run on an Apple Silicon (arm64) Mac."
  exit 1
fi

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <git-repo-url> [fourdst-wheels-dir]"
  echo "  fourdst-wheels-dir: optional local directory of fourdst wheels to"
  echo "  install from instead of PyPI (for bootstrapping a new fourdst+gridfire pair)"
  exit 1
fi

# 2. Setup Directories
REPO_URL="$1"
LOCAL_FOURDST_WHEELS="${2:-}"
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

# 3. Build Configuration
# NOTE: must match the value used to build the fourdst wheels — the two
# packages are one ABI unit.
export MACOSX_DEPLOYMENT_TARGET=15.0

# Does this project link against the fourdst wheel? Derive the pin from
# pyproject.toml so there is a single source of truth.
FOURDST_PIN="$(grep -oE 'fourdst==[0-9][0-9a-zA-Z.]*' pyproject.toml | head -n1 || true)"
if [[ -n "${FOURDST_PIN}" ]]; then
  echo "➤ Project depends on ${FOURDST_PIN}; wheel repair will exclude fourdst libraries"
fi

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

    # Install build deps explicitly so we can skip build isolation.
    # IMPORTANT: with --no-build-isolation, EVERYTHING in
    # build-system.requires must be installed here by hand — including
    # fourdst, otherwise the meson probe (`import fourdst`) fails.
    "$PY" -m pip install --upgrade pip setuptools wheel meson-python delocate
    "$PY" -m pip install meson==1.9.1

    if [[ -n "${FOURDST_PIN}" ]]; then
      if [[ -n "${LOCAL_FOURDST_WHEELS}" ]]; then
        "$PY" -m pip install --force-reinstall \
          --find-links "${LOCAL_FOURDST_WHEELS}" "${FOURDST_PIN}"
      else
        "$PY" -m pip install --force-reinstall "${FOURDST_PIN}"
      fi
    fi

    echo "➤ Building wheel with ccache enabled"
    echo "➤ Found meson version $(meson --version)"

    CC="ccache clang" CXX="ccache clang++" \
      "$PY" -m pip wheel . --no-build-isolation -w "${WHEEL_DIR}" -v

    # We expect exactly one new wheel in the tmp dir per iteration
    CURRENT_WHEEL=$(find "${WHEEL_DIR}" -name "*.whl" | head -n 1)

    echo "➤ Repairing wheel with delocate"
    if [[ -n "${FOURDST_PIN}" ]]; then
      # Resolve @rpath references against the installed fourdst wheel,
      # but EXCLUDE its libraries from being grafted into this wheel:
      # they must stay a runtime dependency, or cross-package pybind11
      # type compatibility breaks.
      FOURDST_LIB_PATH="$("$PY" -c 'import fourdst, os; print(os.pathsep.join(fourdst.get_lib_dirs()))')"
      DYLD_LIBRARY_PATH="${FOURDST_LIB_PATH}" \
        delocate-wheel --require-archs arm64 \
          -e composition -e logging -e const -e reflect_cpp \
          -w "${FINAL_WHEEL_DIR}" -v "$CURRENT_WHEEL"
    else
      delocate-wheel --require-archs arm64 -w "${FINAL_WHEEL_DIR}" -v "$CURRENT_WHEEL"
    fi

    # Post-repair sanity check: import the wheel in a throwaway env and
    # make sure no fourdst library snuck back in.
    if [[ -n "${FOURDST_PIN}" ]]; then
      REPAIRED_WHEEL=$(find "${FINAL_WHEEL_DIR}" -name "*.whl" -newer "$CURRENT_WHEEL" | head -n 1)
      if [[ -n "${REPAIRED_WHEEL}" ]] && unzip -l "${REPAIRED_WHEEL}" | grep -E 'libcomposition|liblogging|libconst|libreflect_cpp' ; then
        echo "ERROR: repaired wheel contains vendored fourdst libraries"
        exit 1
      fi
    fi

    # Clean up the intermediate wheel from this iteration so it doesn't confuse the next
    rm "$CURRENT_WHEEL"
  )
done

# Cleanup
rm -rf "${TMPDIR}"
rm -rf "${WHEEL_DIR}"

echo "✅ All builds complete. Artifacts in ${FINAL_WHEEL_DIR}"