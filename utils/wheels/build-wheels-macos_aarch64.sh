#!/usr/bin/env bash
set -euo pipefail

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

for TOOL in pyenv cmake git; do
  if ! command -v "$TOOL" &> /dev/null; then
    echo "Error: ${TOOL} not found."
    exit 1
  fi
done

REPO_URL="$1"
LOCAL_FOURDST_WHEELS="${2:-}"
WORK_DIR="$(pwd)"
WHEEL_DIR="${WORK_DIR}/wheels_macos_aarch64_tmp"
FINAL_WHEEL_DIR="${WORK_DIR}/wheels_macos_aarch64"

echo "➤ Creating wheel output directories"
mkdir -p "${WHEEL_DIR}"
mkdir -p "${FINAL_WHEEL_DIR}"

export MACOSX_DEPLOYMENT_TARGET=15.0

LLVM_VER="17.0.6"
LIBOMP_PREFIX="${WORK_DIR}/libomp-${MACOSX_DEPLOYMENT_TARGET}-${LLVM_VER}"

if [[ ! -f "${LIBOMP_PREFIX}/lib/libomp.dylib" ]]; then
  echo "➤ Building libomp ${LLVM_VER} for macOS ${MACOSX_DEPLOYMENT_TARGET} → ${LIBOMP_PREFIX}"
  LIBOMP_SRC_DIR="$(mktemp -d)"
  pushd "${LIBOMP_SRC_DIR}" > /dev/null

  curl -fLO "https://github.com/llvm/llvm-project/releases/download/llvmorg-${LLVM_VER}/openmp-${LLVM_VER}.src.tar.xz"
  curl -fLO "https://github.com/llvm/llvm-project/releases/download/llvmorg-${LLVM_VER}/cmake-${LLVM_VER}.src.tar.xz"
  tar xf "openmp-${LLVM_VER}.src.tar.xz"
  tar xf "cmake-${LLVM_VER}.src.tar.xz"
  mv "cmake-${LLVM_VER}.src" cmake

  cmake -S "openmp-${LLVM_VER}.src" -B libomp-build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_OSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET}" \
    -DCMAKE_OSX_ARCHITECTURES=arm64 \
    -DCMAKE_INSTALL_PREFIX="${LIBOMP_PREFIX}" \
    -DLIBOMP_INSTALL_ALIASES=OFF
  cmake --build libomp-build -j "$(sysctl -n hw.ncpu)"
  cmake --install libomp-build

  popd > /dev/null
  rm -rf "${LIBOMP_SRC_DIR}"
fi

LIBOMP_MINOS="$(otool -l "${LIBOMP_PREFIX}/lib/libomp.dylib" | awk '/minos/ {print $2; exit}')"
if [[ "${LIBOMP_MINOS}" != "${MACOSX_DEPLOYMENT_TARGET}"* ]]; then
  echo "Error: cached libomp at ${LIBOMP_PREFIX} targets macOS ${LIBOMP_MINOS}," \
       "expected ${MACOSX_DEPLOYMENT_TARGET}. Delete the directory and re-run."
  exit 1
fi
echo "➤ Using libomp ${LLVM_VER} (min macOS ${LIBOMP_MINOS}) from ${LIBOMP_PREFIX}"

export CPPFLAGS="-I${LIBOMP_PREFIX}/include ${CPPFLAGS:-}"
export LDFLAGS="-L${LIBOMP_PREFIX}/lib ${LDFLAGS:-}"
export LIBRARY_PATH="${LIBOMP_PREFIX}/lib${LIBRARY_PATH:+:${LIBRARY_PATH}}"

TMPDIR="$(mktemp -d)"
echo "➤ Cloning ${REPO_URL} → ${TMPDIR}/project"
git clone --depth 1 "${REPO_URL}" "${TMPDIR}/project"
cd "${TMPDIR}/project"

FOURDST_PIN="$(grep -oE 'fourdst==[0-9][0-9a-zA-Z.]*' pyproject.toml | head -n1 || true)"
if [[ -n "${FOURDST_PIN}" ]]; then
  echo "➤ Project depends on ${FOURDST_PIN}; wheel repair will exclude fourdst libraries"
fi

PYTHON_VERSIONS=("3.9.23" "3.10.18" "3.11.13" "3.12.11" "3.13.5" "3.14.0rc1" "3.14.0rc1t")

eval "$(pyenv init -)"

for PY_VERSION in "${PYTHON_VERSIONS[@]}"; do
  (
    set -e

    pyenv shell "${PY_VERSION}"
    PY="$(pyenv which python)"

    echo "----------------------------------------------------------------"
    echo "➤ Building for $($PY --version) on macOS arm64"
    echo "----------------------------------------------------------------"

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
      "$PY" -m pip wheel . --no-deps --no-build-isolation \
        --config-settings=setup-args="-Dlibomp_prefix=${LIBOMP_PREFIX}" \
        -w "${WHEEL_DIR}" -v
    CURRENT_WHEEL=$(find "${WHEEL_DIR}" -name "*.whl" | head -n 1)

    echo "➤ Repairing wheel with delocate"
    DELOCATE_DYLD_PATH="${LIBOMP_PREFIX}/lib"
    if [[ -n "${FOURDST_PIN}" ]]; then
      FOURDST_LIB_PATH="$("$PY" -c 'import fourdst, os; print(os.pathsep.join(fourdst.get_lib_dirs()))')"
      DELOCATE_DYLD_PATH="${FOURDST_LIB_PATH}:${DELOCATE_DYLD_PATH}"
      DYLD_LIBRARY_PATH="${DELOCATE_DYLD_PATH}" \
        delocate-wheel --require-archs arm64 \
          -e composition -e logging -e const -e reflect_cpp \
          -w "${FINAL_WHEEL_DIR}" -v "$CURRENT_WHEEL"
    else
      DYLD_LIBRARY_PATH="${DELOCATE_DYLD_PATH}" \
        delocate-wheel --require-archs arm64 \
          -w "${FINAL_WHEEL_DIR}" -v "$CURRENT_WHEEL"
    fi

    REPAIRED_WHEEL="${FINAL_WHEEL_DIR}/$(basename "$CURRENT_WHEEL")"

    if [[ -n "${FOURDST_PIN}" ]]; then
      if unzip -l "${REPAIRED_WHEEL}" | grep -E 'libcomposition|liblogging|libconst|libreflect_cpp'; then
        echo "ERROR: repaired wheel contains vendored fourdst libraries"
        exit 1
      fi
    fi
    if unzip -l "${REPAIRED_WHEEL}" | grep -q 'libomp.dylib'; then
      CHECK_DIR="$(mktemp -d)"
      unzip -q "${REPAIRED_WHEEL}" '*/libomp.dylib' -d "${CHECK_DIR}" || true
      BUNDLED_OMP="$(find "${CHECK_DIR}" -name 'libomp.dylib' | head -n1)"
      BUNDLED_MINOS="$(otool -l "${BUNDLED_OMP}" | awk '/minos/ {print $2; exit}')"
      rm -rf "${CHECK_DIR}"
      if [[ "${BUNDLED_MINOS}" != "${MACOSX_DEPLOYMENT_TARGET}"* ]]; then
        echo "ERROR: bundled libomp targets macOS ${BUNDLED_MINOS}, expected ${MACOSX_DEPLOYMENT_TARGET}"
        exit 1
      fi
      echo "➤ Bundled libomp targets macOS ${BUNDLED_MINOS} ✓"
    fi

    rm "$CURRENT_WHEEL"
  )
done

# Cleanup
rm -rf "${TMPDIR}"
rm -rf "${WHEEL_DIR}"

echo "✅ All builds complete. Artifacts in ${FINAL_WHEEL_DIR}"