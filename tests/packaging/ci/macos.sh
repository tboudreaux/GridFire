#!/usr/bin/env bash
# macOS (Apple Silicon) packaging checks, run on the machine itself:
#
#   tests/packaging/ci/macos.sh <platform name> <output dir>
#
# Builds a wheel as a release does (libomp built for MACOSX_DEPLOYMENT_TARGET, delocate with the fourdst libraries
# excluded), uses it from C++, C and Fortran (gf-config, pkg-config, CMake find_package) and from libgomp threads,
# then uses this source tree directly (CMake FetchContent, Meson subproject).
#
# Environment:
#   PYTHON          interpreter to build and test with (default: python3)
#   LIBOMP_LLVM_VERSION  LLVM release whose libomp is bundled (default: the release script's, see utils/wheels)
#   LIBOMP_PREFIX   that libomp, built for the deployment target; built here (into <output dir>) if missing
#   OMP_CC          a GCC with OpenMP for the libgomp thread check (default: the newest Homebrew gcc-NN)
#   FC              Fortran compiler (default: gfortran, else the newest Homebrew gfortran-NN)
#
# Results: <output dir>/results-<platform>.json, logs in <output dir>/logs, the wheel in <output dir>/wheel.
# Exits non-zero if any check failed.
set -uo pipefail

PLATFORM="${1:?usage: $0 <platform name> <output dir>}"
OUT="$(mkdir -p "${2:?usage: $0 <platform name> <output dir>}" && cd "$2" && pwd)"
HERE="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=lib.sh
source "${HERE}/lib.sh"
# shellcheck source=checks.sh
source "${HERE}/checks.sh"
ci_init "${PLATFORM}" "${OUT}"

PYTHON="${PYTHON:-python3}"
WORK="${OUT}/work"
mkdir -p "${WORK}"
export MACOSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET:-15.0}"
LLVM_VER="${LIBOMP_LLVM_VERSION:-19.1.7}"
LIBOMP_PREFIX="${LIBOMP_PREFIX:-${OUT}/libomp-${MACOSX_DEPLOYMENT_TARGET}-${LLVM_VER}}"

newest() { ls /opt/homebrew/bin/"$1"-[0-9]* 2>/dev/null | sort -V | tail -n1; }
OMP_CC="${OMP_CC:-$(newest gcc)}"
if [[ -z "${FC:-}" ]]; then
    if command -v gfortran > /dev/null; then FC=gfortran; else FC="$(newest gfortran)"; fi
fi
export FC CC=clang CXX=clang++
echo "python: $(${PYTHON} --version); gcc for libgomp: ${OMP_CC}; Fortran: ${FC}; libomp: ${LIBOMP_PREFIX}"

# The same libomp a release bundles: LLVM's, built for the deployment target (see utils/wheels)
build_libomp() {
    local src
    src="$(mktemp -d)"
    (
        set -e
        cd "${src}"
        curl -fsSLO "https://github.com/llvm/llvm-project/releases/download/llvmorg-${LLVM_VER}/openmp-${LLVM_VER}.src.tar.xz"
        curl -fsSLO "https://github.com/llvm/llvm-project/releases/download/llvmorg-${LLVM_VER}/cmake-${LLVM_VER}.src.tar.xz"
        tar xf "openmp-${LLVM_VER}.src.tar.xz"
        tar xf "cmake-${LLVM_VER}.src.tar.xz"
        mv "cmake-${LLVM_VER}.src" cmake
        cmake -S "openmp-${LLVM_VER}.src" -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
            -DCMAKE_OSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET}" -DCMAKE_OSX_ARCHITECTURES=arm64 \
            -DCMAKE_INSTALL_PREFIX="${LIBOMP_PREFIX}" -DLIBOMP_INSTALL_ALIASES=OFF
        cmake --build build
        cmake --install build
    )
    local rc=$?
    rm -rf "${src}"
    return ${rc}
}

"${PYTHON}" -m venv "${WORK}/buildenv"
"${WORK}/buildenv/bin/pip" install -q "meson-python>=0.19.0" "meson==1.9.1" ninja "pybind11==3.0.0" \
    "$(ci_fourdst_pin)" delocate "cmake>=3.24" || { echo "could not set up the build environment"; exit 1; }
export PATH="${WORK}/buildenv/bin:${PATH}"

if [[ ! -f "${LIBOMP_PREFIX}/lib/libomp.dylib" ]]; then
    build_libomp || { echo "could not build libomp"; exit 1; }
fi

build_wheel() {
    rm -rf "${WORK}/raw" "${OUT}/wheel"
    (
        cd "${CI_ROOT}" &&
        CPPFLAGS="-I${LIBOMP_PREFIX}/include" LDFLAGS="-L${LIBOMP_PREFIX}/lib" LIBRARY_PATH="${LIBOMP_PREFIX}/lib" \
            pip wheel . --no-deps --no-build-isolation \
            --config-settings=setup-args="-Dlibomp_prefix=${LIBOMP_PREFIX}" -w "${WORK}/raw"
    ) || return 1
    local fourdst_lib
    fourdst_lib="$(python -c 'import fourdst, os; print(os.pathsep.join(fourdst.get_lib_dirs()))')"
    DYLD_LIBRARY_PATH="${fourdst_lib}:${LIBOMP_PREFIX}/lib" delocate-wheel --require-archs arm64 \
        -e composition -e logging -e const -e reflect_cpp -w "${OUT}/wheel" -v "${WORK}"/raw/*.whl || return 1
    ci_assert_no_fourdst "${OUT}"/wheel/*.whl
}

ci_check wheel build_wheel
WHEEL=""
if ci_passed wheel; then WHEEL="$(ls "${OUT}"/wheel/*.whl | head -n1)"; fi
ci_wheel_checks "${PYTHON}" "${WHEEL}" "${WORK}/checks" "${OMP_CC}"
ci_check cmake-fetch check_cmake_fetch "${WORK}/cmake-fetch"
ci_check meson-subproject check_meson_subproject "${WORK}/meson-subproject"
ci_finish
