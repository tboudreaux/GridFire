#!/usr/bin/env bash
# Linux packaging checks, run inside GridFire's manylinux build image (the one the release wheels are built in):
#
#   docker run --rm -v "$PWD:/io/src:ro" -v "$PWD/ci-out:/io/out" tboudreaux/manylinux_2_28_<arch>_boost_1_88_0 \
#       bash /io/src/tests/packaging/ci/linux.sh <wheel|source> <platform name>
#
#   wheel:  build a cp312 wheel and repair it as a release does (auditwheel, fourdst libraries excluded), then use it
#           from C++, C and Fortran (gf-config, pkg-config, CMake find_package) and from libgomp threads
#   source: use this source tree directly: CMake FetchContent, and as a Meson subproject
#
# Results: /io/out/results-<platform>-<mode>.json (which checks passed), logs in /io/out/logs, the wheel in /io/out/wheel.
# Exits non-zero if any check failed.
set -uo pipefail

MODE="${1:?usage: $0 <wheel|source> <platform name>}"
PLATFORM="${2:?usage: $0 <wheel|source> <platform name>}"
OUT=/io/out

# Work on a copy: the mounted source is read-only, and Meson unpacks dependencies into subprojects/
rm -rf /tmp/src
cp -r /io/src /tmp/src
# shellcheck source=lib.sh
source /tmp/src/tests/packaging/ci/lib.sh
# shellcheck source=checks.sh
source /tmp/src/tests/packaging/ci/checks.sh
ci_init "${PLATFORM}" "${OUT}" "${MODE}"

PY=/opt/python/cp312-cp312/bin/python
"${PY}" -m venv /tmp/buildenv
/tmp/buildenv/bin/pip install -q "meson-python>=0.19.0" "meson>=1.9.1" ninja "pybind11==3.0.0" "$(ci_fourdst_pin)" \
    auditwheel "cmake>=3.24" || { echo "could not set up the build environment"; exit 1; }
export PATH=/tmp/buildenv/bin:${PATH}
export FC=gfortran

build_wheel() {
    rm -rf /tmp/raw "${OUT}/wheel"
    (cd /tmp/src && CC=clang CXX=clang++ pip wheel . --no-deps --no-build-isolation -w /tmp/raw) || return 1
    local fourdst_lib
    fourdst_lib="$(python -c 'import fourdst, os; print(os.pathsep.join(fourdst.get_lib_dirs()))')"
    LD_LIBRARY_PATH="${fourdst_lib}" auditwheel repair \
        --exclude "libcomposition.so*" --exclude "liblogging.so*" --exclude "libconst.so*" \
        --exclude "libreflect_cpp.so*" -w "${OUT}/wheel" /tmp/raw/*.whl || return 1
    ci_assert_no_fourdst "${OUT}"/wheel/*.whl
}

case "${MODE}" in
    wheel)
        ci_check wheel build_wheel
        WHEEL=""
        if ci_passed wheel; then WHEEL="$(ls "${OUT}"/wheel/*.whl | head -n1)"; fi
        export CC=gcc CXX=clang++   # consumers: a GCC host program (as DSEP's), C++23 from clang++
        ci_wheel_checks "${PY}" "${WHEEL}" /tmp/checks gcc
        ;;
    source)
        export CC=clang CXX=clang++
        ci_check cmake-fetch check_cmake_fetch /tmp/cmake-fetch
        ci_check meson-subproject check_meson_subproject /tmp/meson-subproject
        ;;
    *)
        echo "unknown mode ${MODE}; expected wheel or source"
        exit 2
        ;;
esac
ci_finish
