#!/usr/bin/env bash
# Build and run the C++, C and Fortran consumers against an installed gridfire wheel, twice: once with the flags
# gf-config prints, once with pkg-config and the wheel's relocatable .pc files.
#
# Usage: tests/packaging/check_installed_wheel.sh <python with gridfire installed> [scratch dir] [route]
#   route: gf-config, pkg-config, or all (default)
# Needs cc and c++ (C++23); gfortran for the Fortran checks and pkg-config for the pkg-config ones (each skipped,
# with a message, if missing). Exits non-zero on the first failure.
set -euo pipefail

PY="${1:?usage: $0 <python> [scratch dir]}"
HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="${2:-$(mktemp -d)}"
ROUTE="${3:-all}"
mkdir -p "${WORK}"
CC="${CC:-cc}"
CXX="${CXX:-c++}"
FC="${FC:-gfortran}"
gf_config() { "${PY}" -c 'import sys, gridfire; sys.exit(gridfire.gf_config_main())' "$@"; }

echo "gridfire $(gf_config --version) from $(gf_config --lib-dir)"

# build_and_run <label> <cpp flags> <c flags> <fortran link flags> <fortran source>
build_and_run() {
    local label="$1" cxxflags="$2" cflags="$3" flibs="$4" fsrc="$5"
    echo "--- ${label}: C++"
    # shellcheck disable=SC2086
    "${CXX}" -std=c++23 -o "${WORK}/consumer_cpp" "${HERE}/consumer.cpp" ${cxxflags}
    "${WORK}/consumer_cpp"
    echo "--- ${label}: C"
    # shellcheck disable=SC2086
    "${CC}" -o "${WORK}/consumer_c" "${HERE}/consumer.c" ${cflags}
    "${WORK}/consumer_c"
    if command -v "${FC}" > /dev/null; then
        echo "--- ${label}: Fortran (module compiled from ${fsrc})"
        test -f "${fsrc}" || { echo "missing ${fsrc}"; exit 1; }
        (cd "${WORK}" && "${FC}" -c "${fsrc}" -o gridfire_mod.o)
        # shellcheck disable=SC2086
        "${FC}" -I"${WORK}" -o "${WORK}/consumer_f" "${HERE}/consumer.f90" "${WORK}/gridfire_mod.o" ${flibs}
        "${WORK}/consumer_f"
    else
        echo "--- ${label}: Fortran skipped (no ${FC})"
    fi
}

if [[ "${ROUTE}" == all || "${ROUTE}" == gf-config ]]; then
    build_and_run "gf-config" \
        "$(gf_config --cflags --libs | tr '\n' ' ')" \
        "$(gf_config --component c --cflags --libs | tr '\n' ' ')" \
        "$(gf_config --component fortran --libs)" \
        "$(gf_config --fortran-source)"
fi

if [[ "${ROUTE}" != all && "${ROUTE}" != pkg-config ]]; then
    :
elif command -v pkg-config > /dev/null; then
    export PKG_CONFIG_PATH
    PKG_CONFIG_PATH="$(gf_config --pkg-config-path)"
    for pc in gridfire gridfire_extern gridfire_fortran; do
        echo "${pc}.pc: version $(pkg-config --modversion "${pc}")"
    done
    build_and_run "pkg-config" \
        "$(pkg-config --cflags --libs gridfire)" \
        "$(pkg-config --cflags --libs gridfire_extern)" \
        "$(pkg-config --libs gridfire_fortran)" \
        "$(pkg-config --variable=fortran_source gridfire_fortran)"
elif [[ "${ROUTE}" == pkg-config ]]; then
    echo "--- pkg-config not found" >&2
    exit 1
else
    echo "--- pkg-config checks skipped (no pkg-config)"
fi
echo "packaging consumers passed"
