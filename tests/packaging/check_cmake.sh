#!/usr/bin/env bash
# Build and run the C++, C and Fortran consumers from CMake (tests/packaging/cmake_consumer), in either mode:
#
#   find_package against an installed gridfire wheel:
#       tests/packaging/check_cmake.sh find <python with gridfire installed> [scratch dir]
#   FetchContent, building this checkout with Meson underneath (the same path a fetch from git takes):
#       tests/packaging/check_cmake.sh fetch [scratch dir]
#
# Needs CMake >= 3.24, Ninja, C/C++23/Fortran compilers; fetch mode also needs Meson >= 1.5 and network access for
# GridFire's dependencies. Exits non-zero on failure.
set -euo pipefail

MODE="${1:?usage: $0 find <python> [scratch dir] | $0 fetch [scratch dir]}"
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "${HERE}/../.." && pwd)"

case "${MODE}" in
    find)
        PY="${2:?usage: $0 find <python> [scratch dir]}"
        WORK="${3:-$(mktemp -d)}"
        PREFIX="$("${PY}" -c 'import gridfire; print(gridfire.gf_get_cmake_prefix())')"
        echo "find_package(GridFire) with CMAKE_PREFIX_PATH=${PREFIX}"
        cmake -S "${HERE}/cmake_consumer" -B "${WORK}" -G Ninja -DCMAKE_PREFIX_PATH="${PREFIX}"
        ;;
    fetch)
        WORK="${2:-$(mktemp -d)}"
        echo "FetchContent of ${ROOT}"
        cmake -S "${HERE}/cmake_consumer" -B "${WORK}" -G Ninja -DGRIDFIRE_CONSUMER_FETCH=ON \
            -DFETCHCONTENT_TRY_FIND_PACKAGE_MODE=NEVER -DFETCHCONTENT_SOURCE_DIR_GRIDFIRE="${ROOT}"
        ;;
    *)
        echo "unknown mode ${MODE}; expected find or fetch" >&2
        exit 2
        ;;
esac
cmake --build "${WORK}"
ctest --test-dir "${WORK}" --output-on-failure
echo "CMake consumers (${MODE}) passed"
