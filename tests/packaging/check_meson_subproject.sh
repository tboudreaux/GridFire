#!/usr/bin/env bash
# Build tests/packaging/meson_consumer with this checkout of GridFire as its subproject, then run its C++, C and
# Fortran tests. The subproject is a clean copy of the repository's files (as a clone of this commit would be: no
# unpacked dependencies), so Meson fetches GridFire's own dependencies into the consumer's subprojects/ directory,
# exactly as for a user whose project ships tests/packaging/meson_consumer/gridfire.wrap.
#
# Usage: tests/packaging/check_meson_subproject.sh [scratch dir]
# Needs git, meson, ninja, network access, C/C++ compilers (C++23) and gfortran. Builds everything from scratch.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "${HERE}/../.." && pwd)"
WORK="${1:-$(mktemp -d)}"
PROJECT="${WORK}/meson_consumer"

rm -rf "${PROJECT}"
mkdir -p "${PROJECT}/subprojects"
cp "${HERE}/meson_consumer/meson.build" "${HERE}/consumer.c" "${HERE}/consumer.cpp" "${HERE}/consumer.f90" "${PROJECT}/"
cp "${HERE}/meson_consumer/gridfire.wrap" "${PROJECT}/subprojects/"
# Tracked files, plus files not yet committed under the directories the build uses (not ignored ones)
mkdir -p "${PROJECT}/subprojects/gridfire"   # present, so meson uses it instead of cloning
(
    cd "${ROOT}"
    { git ls-files; git ls-files --others --exclude-standard -- src build-check build-config build-extra build-python \
        tests subprojects stubs; } | sort -u | tar -cf - -T -
) | tar -xf - -C "${PROJECT}/subprojects/gridfire"

meson setup "${PROJECT}/build" "${PROJECT}" --buildtype=release
meson compile -C "${PROJECT}/build"
meson test -C "${PROJECT}/build" --print-errorlogs
echo "meson subproject consumer passed"
