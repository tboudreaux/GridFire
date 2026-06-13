#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <git-repo-url> [fourdst-wheels-dir]"
  echo "  fourdst-wheels-dir: optional local directory of fourdst wheels to"
  echo "  install from instead of PyPI (for bootstrapping a new fourdst+gridfire pair)"
  exit 1
fi

REPO_URL="$1"
LOCAL_FOURDST_WHEELS="${2:-}"
WORK_DIR="$(pwd)"
WHEEL_DIR="${WORK_DIR}/wheels_linux_aarch64"

echo "➤ Creating wheel output directory at ${WHEEL_DIR}"
mkdir -p "${WHEEL_DIR}"

TMPDIR="$(mktemp -d)"
echo "➤ Cloning ${REPO_URL} → ${TMPDIR}/project"
git clone "${REPO_URL}" "${TMPDIR}/project"

DOCKER_MOUNTS=(-v "${WHEEL_DIR}":/io/wheels -v "${TMPDIR}/project":/io/project)
if [[ -n "${LOCAL_FOURDST_WHEELS}" ]]; then
  DOCKER_MOUNTS+=(-v "${LOCAL_FOURDST_WHEELS}":/io/fourdst-wheels)
fi

for IMAGE in \
  tboudreaux/manylinux_2_28_aarch64_boost_1_88_0:latest
do
  docker run --rm \
    "${DOCKER_MOUNTS[@]}" \
    "${IMAGE}" \
    /bin/bash -uxo pipefail -c '
      cd /io/project

      # ----------------------------------------------------------------
      # Project identity: package name from pyproject.toml, version from
      # meson (the version pip will stamp into the wheel). Used both for
      # the skip-if-already-built check and the post-repair checks, so a
      # stale wheel from an OLDER project version never causes a skip.
      # ----------------------------------------------------------------
      PKG="$(sed -n "s/^name *= *\"\(.*\)\"/\1/p" pyproject.toml | head -n1)"
      PKG="${PKG//-/_}"   # wheel filename normalization

      BOOT_PY=/opt/python/cp312-cp312/bin/python
      "$BOOT_PY" -m pip install --quiet meson
      VERSION="$("$BOOT_PY" -c "
import json, subprocess, sys
out = subprocess.check_output(
    [sys.executable, \"-m\", \"mesonbuild.mesonmain\", \"introspect\",
     \"meson.build\", \"--projectinfo\"])
print(json.loads(out)[\"version\"])
" 2>/dev/null || true)"
      if [ -z "$VERSION" ]; then
        # fallback: literal version in project()
        VERSION="$(grep -oE "version *: *.[0-9][0-9a-zA-Z.+-]*" meson.build | head -n1 | grep -oE "[0-9][0-9a-zA-Z.+-]*" || true)"
      fi
      if [ -z "$VERSION" ]; then
        echo "ERROR: could not determine project version; refusing to guess for skip logic"
        exit 1
      fi
      echo "➤ Building ${PKG} ${VERSION}"

      # Does this project link against the fourdst wheel? Single source of
      # truth: the pin in pyproject.toml.
      FOURDST_PIN="$(grep -oE "fourdst==[0-9][0-9a-zA-Z.]*" pyproject.toml | head -n1 || true)"

      # If a local fourdst wheel dir was mounted, let pip (including the
      # isolated build env) resolve fourdst from it.
      if [ -d /io/fourdst-wheels ]; then
        export PIP_FIND_LINKS=/io/fourdst-wheels
      fi

      build_one() {
        # Runs the full build+repair for one interpreter. Returns nonzero
        # on any failure; never exits the whole script (errexit is off in
        # the caller around this function).
        set -e
        local PY="$1" PYTAG="$2"

        "$PY" -m pip install --upgrade pip setuptools wheel meson meson-python

        # Build into a per-iteration temp dir so we repair exactly the
        # wheel we just built.
        local BUILD_WHEEL_DIR
        BUILD_WHEEL_DIR="$(mktemp -d)"
        CC=clang CXX=clang++ "$PY" -m pip wheel . --no-deps \
          -w "$BUILD_WHEEL_DIR" -vv

        local CURRENT_WHEEL
        CURRENT_WHEEL="$(find "$BUILD_WHEEL_DIR" -name "*.whl" | head -n1)"

        if [ -n "$FOURDST_PIN" ]; then
          # Install fourdst for THIS interpreter so auditwheel can resolve
          # the libraries it must NOT graft. Excluding them keeps fourdst a
          # runtime dependency: grafting copies would break cross-package
          # pybind11 type compatibility.
          "$PY" -m pip install --force-reinstall "$FOURDST_PIN"
          local FOURDST_LIB_PATH
          FOURDST_LIB_PATH="$("$PY" -c "import fourdst, os; print(os.pathsep.join(fourdst.get_lib_dirs()))")"
          LD_LIBRARY_PATH="$FOURDST_LIB_PATH" auditwheel repair \
            --exclude "libcomposition.so*" \
            --exclude "liblogging.so*" \
            --exclude "libconst.so*" \
            --exclude "libreflect_cpp.so*" \
            -w /io/wheels "$CURRENT_WHEEL"

          # Post-repair sanity check on the wheel we just produced
          local REPAIRED
          REPAIRED="$(find /io/wheels -name "${PKG}-${VERSION}-${PYTAG}-*manylinux*.whl" | head -n1)"
          if [ -z "$REPAIRED" ]; then
            echo "ERROR: repaired wheel for ${PYTAG} not found after auditwheel"
            return 1
          fi
          if unzip -l "$REPAIRED" | grep -E "libcomposition|liblogging|libconst[^a-z]|libreflect_cpp"; then
            echo "ERROR: repaired wheel contains vendored fourdst libraries"
            rm -f "$REPAIRED"   # do not leave a poisoned wheel that would be skipped next run
            return 1
          fi
        else
          auditwheel repair -w /io/wheels "$CURRENT_WHEEL"
        fi

        rm -rf "$BUILD_WHEEL_DIR"
      }

      FAILED_TAGS=""
      SKIPPED_TAGS=""
      BUILT_TAGS=""

      for PY in /opt/python/*/bin/python; do
        # /opt/python/<pythontag>-<abitag>/bin/python — the directory name
        # is exactly the {python tag}-{abi tag} pair used in wheel filenames
        PYTAG="$(basename "$(dirname "$(dirname "$PY")")")"

        # ------------------------------------------------------------
        # 1. Skip if a repaired wheel for THIS name+version+interpreter
        #    already exists (a wheel from an older version will not match
        #    because VERSION is part of the pattern).
        # ------------------------------------------------------------
        if compgen -G "/io/wheels/${PKG}-${VERSION}-${PYTAG}-*manylinux*.whl" > /dev/null; then
          echo "➤ ${PYTAG}: wheel for ${PKG} ${VERSION} already present — skipping"
          SKIPPED_TAGS="${SKIPPED_TAGS} ${PYTAG}"
          continue
        fi

        # ------------------------------------------------------------
        # 2. Build; on failure, record and continue with the next python
        # ------------------------------------------------------------
        echo "================================================================"
        echo "➤ ${PYTAG}: building ${PKG} ${VERSION}"
        echo "================================================================"
        if ( build_one "$PY" "$PYTAG" ); then
          BUILT_TAGS="${BUILT_TAGS} ${PYTAG}"
        else
          echo "✗ ${PYTAG}: BUILD FAILED — continuing with remaining versions"
          FAILED_TAGS="${FAILED_TAGS} ${PYTAG}"
        fi
      done

      echo "================================================================"
      echo "Summary for ${PKG} ${VERSION}:"
      echo "  built:  ${BUILT_TAGS:- none}"
      echo "  skipped:${SKIPPED_TAGS:- none}"
      echo "  failed: ${FAILED_TAGS:- none}"
      echo "================================================================"

      if [ -n "$FAILED_TAGS" ]; then
        echo "✗ Some builds failed:${FAILED_TAGS}"
        exit 1
      fi
      echo "✅ Linux wheels ready in /io/wheels"
    '
done