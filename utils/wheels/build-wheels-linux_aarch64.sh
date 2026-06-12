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
    /bin/bash -eux -c '
      cd /io/project

      FOURDST_PIN="$(grep -oE "fourdst==[0-9][0-9a-zA-Z.]*" pyproject.toml | head -n1 || true)"

      if [ -d /io/fourdst-wheels ]; then
        export PIP_FIND_LINKS=/io/fourdst-wheels
      fi

      for PY in /opt/python/*/bin/python; do
        "$PY" -m pip install --upgrade pip setuptools wheel meson meson-python

        BUILD_WHEEL_DIR="$(mktemp -d)"
        CC=clang CXX=clang++ "$PY" -m pip wheel . \
          --no-deps --config-settings=setup-args=-Dunity=on \
          -w "$BUILD_WHEEL_DIR" -vv

        CURRENT_WHEEL="$(find "$BUILD_WHEEL_DIR" -name "*.whl" | head -n1)"

        if [ -n "$FOURDST_PIN" ]; then
          "$PY" -m pip install --force-reinstall "$FOURDST_PIN"
          FOURDST_LIB_PATH="$("$PY" -c "import fourdst, os; print(os.pathsep.join(fourdst.get_lib_dirs()))")"
          LD_LIBRARY_PATH="$FOURDST_LIB_PATH" auditwheel repair \
            --exclude "libcomposition.so*" \
            --exclude "liblogging.so*" \
            --exclude "libconst.so*" \
            --exclude "libreflect_cpp.so*" \
            -w /io/wheels "$CURRENT_WHEEL"

          REPAIRED="$(ls -t /io/wheels/*.whl | head -n1)"
          if unzip -l "$REPAIRED" | grep -E "libcomposition|liblogging|libconst[^a-z]|libreflect_cpp"; then
            echo "ERROR: repaired wheel contains vendored fourdst libraries"
            exit 1
          fi
        else
          auditwheel repair -w /io/wheels "$CURRENT_WHEEL"
        fi

        rm -rf "$BUILD_WHEEL_DIR"
      done

      echo "Linux wheels ready in /io/wheels"
    '
done
