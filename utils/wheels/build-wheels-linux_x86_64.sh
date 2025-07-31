#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <git-repo-url>"
  exit 1
fi

REPO_URL="$1"
WORK_DIR="$(pwd)"
WHEEL_DIR="${WORK_DIR}/wheels_linux_x86_64"

echo "➤ Creating wheel output directory at ${WHEEL_DIR}"
mkdir -p "${WHEEL_DIR}"

TMPDIR="$(mktemp -d)"
echo "➤ Cloning ${REPO_URL} → ${TMPDIR}/project"
git clone "${REPO_URL}" "${TMPDIR}/project"

for IMAGE in \
  tboudreaux/manylinux_2_28_x86_64_boost_1_88_0:latest
do
  docker run --rm \
    -v "${WHEEL_DIR}":/io/wheels \
    -v "${TMPDIR}/project":/io/project \
    "${IMAGE}" \
    /bin/bash -eux -c '
      cd /io/project
      for PY in /opt/python/*/bin/python; do
        "$PY" -m pip install --upgrade pip setuptools wheel meson meson-python
        CC=clang CXX=clang++ "$PY" -m pip wheel . --config-settings=setup-args=-Dunity=on -w /io/wheels -vv
        auditwheel repair /io/wheels/*.whl -w /io/wheels
      done

      echo "✅ Linux wheels ready in /io/wheels"
    '
done