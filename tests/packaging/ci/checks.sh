# Checks shared by the Linux and macOS packaging CI scripts (source after lib.sh). Each runs one way of using GridFire
# and is run through ci_check, so its outcome is recorded.

CI_HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CI_PACKAGING="$(cd "${CI_HERE}/.." && pwd)"
CI_ROOT="$(cd "${CI_PACKAGING}/../.." && pwd)"

# The fourdst version GridFire's wheel is built and used with (pinned in pyproject.toml)
ci_fourdst_pin() {
    grep -oE 'fourdst==[0-9][0-9a-zA-Z.]*' "${CI_ROOT}/pyproject.toml" | head -n1
}

# Fail if a repaired wheel carries its own copy of the fourdst libraries (it must use the fourdst wheel's)
ci_assert_no_fourdst() {
    if unzip -l "$1" | grep -E 'libcomposition|liblogging|libconst[^a-z]|libreflect_cpp'; then
        echo "ERROR: $1 contains vendored fourdst libraries"
        return 1
    fi
}

# A fresh environment, like a user's: fourdst from PyPI, then the wheel; GridFire imports and gf-config answers
check_python_install() {
    local py="$1" env="$2" wheel="$3"
    rm -rf "${env}"
    "${py}" -m venv "${env}"
    "${env}/bin/pip" install -q "$(ci_fourdst_pin)" "${wheel}"
    "${env}/bin/python" -c 'import gridfire; print("gridfire", gridfire.gf_version(), "from", gridfire.__file__)'
    "${env}/bin/gf-config" --version --cmake-prefix --pkg-config-path
}

# Contexts used concurrently from a host program built with GCC's OpenMP (libgomp), as DSEP built with
# gfortran -fopenmp would be, while GridFire carries its own LLVM libomp
check_threads_libgomp() {
    local env="$1" cc="$2" work="$3"
    mkdir -p "${work}"
    # shellcheck disable=SC2046
    "${cc}" -fopenmp -o "${work}/threads_gomp" "${CI_ROOT}/tests/extern/C/gridfire_threads.c" \
        $("${env}/bin/gf-config" --component c --cflags --libs | tr '\n' ' ') -lm
    local n
    for n in 8 16; do
        OMP_NUM_THREADS=${n} "${work}/threads_gomp"
    done
}

# The C++, C and Fortran consumers against the installed wheel, through gf-config or pkg-config
check_wheel_consumers() {
    local env="$1" work="$2" route="$3"
    bash "${CI_PACKAGING}/check_installed_wheel.sh" "${env}/bin/python" "${work}" "${route}"
}

check_cmake_find() {
    local env="$1" work="$2"
    rm -rf "${work}"
    bash "${CI_PACKAGING}/check_cmake.sh" find "${env}/bin/python" "${work}"
}

check_cmake_fetch() {
    rm -rf "$1"
    bash "${CI_PACKAGING}/check_cmake.sh" fetch "$1"
}

check_meson_subproject() {
    rm -rf "$1"
    bash "${CI_PACKAGING}/check_meson_subproject.sh" "$1"
}

# Run every check that needs the installed wheel, or record them all as skipped if there is no wheel
ci_wheel_checks() {
    local py="$1" wheel="$2" work="$3" omp_cc="$4"
    local names=(python-install gf-config pkg-config threads-libgomp cmake-find)
    if [[ -z "${wheel}" || ! -f "${wheel}" ]]; then
        local n
        for n in "${names[@]}"; do ci_skip "${n}" "no wheel"; done
        return
    fi
    local env="${work}/userenv"
    ci_check python-install check_python_install "${py}" "${env}" "${wheel}"
    if ! ci_passed python-install; then
        local n
        for n in "${names[@]:1}"; do ci_skip "${n}" "the wheel did not install"; done
        return
    fi
    ci_check gf-config check_wheel_consumers "${env}" "${work}/consumers-gf-config" gf-config
    ci_check pkg-config check_wheel_consumers "${env}" "${work}/consumers-pkg-config" pkg-config
    ci_check threads-libgomp check_threads_libgomp "${env}" "${omp_cc}" "${work}/threads"
    ci_check cmake-find check_cmake_find "${env}" "${work}/cmake-find"
}
