#!/usr/bin/env bash
#
# install.sh -- GridFire build-dependency checker.
#
# This script does NOT install anything. It detects your operating system,
# checks each build dependency, and prints the exact commands to run if
# something is missing or too old. Nothing on your system is modified.
#
# Dependencies checked:
#   1. Meson      >= 1.6.0
#   2. Ninja      (any recent)
#   3. CMake      >= 3.16
#   4. pkg-config (or pkgconf)
#   5. A C and C++ toolchain: GCC >= 14  OR  Clang >= 17  (either is fine),
#      with a working C++23 standard library
#   6. Fortran (gfortran >= 14) -- only when --fortran is passed
#
# Supported platforms: Linux (apt / dnf / pacman / zypper) and macOS on
# Apple Silicon. macOS x86_64 and Windows are NOT supported.
#
# Usage:
#   ./install.sh            check core dependencies
#   ./install.sh --fortran  also check for a Fortran compiler
#   ./install.sh --help     show this help

set -u

MESON_MIN="1.6.0"
CMAKE_MIN="3.16"
GCC_MIN="14"
CLANG_MIN="17"
GFORTRAN_MIN="14"

if [ -t 1 ]; then
  RED=$'\033[0;31m'; GREEN=$'\033[0;32m'; YELLOW=$'\033[1;33m'
  BLUE=$'\033[0;34m'; BOLD=$'\033[1m'; NC=$'\033[0m'
else
  RED=""; GREEN=""; YELLOW=""; BLUE=""; BOLD=""; NC=""
fi

ok()    { printf '  %s[ OK ]%s %s\n' "$GREEN" "$NC" "$1"; }
bad()   { printf '  %s[FAIL]%s %s\n' "$RED" "$NC" "$1"; }
warn()  { printf '  %s[WARN]%s %s\n' "$YELLOW" "$NC" "$1"; }
info()  { printf '  %s[INFO]%s %s\n' "$BLUE" "$NC" "$1"; }
hdr()   { printf '\n%s%s%s\n' "$BOLD" "$1" "$NC"; }

FIX_NOTES=""
FAIL_COUNT=0
add_fix() { FIX_NOTES="${FIX_NOTES}$1\n"; }
note_fail() { FAIL_COUNT=$((FAIL_COUNT + 1)); }

WANT_FORTRAN=0
for arg in "$@"; do
  case "$arg" in
    --fortran) WANT_FORTRAN=1 ;;
    -h|--help) sed -n '2,/^$/p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) printf '%sUnknown argument: %s%s\n' "$RED" "$arg" "$NC"; exit 2 ;;
  esac
done

version_ge() {
  local a b ai bi n i
  IFS='.' read -ra a <<< "$1"
  IFS='.' read -ra b <<< "$2"
  n=${#a[@]}; [ ${#b[@]} -gt "$n" ] && n=${#b[@]}
  for ((i=0; i<n; i++)); do
    ai=${a[i]:-0}; bi=${b[i]:-0}
    ai=${ai%%[!0-9]*}; bi=${bi%%[!0-9]*}
    ai=${ai:-0}; bi=${bi:-0}
    if ((10#$ai > 10#$bi)); then return 0; fi
    if ((10#$ai < 10#$bi)); then return 1; fi
  done
  return 0
}

extract_version() {
  printf '%s' "$1" | grep -oE '[0-9]+(\.[0-9]+){0,3}' | head -n1
}

OS="$(uname -s)"
ARCH="$(uname -m)"
PKG_MGR=""
DISTRO_NAME=""

detect_platform() {
  case "$OS" in
    Darwin)
      if [ "$ARCH" != "arm64" ]; then
        printf '%sUnsupported platform.%s\n' "$RED" "$NC"
        printf 'GridFire supports macOS on Apple Silicon (arm64) only; detected macOS on %s.\n' "$ARCH"
        printf 'macOS x86_64 (Intel) is not supported.\n'
        exit 1
      fi
      DISTRO_NAME="macOS (Apple Silicon)"
      if command -v brew >/dev/null 2>&1; then
        PKG_MGR="brew"
      else
        PKG_MGR="brew-missing"
      fi
      ;;
    Linux)
      if command -v apt-get >/dev/null 2>&1; then PKG_MGR="apt"
      elif command -v dnf >/dev/null 2>&1; then PKG_MGR="dnf"
      elif command -v pacman >/dev/null 2>&1; then PKG_MGR="pacman"
      elif command -v zypper >/dev/null 2>&1; then PKG_MGR="zypper"
      else PKG_MGR="unknown"; fi
      if [ -r /etc/os-release ]; then
        # shellcheck disable=SC1091
        DISTRO_NAME="$(. /etc/os-release; printf '%s' "${PRETTY_NAME:-Linux}")"
      else
        DISTRO_NAME="Linux"
      fi
      ;;
    *)
      printf '%sUnsupported platform: %s.%s\n' "$RED" "$OS" "$NC"
      printf 'GridFire supports Linux and macOS (Apple Silicon). Windows is not supported.\n'
      printf 'On Windows, use WSL2 with a supported Linux distribution.\n'
      exit 1
      ;;
  esac
}

pkg_cmd() {
  case "$PKG_MGR" in
    apt)
      case "$1" in
        ninja)     echo "sudo apt install ninja-build" ;;
        cmake)     echo "sudo apt install cmake" ;;
        pkgconfig) echo "sudo apt install pkg-config" ;;
        gcc)       echo "sudo apt install g++-14 gcc-14" ;;
        clang)     echo "sudo apt install clang-18 libstdc++-14-dev" ;;
        gfortran)  echo "sudo apt install gfortran-14" ;;
      esac ;;
    dnf)
      case "$1" in
        ninja)     echo "sudo dnf install ninja-build" ;;
        cmake)     echo "sudo dnf install cmake" ;;
        pkgconfig) echo "sudo dnf install pkgconf-pkg-config" ;;
        gcc)       echo "sudo dnf install gcc-c++" ;;
        clang)     echo "sudo dnf install clang" ;;
        gfortran)  echo "sudo dnf install gcc-gfortran" ;;
      esac ;;
    pacman)
      case "$1" in
        ninja)     echo "sudo pacman -S ninja" ;;
        cmake)     echo "sudo pacman -S cmake" ;;
        pkgconfig) echo "sudo pacman -S pkgconf" ;;
        gcc)       echo "sudo pacman -S gcc" ;;
        clang)     echo "sudo pacman -S clang" ;;
        gfortran)  echo "sudo pacman -S gcc-fortran" ;;
      esac ;;
    zypper)
      case "$1" in
        ninja)     echo "sudo zypper install ninja" ;;
        cmake)     echo "sudo zypper install cmake" ;;
        pkgconfig) echo "sudo zypper install pkg-config" ;;
        gcc)       echo "sudo zypper install gcc14-c++" ;;
        clang)     echo "sudo zypper install clang" ;;
        gfortran)  echo "sudo zypper install gcc14-fortran" ;;
      esac ;;
    brew|brew-missing)
      case "$1" in
        ninja)     echo "brew install ninja" ;;
        cmake)     echo "brew install cmake" ;;
        pkgconfig) echo "brew install pkg-config" ;;
        gcc)       echo "brew install gcc" ;;
        clang)     echo "brew install llvm   # then use /opt/homebrew/opt/llvm/bin/clang++" ;;
        gfortran)  echo "brew install gcc    # provides gfortran" ;;
      esac ;;
    *)
      echo "(install '$1' using your system package manager)" ;;
  esac
}

check_meson() {
  if command -v meson >/dev/null 2>&1; then
    local v; v="$(extract_version "$(meson --version 2>/dev/null)")"
    if [ -n "$v" ] && version_ge "$v" "$MESON_MIN"; then
      ok "Meson $v (>= $MESON_MIN)"
    else
      bad "Meson ${v:-?} found but GridFire needs >= $MESON_MIN"
      note_fail
      add_fix "Meson upgrade (isolated, recommended):\n    pipx upgrade meson || pipx install meson\n  or with pip:\n    python3 -m pip install --user --upgrade meson"
    fi
  else
    bad "Meson not found (need >= $MESON_MIN)"
    note_fail
    add_fix "Install Meson (isolated, recommended):\n    pipx install meson\n  If pipx is unavailable:\n    python3 -m pip install --user meson\n  (pipx itself: ${PKG_MGR/brew-missing/brew} install pipx; on Linux: python3 -m pip install --user pipx)"
  fi
}

check_ninja() {
  if command -v ninja >/dev/null 2>&1; then
    ok "Ninja $(extract_version "$(ninja --version 2>/dev/null)")"
  elif command -v samu >/dev/null 2>&1; then
    ok "samurai (ninja-compatible) found"
  else
    bad "Ninja not found"
    note_fail
    add_fix "Install Ninja (isolated, recommended):\n    pipx inject meson ninja   # if meson via pipx\n  or:\n    python3 -m pip install --user ninja\n  or via system packages:\n    $(pkg_cmd ninja)"
  fi
}

check_cmake() {
  if command -v cmake >/dev/null 2>&1; then
    local v; v="$(extract_version "$(cmake --version 2>/dev/null)")"
    if [ -n "$v" ] && version_ge "$v" "$CMAKE_MIN"; then
      ok "CMake $v (>= $CMAKE_MIN)"
    else
      bad "CMake ${v:-?} found but >= $CMAKE_MIN required (SUNDIALS/CppAD subprojects)"
      note_fail
      add_fix "Upgrade CMake:\n    $(pkg_cmd cmake)\n  or isolated:  python3 -m pip install --user --upgrade cmake"
    fi
  else
    bad "CMake not found (required to build the SUNDIALS and CppAD subprojects)"
    note_fail
    add_fix "Install CMake:\n    $(pkg_cmd cmake)\n  or isolated:  python3 -m pip install --user cmake"
  fi
}

check_pkgconfig() {
  if command -v pkg-config >/dev/null 2>&1; then
    ok "pkg-config $(extract_version "$(pkg-config --version 2>/dev/null)")"
  elif command -v pkgconf >/dev/null 2>&1; then
    ok "pkgconf $(extract_version "$(pkgconf --version 2>/dev/null)") (pkg-config compatible)"
  else
    bad "pkg-config not found"
    note_fail
    add_fix "Install pkg-config:\n    $(pkg_cmd pkgconfig)"
  fi
}

probe_cxx23() {
  local cxx="$1" dir src out rc
  command -v "$cxx" >/dev/null 2>&1 || return 2
  dir="$(mktemp -d 2>/dev/null)" || return 3
  src="$dir/gf_probe.cpp"
  out="$dir/gf_probe.out"
  cat > "$src" <<'EOF'
#include <print>
#include <format>
int main() { std::println("ok {}", 23); }
EOF
  "$cxx" -std=c++23 "$src" -o "$out" >/dev/null 2>&1
  rc=$?
  rm -rf "$dir"
  return $rc
}
check_toolchain() {
  local found_ok=0 detail=""

  local candidates=()
  [ -n "${CXX:-}" ] && candidates+=("$CXX")
  candidates+=(clang++ clang++-21 clang++-20 clang++-19 clang++-18 clang++-17 g++ g++-16 g++-15 g++-14)
  if [ "$OS" = "Darwin" ]; then
    candidates+=(/opt/homebrew/opt/llvm/bin/clang++ /usr/local/opt/llvm/bin/clang++)
  fi

  local seen="" cxx ver id minreq
  for cxx in "${candidates[@]}"; do
    command -v "$cxx" >/dev/null 2>&1 || continue
    local rp; rp="$(command -v "$cxx")"
    case " $seen " in *" $rp "*) continue ;; esac
    seen="$seen $rp"

    ver="$("$cxx" --version 2>/dev/null | head -n1)"
    if printf '%s' "$ver" | grep -qiE 'clang'; then
      id="clang"; minreq="$CLANG_MIN"
    elif printf '%s' "$ver" | grep -qiE 'g\+\+|gcc|free software foundation'; then
      id="gcc"; minreq="$GCC_MIN"
    else
      id="unknown"; minreq=""
    fi
    local vnum; vnum="$(extract_version "$ver")"

    if probe_cxx23 "$cxx"; then
      if [ -n "$minreq" ] && [ -n "$vnum" ] && ! version_ge "$vnum" "$minreq"; then
        warn "$cxx ($id $vnum): C++23 probe passed but below tested floor ($minreq) -- usable, upgrade if you hit errors"
      else
        ok "$cxx ($id ${vnum:-?}): C++23 standard library working"
      fi
      found_ok=1
      [ -z "${WORKING_CXX:-}" ] && WORKING_CXX="$rp"
    else
      detail="${detail}\n    - $cxx ($id ${vnum:-?}): present but failed the C++23 compile probe"
    fi
  done

  local have_cc=0
  for c in "${CC:-}" cc gcc clang; do
    [ -n "$c" ] && command -v "$c" >/dev/null 2>&1 && { have_cc=1; break; }
  done

  if [ "$found_ok" -eq 1 ] && [ "$have_cc" -eq 1 ]; then
    ok "C toolchain present"
    [ -n "${WORKING_CXX:-}" ] && info "Suggested: CXX=$WORKING_CXX meson setup build"
    return
  fi

  if [ "$found_ok" -ne 1 ]; then
    bad "No C++ compiler with working C++23 support found (need GCC >= $GCC_MIN or Clang >= $CLANG_MIN)"
    note_fail
    [ -n "$detail" ] && printf '       Compilers seen but rejected:%b\n' "$detail"
  fi
  if [ "$have_cc" -ne 1 ]; then
    bad "No C compiler found"
    note_fail
  fi


  local fix="Install a C/C++ toolchain (Clang >= $CLANG_MIN or GCC >= $GCC_MIN):\n    $(pkg_cmd clang)\n  or GCC:\n    $(pkg_cmd gcc)"
  # local fix="Install a C/C++ toolchain (GCC >= $GCC_MIN or Clang >= $CLANG_MIN):\n    $(pkg_cmd gcc)\n  or Clang:\n    $(pkg_cmd clang)"
  if [ "$OS" = "Linux" ]; then
    fix="$fix\n  Note: on Linux, Clang uses the system libstdc++ -- if a new Clang fails the\n  C++23 probe, install a newer libstdc++ (e.g. via g++-14) and Clang will use it."
  elif [ "$OS" = "Darwin" ]; then
    fix="$fix\n  Apple Clang: update Xcode / Command Line Tools to Xcode 16+\n    xcode-select --install   (or update via the App Store)"
  fi
  add_fix "$fix"
}

check_fortran() {
  [ "$WANT_FORTRAN" -eq 1 ] || return 0
  if command -v gfortran >/dev/null 2>&1; then
    local v; v="$(extract_version "$(gfortran --version 2>/dev/null | head -n1)")"
    if [ -n "$v" ] && version_ge "$v" "$GFORTRAN_MIN"; then
      ok "gfortran $v (>= $GFORTRAN_MIN)"
    else
      bad "gfortran ${v:-?} found but >= $GFORTRAN_MIN required"
      note_fail
      add_fix "Upgrade gfortran:\n    $(pkg_cmd gfortran)"
    fi
  else
    bad "gfortran not found (required by --fortran / -Dbuild_fortran=true)"
    note_fail
    add_fix "Install gfortran:\n    $(pkg_cmd gfortran)\n  (GridFire only supports gfortran; other Fortran compilers need -Dunsafe_fortran=true)"
  fi
}

printf '%s======================================================%s\n' "$BOLD" "$NC"
printf '%s  GridFire build-dependency check%s\n' "$BOLD" "$NC"
printf '%s======================================================%s\n' "$BOLD" "$NC"

detect_platform
info "Platform: $DISTRO_NAME ($ARCH)"
case "$PKG_MGR" in
  brew-missing)
    warn "Homebrew not found. Install it first: https://brew.sh"
    warn "All macOS remediation commands below assume Homebrew." ;;
  unknown)
    warn "Could not identify a supported package manager (apt/dnf/pacman/zypper)." ;;
  *) info "Package manager: $PKG_MGR" ;;
esac

hdr "Build tools"
check_meson
check_ninja
check_cmake
check_pkgconfig

hdr "Compiler toolchain"
check_toolchain

if [ "$WANT_FORTRAN" -eq 1 ]; then
  hdr "Fortran (requested via --fortran)"
  check_fortran
fi

printf '\n%s======================================================%s\n' "$BOLD" "$NC"
if [ "$FAIL_COUNT" -eq 0 ]; then
  printf '%s  All dependencies satisfied.%s\n' "$GREEN" "$NC"
  printf '%s======================================================%s\n' "$BOLD" "$NC"
  printf '\nYou can build GridFire with:\n\n'
  if [ -n "${WORKING_CXX:-}" ]; then
    printf '    CXX=%s meson setup build\n' "$WORKING_CXX"
  else
    printf '    meson setup build\n'
  fi
  printf '    meson compile -C build\n'
  printf '    sudo meson install -C build\n\n'
  exit 0
else
  printf '%s  %d dependency problem(s) found -- see below.%s\n' "$RED" "$FAIL_COUNT" "$NC"
  printf '%s======================================================%s\n' "$BOLD" "$NC"
  printf '\n%sThis script does not install anything. Run the commands below yourself:%s\n\n' "$BOLD" "$NC"
  printf '%b\n' "$FIX_NOTES"
  printf 'Re-run %s./install.sh%s after installing to confirm.\n\n' "$BOLD" "$NC"
  exit 1
fi