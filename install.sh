#!/usr/bin/env bash
#
# install.sh - Comprehensive interactive installation script for GridFire.
#
# This script performs the following actions:
# 1. Checks for essential system dependencies against required minimum versions.
# 2. If run with the --tui flag, it provides a comprehensive text-based user interface
#    to select and install dependencies, configure the build, and run build steps.
# 3. If run without flags, it prompts the user interactively for each missing dependency.
# 4. Provides detailed installation instructions for various Linux distributions and macOS.
# 5. Once all dependencies are met, it can run the meson setup, compile, install, and test commands.
# 6. Logs all operations to a file for easy debugging.

set -o pipefail

# --- Configuration ---
LOGFILE="GridFire_Installer.log"
NOTES_FILE="notes.txt"
MIN_GCC_VER="13.0.0"
MIN_CLANG_VER="16.0.0"
MIN_MESON_VER="1.5.0"

# --- Build Configuration Globals ---
BUILD_DIR="build"
INSTALL_PREFIX="/usr/local"
MESON_BUILD_TYPE="debug" # Default to debug to match Meson's default and avoid optimization warnings
MESON_LOG_LEVEL="info"
MESON_PKG_CONFIG="true"
MESON_NUM_CORES=$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1)
C_COMPILER=""
CC_COMPILER=""
FC_COMPILER=""
declare -A VALID_COMPILERS

# --- ANSI Color Codes ---
RED="\033[0;31m"
GREEN="\033[0;32m"
YELLOW="\033[0;33m"
BLUE="\033[0;34m"
NC="\033[0m" # No Color

# --- Logging and Utility Functions ---

# Log a message to stdout and to the log file.
log() {
  local message="$1"
  # Print the colored message to the terminal.
  echo -e "${message}"
  # Strip ANSI color codes and append to the log file.
  echo -e "${message}" | sed 's/\x1B\[[0-9;]*[mK]//g' >> "$LOGFILE"
}

# Version comparison function. Returns 0 if v1=v2, 1 if v1>v2, 2 if v1<v2
vercomp() {
    local v1=$1 v2=$2
    if [[ "$v1" == "$v2" ]]; then return 0; fi
    local IFS=.
    read -ra v1_parts <<< "$v1"
    read -ra v2_parts <<< "$v2"

    local i
    # Pad the shorter version array with zeros
    local len1=${#v1_parts[@]}
    local len2=${#v2_parts[@]}
    if (( len1 > len2 )); then
        for ((i=len2; i<len1; i++)); do v2_parts[i]=0; done
    elif (( len2 > len1 )); then
        for ((i=len1; i<len2; i++)); do v1_parts[i]=0; done
    fi

    for i in "${!v1_parts[@]}"; do
        local p1=${v1_parts[i]}
        local p2=${v2_parts[i]}

        # Ensure parts are integers before comparison
        if ! [[ "$p1" =~ ^[0-9]+$ ]]; then p1=0; fi
        if ! [[ "$p2" =~ ^[0-9]+$ ]]; then p2=0; fi

        if (( p1 > p2 )); then return 1; fi
        if (( p1 < p2 )); then return 2; fi
    done

    return 0
}


# Show the help message and exit.
show_help() {
  echo "Usage: ./install.sh [options]"
  echo
  echo "This script checks for dependencies, installs them, and builds the project."
  echo
  echo "Options:"
  echo "  --tui         Run in Text-based User Interface mode for interactive dependency installation and build control."
  echo "  --help, -h    Show this help message and exit."
  echo "  --clean       Remove the build directory and log file before starting."
  echo
  echo "The script will automatically detect your OS and suggest the correct package manager commands."
  echo "All output is logged to ${LOGFILE}."
}

# --- OS Detection ---
OS_NAME=$(uname -s)
DISTRO_ID="unknown"
if [ -f /etc/os-release ]; then
  . /etc/os-release
  DISTRO_ID=$ID
fi
if [[ "$OS_NAME" == "Darwin" ]]; then
  OS_NAME="macOS"
fi

# --- Dependency Check Functions ---

check_command() {
  command -v "$1" &>/dev/null
}

set_compilers() {
    if [[ "$CC_COMPILER" == *"clang++"* ]]; then
        C_COMPILER=$(echo "$CC_COMPILER" | sed 's/clang++/clang/')
        # Try to find a corresponding flang or fallback to gfortran
        local fc_ver=$(echo "$C_COMPILER" | grep -oE '[0-9]+')
        if check_command "flang-$fc_ver"; then
            FC_COMPILER="flang-$fc_ver"
        elif check_command "flang"; then
            FC_COMPILER="flang"
        elif check_command "gfortran"; then
            FC_COMPILER="gfortran"
        else
            FC_COMPILER=""
        fi
    elif [[ "$CC_COMPILER" == *"g++"* ]]; then
        C_COMPILER=$(echo "$CC_COMPILER" | sed 's/g++/gcc/')
        FC_COMPILER=$(echo "$CC_COMPILER" | sed 's/g++/gfortran/')
        if ! check_command "$FC_COMPILER"; then
             # Fallback to generic gfortran if versioned one not found
             if check_command "gfortran"; then
                 FC_COMPILER="gfortran"
             else
                 FC_COMPILER=""
             fi
        fi
    else
        C_COMPILER=""
        FC_COMPILER=""
    fi
}

check_compiler() {
    VALID_COMPILERS=()
    local potential_gccs=("g++-14" "g++-13" "g++")
    local potential_clangs=("clang++-17" "clang++-16" "clang++")

    for cc in "${potential_gccs[@]}"; do
        if check_command "$cc"; then
            local ver; ver=$($cc -dumpversion | grep -oE '[0-9]+(\.[0-9]+)*' | head -n1)
            if [[ -n "$ver" ]]; then
                vercomp "$ver" "$MIN_GCC_VER"
                if [[ $? -ne 2 ]]; then
                    VALID_COMPILERS["g++ ($ver)"]="$cc"
                    break # Found a valid one, stop searching for g++
                fi
            fi
        fi
    done

    for cc in "${potential_clangs[@]}"; do
        if check_command "$cc"; then
            local ver; ver=$($cc --version | head -n1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -n1)
            if [[ -n "$ver" ]]; then
                vercomp "$ver" "$MIN_CLANG_VER"
                if [[ $? -ne 2 ]]; then
                    VALID_COMPILERS["clang++ ($ver)"]="$cc"
                    break # Found a valid one
                fi
            fi
        fi
    done

    if [ ${#VALID_COMPILERS[@]} -gt 0 ]; then
        if ! [[ " ${VALID_COMPILERS[*]} " =~ " ${CC_COMPILER} " ]]; then
            # Correctly get the first value from the associative array
            local keys=("${!VALID_COMPILERS[@]}")
            local first_key="${keys[0]}"
            CC_COMPILER="${VALID_COMPILERS[$first_key]}"
        fi
        set_compilers
        log "${GREEN}[OK] Found valid C++ compiler(s). Using C: '${C_COMPILER}', C++: '${CC_COMPILER}', FC: '${FC_COMPILER}'.${NC}"
        return 0
    else
        log "${RED}[FAIL] No valid C++ compiler found (GCC >= ${MIN_GCC_VER} or Clang >= ${MIN_CLANG_VER}).${NC}"
        CC_COMPILER=""
        C_COMPILER=""
        FC_COMPILER=""
        return 1
    fi
}

check_python_dev() {
  if check_command python3 && python3-config --includes &>/dev/null; then
    log "${GREEN}[OK] Found Python 3 development headers.${NC}"
    return 0
  else
    log "${RED}[FAIL] Python 3 development headers not found.${NC}"
    return 1
  fi
}

check_meson_python() {
    if python3 -c "import mesonpy" &>/dev/null; then
        log "${GREEN}[OK] Found meson-python package.${NC}"
        return 0
    else
        log "${RED}[FAIL] meson-python package not found.${NC}"
        return 1
    fi
}

check_cmake() {
  if check_command cmake; then
    log "${GREEN}[OK] Found CMake: $(cmake --version | head -n1)${NC}"
    return 0
  else
    log "${RED}[FAIL] CMake not found.${NC}"
    return 1
  fi
}

check_meson() {
    if check_command meson; then
        local ver; ver=$(meson --version)
        vercomp "$ver" "$MIN_MESON_VER"
        if [[ $? -ne 2 ]]; then
            log "${GREEN}[OK] Found Meson ${ver}.${NC}"
            return 0
        else
            log "${RED}[FAIL] Meson version ${ver} is too old. Need >= ${MIN_MESON_VER}.${NC}"
            return 1
        fi
    else
        log "${RED}[FAIL] Meson not found.${NC}"
        return 1
    fi
}

check_boost() {
  log "${BLUE}[Info] Checking for Boost dependency using Meson...${NC}"
  local test_dir="meson-boost-test"
  rm -rf "$test_dir" && mkdir -p "$test_dir"
  cat > "$test_dir/meson.build" <<EOF
project('boost-check', 'cpp')
boost_dep = dependency('boost', required: true)
EOF
  if meson setup "$test_dir/build" "$test_dir" &>/dev/null; then
    log "${GREEN}[OK] Found Boost libraries.${NC}"
    rm -rf "$test_dir"
    return 0
  else
    log "${RED}[FAIL] Boost libraries not found by Meson.${NC}"
    rm -rf "$test_dir"
    return 1
  fi
}

# --- Dependency Installation Functions ---

get_compiler_install_cmd() {
    local compiler_to_install="$1" # e.g., g++ or clang++
    local cmd=""
    case "$OS_NAME" in
        "macOS")
            local brew_cmd; brew_cmd=$(command -v brew)
            case "$compiler_to_install" in
                "g++") cmd="$brew_cmd install gcc" ;;
                "clang++") cmd="$brew_cmd install llvm" ;; # llvm provides clang
            esac
            ;;
        "Linux")
            case "$DISTRO_ID" in
                "ubuntu"|"debian"|"linuxmint")
                    case "$compiler_to_install" in
                        "g++") cmd="sudo apt-get install -y g++-13 gfortran-13" ;;
                        "clang++") cmd="sudo apt-get install -y clang-16" ;;
                    esac
                    ;;
                "fedora")
                    # Fedora usually has recent versions in main repos
                    case "$compiler_to_install" in
                        "g++") cmd="sudo dnf install -y gcc-c++ gcc-gfortran" ;;
                        "clang++") cmd="sudo dnf install -y clang" ;;
                    esac
                    ;;
                "arch"|"manjaro")
                     # Arch is rolling release, should be fine
                    case "$compiler_to_install" in
                        "g++") cmd="sudo pacman -S --noconfirm gcc gcc-fortran" ;;
                        "clang++") cmd="sudo pacman -S --noconfirm clang" ;;
                    esac
                    ;;
            esac
            ;;
    esac
    echo "$cmd"
}


get_install_cmd() {
  local dep_name="$1"
  local cmd=""
  case "$OS_NAME" in
    "macOS")
      [ ! -f /opt/homebrew/bin/brew ] && [ ! -f /usr/local/bin/brew ] && log "${YELLOW}Homebrew not found. Please install from https://brew.sh${NC}"
      local brew_cmd
      brew_cmd=$(command -v brew)
      case "$dep_name" in
        "compiler") cmd="$brew_cmd install gcc llvm" ;; # Install both
        "python-dev") cmd="$brew_cmd install python3" ;;
        "meson-python") cmd="python3 -m pip install meson-python" ;;
        "meson") cmd="python3 -m pip install --upgrade meson" ;;
        "cmake") cmd="$brew_cmd install cmake" ;;
        "boost") cmd="$brew_cmd install boost" ;;
        "dialog") cmd="$brew_cmd install dialog" ;;
      esac
      ;;
    "Linux")
      case "$DISTRO_ID" in
        "ubuntu"|"debian"|"linuxmint")
          case "$dep_name" in
            "compiler") cmd="sudo apt-get install -y g++-13 gfortran-13 clang-16" ;;
            "python-dev") cmd="sudo apt-get install -y python3-dev" ;;
            "meson-python") cmd="python3 -m pip install meson-python" ;;
            "meson") cmd="python3 -m pip install --upgrade meson" ;;
            "cmake") cmd="sudo apt-get install -y cmake" ;;
            "boost") cmd="sudo apt-get install -y libboost-all-dev" ;;
            "dialog") cmd="sudo apt-get install -y dialog" ;;
          esac
          ;;
        "fedora")
          case "$dep_name" in
            "compiler") cmd="sudo dnf install -y gcc-c++ gcc-gfortran clang" ;;
            "python-dev") cmd="sudo dnf install -y python3-devel" ;;
            "meson-python") cmd="python3 -m pip install meson-python" ;;
            "meson") cmd="python3 -m pip install --upgrade meson" ;;
            "cmake") cmd="sudo dnf install -y cmake" ;;
            "boost") cmd="sudo dnf install -y boost-devel" ;;
            "dialog") cmd="sudo dnf install -y dialog" ;;
          esac
          ;;
        "arch"|"manjaro")
          case "$dep_name" in
            "compiler") cmd="sudo pacman -S --noconfirm gcc gcc-fortran clang" ;;
            "python-dev") cmd="sudo pacman -S --noconfirm python" ;;
            "meson-python") cmd="python3 -m pip install meson-python" ;;
            "meson") cmd="python3 -m pip install --upgrade meson" ;;
            "cmake") cmd="sudo pacman -S --noconfirm cmake" ;;
            "boost") cmd="sudo pacman -S --noconfirm boost" ;;
            "dialog") cmd="sudo pacman -S --noconfirm dialog" ;;
          esac
          ;;
        *) log "${YELLOW}[Warn] Unsupported Linux distribution: ${DISTRO_ID}.${NC}" ;;
      esac
      ;;
    *) log "${YELLOW}[Warn] Unsupported OS: ${OS_NAME}.${NC}" ;;
  esac
  echo "$cmd"
}

# --- Build Functions ---

run_meson_setup() {
    log "\n${BLUE}--- Configuring Meson Build ---${NC}"
    if [ ! -f "meson.build" ]; then
        log "${RED}[FATAL] meson.build file not found. Cannot proceed.${NC}"; return 1;
    fi
    if [ -z "$CC_COMPILER" ] || [ -z "$C_COMPILER" ]; then
        log "${RED}[FATAL] No valid C/C++ compiler selected. Configure one first.${NC}"; return 1;
    fi
    local reconfigure_flag=""
    if [ -d "$BUILD_DIR" ]; then
        reconfigure_flag="--reconfigure"
        log "${YELLOW}[Info] Existing build directory found. Will reconfigure.${NC}"
    fi

    local meson_opts=()
    meson_opts+=("-Dbuildtype=${MESON_BUILD_TYPE}")
    meson_opts+=("-Dlog_level=${MESON_LOG_LEVEL}")
    meson_opts+=("-Dpkg-config=${MESON_PKG_CONFIG}")
    meson_opts+=("--prefix=${INSTALL_PREFIX}")

    log "${BLUE}[Info] Using C compiler:       ${C_COMPILER}${NC}"
    log "${BLUE}[Info] Using C++ compiler:     ${CC_COMPILER}${NC}"
    log "${BLUE}[Info] Using Fortran compiler: ${FC_COMPILER}${NC}"
    log "${BLUE}[Info] Running meson setup with options: ${meson_opts[*]}${NC}"
    # Set CC, CXX, and FC environment variables for the meson command
    if ! CC="${C_COMPILER}" CXX="${CC_COMPILER}" FC="${FC_COMPILER}" meson setup "${BUILD_DIR}" "${meson_opts[@]}" ${reconfigure_flag}; then
        log "${RED}[FATAL] Meson setup failed. See log for details.${NC}"; return 1;
    fi
    log "${GREEN}[Success] Meson setup complete.${NC}"
}

run_meson_compile() {
    log "\n${BLUE}--- Compiling Project ---${NC}"
    if [ ! -d "$BUILD_DIR" ]; then
        log "${RED}[FATAL] Build directory not found. Run setup first.${NC}"; return 1;
    fi
    log "${BLUE}[Info] Running meson compile with ${MESON_NUM_CORES} cores...${NC}"
    if ! meson compile -C "${BUILD_DIR}" -j "${MESON_NUM_CORES}"; then
        log "${RED}[FATAL] Meson compile failed. See log for details.${NC}"; return 1;
    fi
    log "${GREEN}[Success] Meson compile complete.${NC}"
}

run_meson_install() {
    log "\n${BLUE}--- Installing Project ---${NC}"
    if [ ! -d "$BUILD_DIR" ]; then
        log "${RED}[FATAL] Build directory not found. Run setup and compile first.${NC}"; return 1;
    fi
    log "${BLUE}[Info] Running meson install (prefix: ${INSTALL_PREFIX})...${NC}"
    if ! sudo meson install -C "${BUILD_DIR}"; then
        log "${RED}[FATAL] Meson install failed. See log for details.${NC}"; return 1;
    fi
    log "${GREEN}[Success] Meson install complete.${NC}"
}

run_meson_tests() {
    log "\n${BLUE}--- Running Tests ---${NC}"
    if [ ! -d "$BUILD_DIR" ]; then
        log "${RED}[FATAL] Build directory not found. Run setup and compile first.${NC}"; return 1;
    fi
    log "${BLUE}[Info] Running meson test...${NC}"
    if ! meson test -C "${BUILD_DIR}"; then
        log "${RED}[FATAL] Meson tests failed. See log for details.${NC}"; return 1;
    fi
    log "${GREEN}[Success] Tests passed.${NC}"
}

# --- TUI Functions ---

check_dialog_installed() {
  if ! check_command dialog; then
    log "${YELLOW}[Warn] The 'dialog' utility is required for TUI mode.${NC}"
    local install_cmd
    install_cmd=$(get_install_cmd "dialog")
    if [ -n "$install_cmd" ]; then
      if prompt_yes_no "Attempt to install it now? (y/n):"; then
        eval "$install_cmd"
        if ! check_command dialog; then
          log "${RED}[FATAL] Failed to install 'dialog'. Cannot run in TUI mode.${NC}"
          return 1
        fi
      else
        log "${RED}[FATAL] Cannot run in TUI mode without 'dialog'. Exiting.${NC}"
        return 1
      fi
    else
      log "${RED}[FATAL] Cannot auto-install 'dialog'. Please install it manually.${NC}"
      return 1
    fi
  fi
  return 0
}

run_dependency_installer_tui() {
  # This function now just calls the check functions to populate status
  declare -A DEP_STATUS
  check_compiler >/dev/null; DEP_STATUS[compiler]=$?
  check_python_dev >/dev/null; DEP_STATUS[python-dev]=$?
  check_meson_python >/dev/null; DEP_STATUS[meson-python]=$?
  check_cmake >/dev/null; DEP_STATUS[cmake]=$?
  check_meson >/dev/null; DEP_STATUS[meson]=$?
  check_boost >/dev/null; DEP_STATUS[boost]=$?

  local choices
  choices=$(dialog --clear --backtitle "Project Dependency Installer" \
    --title "Install System Dependencies" \
    --checklist "Select dependencies to install. Already found dependencies are unchecked." 20 70 6 \
    "compiler" "C++ Compilers (g++, clang++)" "$([[ ${DEP_STATUS[compiler]} -ne 0 ]] && echo "on" || echo "off")" \
    "python-dev" "Python 3 Dev Headers" "$([[ ${DEP_STATUS[python-dev]} -ne 0 ]] && echo "on" || echo "off")" \
    "meson-python" "meson-python (for Python bindings)" "$([[ ${DEP_STATUS[meson-python]} -ne 0 ]] && echo "on" || echo "off")" \
    "cmake" "CMake" "$([[ ${DEP_STATUS[cmake]} -ne 0 ]] && echo "on" || echo "off")" \
    "meson" "Meson Build System (>=${MIN_MESON_VER})" "$([[ ${DEP_STATUS[meson]} -ne 0 ]] && echo "on" || echo "off")" \
    "boost" "Boost Libraries" "$([[ ${DEP_STATUS[boost]} -ne 0 ]] && echo "on" || echo "off")" \
    3>&1 1>&2 2>&3)

  clear
  if [ -z "$choices" ]; then log "${YELLOW}[Info] No dependencies selected.${NC}"; return; fi

  for choice in $choices; do
    local dep; dep=$(echo "$choice" | tr -d '"')
    log "\n${BLUE}--- Installing ${dep} ---${NC}"
    local install_cmd; install_cmd=$(get_install_cmd "$dep")
    if [ -n "$install_cmd" ]; then
      eval "$install_cmd" 2>&1 | tee -a "$LOGFILE"
    else
      log "${RED}[Error] No automatic installation command for '${dep}'. Please install manually.${NC}"
      dialog --msgbox "Could not find an automatic installation command for '${dep}' on your system. Please install it manually." 8 60
    fi
  done
  # Re-run check to update status
  check_compiler
}

run_python_bindings_tui() {
    if ! check_meson_python; then
        dialog --msgbox "The 'meson-python' package is required to build Python bindings. Please install it from the 'Install System Dependencies' menu first." 10 70
        return
    fi
    if [ -z "$CC_COMPILER" ]; then
        dialog --msgbox "No C++ compiler is selected. Please configure one from the 'Configure Build Options' menu before building the Python bindings." 10 70
        return
    fi

    local python_exec
    python_exec=$(command -v python3)
    if [ -z "$python_exec" ]; then
        dialog --msgbox "Could not find 'python3' executable. Please ensure Python 3 is installed and in your PATH." 8 60
        return
    fi

    local choice
    choice=$(dialog --clear --backtitle "Python Bindings Installer" \
        --title "Install Python Bindings" \
        --menu "Using Python: ${python_exec}\nUsing C++ Compiler: ${CC_COMPILER}\n\nSelect installation mode:" 15 70 2 \
        "1" "Developer Mode (pip install -e .)" \
        "2" "User Mode (pip install .)" \
        3>&1 1>&2 2>&3)

    clear
    case "$choice" in
        1)
            log "${BLUE}[Info] Installing Python bindings in Developer Mode...${NC}"
            if ! CC="${C_COMPILER}" CXX="${CC_COMPILER}" FC="${FC_COMPILER}" pip install -e . --no-build-isolation -vv; then
                log "${RED}[Error] Failed to install Python bindings in developer mode.${NC}"
                dialog --msgbox "Developer mode installation failed. Check the log for details." 8 60
            else
                log "${GREEN}[Success] Python bindings installed in developer mode.${NC}"
                dialog --msgbox "Successfully installed Python bindings in developer mode." 8 60
            fi
            ;;
        2)
            log "${BLUE}[Info] Installing Python bindings in User Mode...${NC}"
            if ! CC="${C_COMPILER}" CXX="${CC_COMPILER}" FC="${FC_COMPILER}" pip install .; then
                log "${RED}[Error] Failed to install Python bindings in user mode.${NC}"
                dialog --msgbox "User mode installation failed. Check the log for details." 8 60
            else
                log "${GREEN}[Success] Python bindings installed in user mode.${NC}"
                dialog --msgbox "Successfully installed Python bindings in user mode." 8 60
            fi
            ;;
        *)
            log "${YELLOW}[Info] Python binding installation cancelled.${NC}"
            ;;
    esac
}

run_compiler_selection_tui() {
    local gpp_ver; gpp_ver=$(g++ -dumpversion 2>/dev/null | grep -oE '[0-9]+(\.[0-9]+)*' | head -n1)
    local clang_ver; clang_ver=$(clang++ --version 2>/dev/null | head -n1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -n1)

    local gpp_ok=false
    local clang_ok=false
    if [[ -n "$gpp_ver" ]]; then vercomp "$gpp_ver" "$MIN_GCC_VER"; [[ $? -ne 2 ]] && gpp_ok=true; fi
    if [[ -n "$clang_ver" ]]; then vercomp "$clang_ver" "$MIN_CLANG_VER"; [[ $? -ne 2 ]] && clang_ok=true; fi

    if ! $gpp_ok && ! $clang_ok; then
        # No valid compilers found
        local choices
        choices=$(dialog --title "Compiler Installation" --checklist "No valid C++ compiler found. Please select which to install:" 15 70 2 \
            "g++" "GNU C++ Compiler (>=${MIN_GCC_VER})" "on" \
            "clang++" "Clang C++ Compiler (>=${MIN_CLANG_VER}, often faster)" "off" 3>&1 1>&2 2>&3)
        if [ -n "$choices" ]; then
            for choice in $choices; do
                local compiler_to_install; compiler_to_install=$(echo "$choice" | tr -d '"')
                local install_cmd; install_cmd=$(get_compiler_install_cmd "$compiler_to_install")
                if [ -n "$install_cmd" ]; then eval "$install_cmd" 2>&1 | tee -a "$LOGFILE"; fi
            done
        fi
    elif ! $gpp_ok && [[ -n "$gpp_ver" ]]; then
        # g++ found but too old
        if dialog --title "Compiler Update" --yesno "Found g++ version ${gpp_ver}, but require >= ${MIN_GCC_VER}.\n\nAttempt to install a compatible version?" 10 70; then
            local install_cmd; install_cmd=$(get_compiler_install_cmd "g++")
            if [ -n "$install_cmd" ]; then eval "$install_cmd" 2>&1 | tee -a "$LOGFILE"; fi
        fi
    elif ! $clang_ok && [[ -n "$clang_ver" ]]; then
        # clang++ found but too old
        if dialog --title "Compiler Update" --yesno "Found clang++ version ${clang_ver}, but require >= ${MIN_CLANG_VER}.\n\nAttempt to install a compatible version?" 10 70; then
            local install_cmd; install_cmd=$(get_compiler_install_cmd "clang++")
            if [ -n "$install_cmd" ]; then eval "$install_cmd" 2>&1 | tee -a "$LOGFILE"; fi
        fi
    fi

    # Re-check compilers and let user choose if multiple are available
    check_compiler
    if [ ${#VALID_COMPILERS[@]} -gt 0 ]; then
        local menu_items=()
        for name in "${!VALID_COMPILERS[@]}"; do
             menu_items+=("$name" "")
        done
        local compiler_choice_key
        compiler_choice_key=$(dialog --title "Select C++ Compiler" --menu "Select the C++ compiler to use:" 15 70 ${#VALID_COMPILERS[@]} "${menu_items[@]}" 3>&1 1>&2 2>&3)
        if [ -n "$compiler_choice_key" ]; then
            CC_COMPILER="${VALID_COMPILERS[$compiler_choice_key]}"
            set_compilers
            log "${BLUE}[Config] Set C compiler to: ${C_COMPILER}${NC}"
            log "${BLUE}[Config] Set C++ compiler to: ${CC_COMPILER}${NC}"
            log "${BLUE}[Config] Set Fortran compiler to: ${FC_COMPILER}${NC}"
        fi
    else
        dialog --msgbox "No valid C++ compiler could be found or installed. Please install one manually that meets the version requirements." 8 70
    fi
}

run_build_config_tui() {
    local choice
    choice=$(dialog --clear --backtitle "Build Configuration" \
        --title "Configure Build Options" \
        --menu "Select an option to configure:" 20 70 7 \
        "1" "Build Directory (current: ${BUILD_DIR})" \
        "2" "Install Prefix (current: ${INSTALL_PREFIX})" \
        "3" "Manage & Select C/C++/Fortran Compiler" \
        "4" "Build Type (current: ${MESON_BUILD_TYPE})" \
        "5" "Log Level (current: ${MESON_LOG_LEVEL})" \
        "6" "Generate pkg-config (current: ${MESON_PKG_CONFIG})" \
        "7" "Number of Cores (current: ${MESON_NUM_CORES})" \
        3>&1 1>&2 2>&3)

    clear
    case "$choice" in
        1)
            local new_dir
            new_dir=$(dialog --title "Set Build Directory" --inputbox "Enter new build directory name:" 10 60 "${BUILD_DIR}" 3>&1 1>&2 2>&3)
            if [ -n "$new_dir" ]; then
                BUILD_DIR="$new_dir"
                log "${BLUE}[Config] Set build directory to: ${BUILD_DIR}${NC}"
            fi
            ;;
        2)
            local new_prefix
            new_prefix=$(dialog --title "Set Install Prefix" --inputbox "Enter absolute path for installation prefix:" 10 60 "${INSTALL_PREFIX}" 3>&1 1>&2 2>&3)
            if [ -n "$new_prefix" ]; then
                INSTALL_PREFIX="$new_prefix"
                log "${BLUE}[Config] Set install prefix to: ${INSTALL_PREFIX}${NC}"
            fi
            ;;
        3)
            run_compiler_selection_tui
            ;;
        4)
            local build_type_choice
            build_type_choice=$(dialog --title "Select Build Type" --menu "" 15 70 4 \
                "debug" "No optimizations, with debug symbols" \
                "release" "Optimized for performance" \
                "debugoptimized" "With debug symbols and optimization" \
                "plain" "Custom flags only" \
                3>&1 1>&2 2>&3)
            if [ -n "$build_type_choice" ]; then
                MESON_BUILD_TYPE="$build_type_choice"
                log "${BLUE}[Config] Set build type to: ${MESON_BUILD_TYPE}${NC}"
            fi
            ;;
        5)
            local log_level_choice
            log_level_choice=$(dialog --title "Select Log Level" --menu "" 15 70 8 \
                "traceL3" "" "traceL2" "" "traceL1" "" "debug" "" "info" "" "warning" "" "error" "" "critical" "" \
                3>&1 1>&2 2>&3)
            if [ -n "$log_level_choice" ]; then
                MESON_LOG_LEVEL="$log_level_choice"
                log "${BLUE}[Config] Set log level to: ${MESON_LOG_LEVEL}${NC}"
            fi
            ;;
        6)
            if dialog --title "Generate pkg-config" --yesno "Generate gridfire.pc file?" 7 60; then
                MESON_PKG_CONFIG="true"
            else
                MESON_PKG_CONFIG="false"
            fi
            log "${BLUE}[Config] Set pkg-config generation to: ${MESON_PKG_CONFIG}${NC}"
            ;;
        7)
            local max_cores; max_cores=$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1)
            local core_choice
            core_choice=$(dialog --title "Set Number of Cores" --inputbox "Enter number of cores for compilation.\nAvailable: ${max_cores}" 10 60 "${MESON_NUM_CORES}" 3>&1 1>&2 2>&3)
            if [[ "$core_choice" =~ ^[0-9]+$ ]] && [ "$core_choice" -gt 0 ]; then
                MESON_NUM_CORES="$core_choice"
                log "${BLUE}[Config] Set number of cores to: ${MESON_NUM_CORES}${NC}"
            elif [ -n "$core_choice" ]; then
                 dialog --msgbox "Invalid input. Please enter a positive number." 6 40
            fi
            ;;
    esac
}

run_notes_tui() {
    if [ ! -f "$NOTES_FILE" ]; then
        dialog --msgbox "Notes file '${NOTES_FILE}' not found." 8 50
        return
    fi

    local notes_content=""
    local counter=1
    # Read file, filter comments, and process non-empty lines
    while IFS= read -r line; do
        # Skip empty or comment lines
        if [[ -z "$line" || "$line" =~ ^[[:space:]]*# ]]; then
            continue
        fi
        echo "$line"
        notes_content+="${counter}. ${line}\n\n"
        ((counter++))
    done < "$NOTES_FILE"

    if [ -z "$notes_content" ]; then
        dialog --msgbox "No notes found in '${NOTES_FILE}'." 8 50
    else
        dialog --title "Installer Notes" --msgbox "$notes_content" 20 70
    fi
}

run_main_tui() {
    if ! check_dialog_installed; then return 1; fi
    # Initial check to populate compiler list and set a default
    check_compiler

    local sudo_status="User Mode"
    if [ "$EUID" -eq 0 ]; then
        sudo_status="Root/Sudo Mode"
    fi

    while true; do
        local choice
        choice=$(dialog --clear --backtitle "GridFire Installer - [${sudo_status}]" \
            --title "Main Menu" \
            --menu "C: ${C_COMPILER:-N/A} C++: ${CC_COMPILER:-N/A} FC: ${FC_COMPILER:-N/A}\nDIR: ${BUILD_DIR} | TYPE: ${MESON_BUILD_TYPE} | CORES: ${MESON_NUM_CORES}\nPREFIX: ${INSTALL_PREFIX}\nLOG: ${MESON_LOG_LEVEL} | PKG-CONFIG: ${MESON_PKG_CONFIG}" 22 78 11 \
            "1" "Install System Dependencies" \
            "2" "Configure Build Options" \
            "3" "Install Python Bindings" \
            "4" "Run Full Build (Setup + Compile)" \
            "5" "Run Meson Setup/Reconfigure" \
            "6" "Run Meson Compile" \
            "7" "Run Meson Install (requires sudo)" \
            "8" "Run Tests" \
            "9" "View Notes" \
            "10" "Exit" \
            3>&1 1>&2 2>&3)

        clear
        case "$choice" in
            1) run_dependency_installer_tui ;;
            2) run_build_config_tui ;;
            3) run_python_bindings_tui ;;
            4) run_meson_setup && run_meson_compile ;;
            5) run_meson_setup ;;
            6) run_meson_compile ;;
            7) run_meson_install ;;
            8) run_meson_tests ;;
            9) run_notes_tui ;;
            10) break ;;
            *) log "${YELLOW}[Info] TUI cancelled.${NC}"; break ;;
        esac
    done
    clear
}

# --- Script Entry Point ---
main() {
  if [[ " $@ " =~ " --help " ]] || [[ " $@ " =~ " -h " ]]; then show_help; exit 0; fi
  if [[ " $@ " =~ " --clean " ]]; then log "${BLUE}[Info] Cleaning up...${NC}"; rm -rf "$BUILD_DIR" "$LOGFILE"; fi

  echo "" > "$LOGFILE" # Clear log file
  log "--- GridFire Installation Log ---"
  log "Date: $(date)"
  log "OS: ${OS_NAME}, Distro: ${DISTRO_ID}"

  if [[ " $@ " =~ " --tui " ]]; then
    run_main_tui
    log "${GREEN}Exited TUI mode.${NC}"
    exit 0
  fi

  # --- Non-TUI path ---
  log "\n${BLUE}--- Checking System Dependencies (CLI Mode) ---${NC}"
  # Run check_compiler first to set the default
  check_compiler
  declare -A CHECKS=(
      [python-dev]="check_python_dev" [meson-python]="check_meson_python" [cmake]="check_cmake"
      [meson]="check_meson" [boost]="check_boost"
  )
  if ! check_compiler; then
      # Handle case where no compiler is found
      local install_cmd; install_cmd=$(get_install_cmd "compiler")
      if [ -n "$install_cmd" ]; then
          if prompt_yes_no "Dependency 'compiler' is missing. Attempt to install? (y/n):"; then
              eval "$install_cmd" 2>&1 | tee -a "$LOGFILE"
          fi
      else
          log "${RED}[Error] No automatic installation for 'compiler'. Please install manually.${NC}"
      fi
  fi

  for dep in "${!CHECKS[@]}"; do
      if ! ${CHECKS[$dep]}; then
          local install_cmd; install_cmd=$(get_install_cmd "$dep")
          if [ -n "$install_cmd" ]; then
              if prompt_yes_no "Dependency '${dep}' is missing. Attempt to install? (y/n):"; then
                  eval "$install_cmd" 2>&1 | tee -a "$LOGFILE"
              fi
          else
              log "${RED}[Error] No automatic installation for '${dep}'. Please install manually.${NC}"
          fi
      fi
  done

  log "\n${BLUE}--- Re-checking all dependencies ---${NC}"
  local final_fail=false
  # Re-add compiler check to the list for final verification
  CHECKS[compiler]="check_compiler"
  for dep in "${!CHECKS[@]}"; do
      if ! ${CHECKS[$dep]}; then
          log "${RED}[FATAL] Dependency still missing: ${dep}${NC}"
          final_fail=true
      fi
  done

  if $final_fail; then
      log "${RED}Please install missing dependencies and re-run.${NC}"
      exit 1
  fi

  log "\n${GREEN}--- All dependencies met. Proceeding with build. ---${NC}"
  run_meson_setup && run_meson_compile
}

main "$@"
