#!/usr/bin/env bash
#
# install.sh - Comprehensive interactive installation script for GridFire.
#
# This script performs the following actions:
# 1. Checks for essential system dependencies.
# 2. If run with the --tui flag, it provides a comprehensive text-based user interface
#    to select and install dependencies, configure the build, and run build steps.
# 3. If run without flags, it prompts the user interactively for each missing dependency.
# 4. Provides detailed installation instructions for various Linux distributions and macOS.
# 5. Once all dependencies are met, it can run the meson setup, compile, install, and test commands.
# 6. Logs all operations to a file for easy debugging.

set -o pipefail

# --- Configuration ---
LOGFILE="GridFire_Installer.log"

# --- Build Configuration Globals ---
BUILD_DIR="build"
INSTALL_PREFIX="/usr/local"
MESON_BUILD_TYPE="release"
MESON_LOG_LEVEL="info"
MESON_PKG_CONFIG="true"
MESON_NUM_CORES=$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1)

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

# Prompt the user for a yes/no answer.
prompt_yes_no() {
  local prompt_msg="$1"
  local answer
  while true; do
    read -p "$(echo -e "${YELLOW}${prompt_msg}${NC}") " answer
    case "$answer" in
      [Yy]* ) return 0;;
      [Nn]* ) return 1;;
      * ) echo "Please answer yes or no.";;
    esac
  done
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

check_compiler() {
  if check_command g++; then
    log "${GREEN}[OK] Found C++ compiler: g++${NC}"
    return 0
  elif check_command clang++; then
    log "${GREEN}[OK] Found C++ compiler: clang++${NC}"
    return 0
  else
    log "${RED}[FAIL] No C++ compiler (g++ or clang++) found.${NC}"
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
    log "${GREEN}[OK] Found Meson: $(meson --version)${NC}"
    return 0
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

run_install_cmd() {
  local cmd="$1"
  log "${BLUE}[Exec] Running: ${cmd}${NC}"
  if prompt_yes_no "Execute this command? (y/n):"; then
    eval "${cmd}" 2>&1 | tee -a "$LOGFILE"
    return $?
  else
    log "${YELLOW}[Skip] User skipped installation command.${NC}"
    return 1
  fi
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
        "compiler") cmd="xcode-select --install" ;;
        "python-dev") cmd="$brew_cmd install python3" ;;
        "cmake") cmd="$brew_cmd install cmake" ;;
        "meson") cmd="$brew_cmd install meson" ;;
        "boost") cmd="$brew_cmd install boost" ;;
        "dialog") cmd="$brew_cmd install dialog" ;;
      esac
      ;;
    "Linux")
      case "$DISTRO_ID" in
        "ubuntu"|"debian"|"linuxmint")
          case "$dep_name" in
            "compiler") cmd="sudo apt-get install -y build-essential" ;;
            "python-dev") cmd="sudo apt-get install -y python3-dev" ;;
            "cmake") cmd="sudo apt-get install -y cmake" ;;
            "meson") cmd="sudo apt-get install -y meson" ;;
            "boost") cmd="sudo apt-get install -y libboost-all-dev" ;;
            "dialog") cmd="sudo apt-get install -y dialog" ;;
          esac
          ;;
        "fedora")
          case "$dep_name" in
            "compiler") cmd="sudo dnf install -y gcc-c++" ;;
            "python-dev") cmd="sudo dnf install -y python3-devel" ;;
            "cmake") cmd="sudo dnf install -y cmake" ;;
            "meson") cmd="sudo dnf install -y meson" ;;
            "boost") cmd="sudo dnf install -y boost-devel" ;;
            "dialog") cmd="sudo dnf install -y dialog" ;;
          esac
          ;;
        "arch"|"manjaro")
          case "$dep_name" in
            "compiler") cmd="sudo pacman -S --noconfirm base-devel" ;;
            "python-dev") cmd="sudo pacman -S --noconfirm python" ;;
            "cmake") cmd="sudo pacman -S --noconfirm cmake" ;;
            "meson") cmd="sudo pacman -S --noconfirm meson" ;;
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

    log "${BLUE}[Info] Running meson setup with options: ${meson_opts[*]}${NC}"
    if ! meson setup "${BUILD_DIR}" "${meson_opts[@]}" ${reconfigure_flag}; then
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
  declare -A DEP_STATUS
  check_compiler; DEP_STATUS[compiler]=$?
  check_python_dev; DEP_STATUS[python-dev]=$?
  check_cmake; DEP_STATUS[cmake]=$?
  check_meson; DEP_STATUS[meson]=$?
  check_boost; DEP_STATUS[boost]=$?

  local choices
  choices=$(dialog --clear --backtitle "Project Dependency Installer" \
    --title "Install System Dependencies" \
    --checklist "Select dependencies to install. Already found dependencies are unchecked." 20 70 5 \
    "compiler" "C++ Compiler (g++ or clang++)" "$([[ ${DEP_STATUS[compiler]} -ne 0 ]] && echo "on" || echo "off")" \
    "python-dev" "Python 3 Dev Headers" "$([[ ${DEP_STATUS[python-dev]} -ne 0 ]] && echo "on" || echo "off")" \
    "cmake" "CMake" "$([[ ${DEP_STATUS[cmake]} -ne 0 ]] && echo "on" || echo "off")" \
    "meson" "Meson Build System" "$([[ ${DEP_STATUS[meson]} -ne 0 ]] && echo "on" || echo "off")" \
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
}

run_python_bindings_tui() {
    local python_exec
    python_exec=$(command -v python3)
    if [ -z "$python_exec" ]; then
        dialog --msgbox "Could not find 'python3' executable. Please ensure Python 3 is installed and in your PATH." 8 60
        return
    fi

    local choice
    choice=$(dialog --clear --backtitle "Python Bindings Installer" \
        --title "Install Python Bindings" \
        --menu "Using Python: ${python_exec}\n\nSelect installation mode:" 15 70 2 \
        "1" "Developer Mode (pip install -e .)" \
        "2" "User Mode (pip install .)" \
        3>&1 1>&2 2>&3)

    clear
    case "$choice" in
        1)
            log "${BLUE}[Info] Installing Python bindings in Developer Mode...${NC}"
            if ! pip install -e . --no-build-isolation -vv; then
                log "${RED}[Error] Failed to install Python bindings in developer mode.${NC}"
                dialog --msgbox "Developer mode installation failed. Check the log for details." 8 60
            else
                log "${GREEN}[Success] Python bindings installed in developer mode.${NC}"
                dialog --msgbox "Successfully installed Python bindings in developer mode." 8 60
            fi
            ;;
        2)
            log "${BLUE}[Info] Installing Python bindings in User Mode...${NC}"
            if ! pip install .; then
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

run_build_config_tui() {
    local choice
    choice=$(dialog --clear --backtitle "Build Configuration" \
        --title "Configure Build Options" \
        --menu "Select an option to configure:" 20 70 6 \
        "1" "Build Directory (current: ${BUILD_DIR})" \
        "2" "Install Prefix (current: ${INSTALL_PREFIX})" \
        "3" "Build Type (current: ${MESON_BUILD_TYPE})" \
        "4" "Log Level (current: ${MESON_LOG_LEVEL})" \
        "5" "Generate pkg-config (current: ${MESON_PKG_CONFIG})" \
        "6" "Number of Cores (current: ${MESON_NUM_CORES})" \
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
            local build_type_choice
            build_type_choice=$(dialog --title "Select Build Type" --menu "" 15 70 3 \
                "release" "Optimized for performance" \
                "debug" "With debug symbols, no optimization" \
                "debugoptimized" "With debug symbols and optimization" \
                3>&1 1>&2 2>&3)
            if [ -n "$build_type_choice" ]; then
                MESON_BUILD_TYPE="$build_type_choice"
                log "${BLUE}[Config] Set build type to: ${MESON_BUILD_TYPE}${NC}"
            fi
            ;;
        4)
            local log_level_choice
            log_level_choice=$(dialog --title "Select Log Level" --menu "" 15 70 8 \
                "traceL3" "" "traceL2" "" "traceL1" "" "debug" "" "info" "" "warning" "" "error" "" "critical" "" \
                3>&1 1>&2 2>&3)
            if [ -n "$log_level_choice" ]; then
                MESON_LOG_LEVEL="$log_level_choice"
                log "${BLUE}[Config] Set log level to: ${MESON_LOG_LEVEL}${NC}"
            fi
            ;;
        5)
            if dialog --title "Generate pkg-config" --yesno "Generate gridfire.pc file?" 7 60; then
                MESON_PKG_CONFIG="true"
            else
                MESON_PKG_CONFIG="false"
            fi
            log "${BLUE}[Config] Set pkg-config generation to: ${MESON_PKG_CONFIG}${NC}"
            ;;
        6)
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

run_main_tui() {
    if ! check_dialog_installed; then return 1; fi
    while true; do
        local build_dir_status="Not Found"
        [ -d "$BUILD_DIR" ] && build_dir_status="Found"

        local choice
        choice=$(dialog --clear --backtitle "GridFire Installation and Build Manager" \
            --title "Main Menu" \
            --menu "DIR: ${BUILD_DIR} | TYPE: ${MESON_BUILD_TYPE} | CORES: ${MESON_NUM_CORES}\nPREFIX: ${INSTALL_PREFIX}\nLOG: ${MESON_LOG_LEVEL} | PKG-CONFIG: ${MESON_PKG_CONFIG}" 20 70 10 \
            "1" "Install System Dependencies" \
            "2" "Configure Build Options" \
            "3" "Install Python Bindings" \
            "4" "Run Full Build (Setup + Compile)" \
            "5" "Run Meson Setup/Reconfigure" \
            "6" "Run Meson Compile" \
            "7" "Run Meson Install (requires sudo)" \
            "8" "Run Tests" \
            "9" "Exit" \
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
            9) break ;;
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
  declare -A CHECKS=(
      [compiler]="check_compiler" [python-dev]="check_python_dev"
      [cmake]="check_cmake" [meson]="check_meson" [boost]="check_boost"
  )
  local all_deps_met=true
  for dep in "${!CHECKS[@]}"; do
      if ! ${CHECKS[$dep]}; then
          all_deps_met=false
          local install_cmd; install_cmd=$(get_install_cmd "$dep")
          if [ -n "$install_cmd" ]; then
              if prompt_yes_no "Dependency '${dep}' is missing. Attempt to install? (y/n):"; then
                  run_install_cmd "$install_cmd"
              fi
          else
              log "${RED}[Error] No automatic installation for '${dep}'. Please install manually.${NC}"
          fi
      fi
  done

  log "\n${BLUE}--- Re-checking all dependencies ---${NC}"
  local final_fail=false
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
