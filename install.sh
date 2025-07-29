#!/usr/bin/env bash
#
# install.sh - Comprehensive interactive installation script for GridFire
#
# This script performs the following actions:
# 1. Checks for essential system dependencies:
#    - A C++ compiler (g++ or clang++)
#    - Python 3 development headers
#    - CMake
#    - Meson
#    - Boost libraries
# 2. If run with the --tui flag, it provides a text-based user interface
#    to select and install missing dependencies.
# 3. If run without flags, it prompts the user interactively for each missing dependency.
# 4. Provides detailed installation instructions for various Linux distributions and macOS.
# 5. Once all dependencies are met, it runs the meson setup and build commands.
# 6. Logs all operations to a file for easy debugging.

set -o pipefail

# --- Configuration ---
LOGFILE="GridFire_Installer.log"
BUILD_DIR="build"

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
  echo "  --tui         Run in Text-based User Interface mode for interactive dependency installation."
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

# Check for a C++ compiler.
check_compiler() {
  if command -v g++ &>/dev/null; then
    COMPILER_CMD="g++"
    log "${GREEN}[OK] Found C++ compiler: $(g++ --version | head -n1)${NC}"
    return 0
  elif command -v clang++ &>/dev/null; then
    COMPILER_CMD="clang++"
    log "${GREEN}[OK] Found C++ compiler: $(clang++ --version | head -n1)${NC}"
    return 0
  else
    log "${RED}[FAIL] No C++ compiler (g++ or clang++) found.${NC}"
    return 1
  fi
}

# Check for Python 3 development headers.
check_python_dev() {
  if python3-config --includes &>/dev/null; then
    log "${GREEN}[OK] Found Python 3 development headers.${NC}"
    return 0
  else
    log "${RED}[FAIL] Python 3 development headers not found.${NC}"
    return 1
  fi
}

# Check for CMake.
check_cmake() {
  if command -v cmake &>/dev/null; then
    log "${GREEN}[OK] Found CMake: $(cmake --version | head -n1)${NC}"
    return 0
  else
    log "${RED}[FAIL] CMake not found.${NC}"
    return 1
  fi
}

# Check for Meson.
check_meson() {
  if command -v meson &>/dev/null; then
    log "${GREEN}[OK] Found Meson: $(meson --version)${NC}"
    return 0
  else
    log "${RED}[FAIL] Meson not found.${NC}"
    return 1
  fi
}

# Check if Boost is installed and available to Meson.
check_boost() {
  log "${BLUE}[Info] Checking for Boost dependency using Meson...${NC}"
  if meson introspect --dependencies -C "$BUILD_DIR" 2>/dev/null | grep -q "boost.*found: YES"; then
      log "${GREEN}[OK] Boost dependency is already met in the Meson build directory.${NC}"
      return 0
  fi

  # A simple meson dependency check as a fallback
  local test_dir="meson-boost-test"
  rm -rf "$test_dir" && mkdir "$test_dir"
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

# Execute an installation command.
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

# Get installation command for a given dependency.
get_install_cmd() {
  local dep_name="$1"
  local cmd=""
  case "$OS_NAME" in
    "macOS")
      case "$dep_name" in
        "compiler") cmd="xcode-select --install" ;;
        "python-dev") cmd="brew install python3" ;;
        "cmake") cmd="brew install cmake" ;;
        "meson") cmd="brew install meson" ;;
        "boost") cmd="brew install boost" ;;
        "dialog") cmd="brew install dialog" ;;
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
        *)
          log "${YELLOW}[Warn] Unsupported Linux distribution: ${DISTRO_ID}. Cannot provide install commands.${NC}"
          ;;
      esac
      ;;
    *)
      log "${YELLOW}[Warn] Unsupported OS: ${OS_NAME}. Cannot provide install commands.${NC}"
      ;;
  esac
  echo "$cmd"
}

# --- Main Logic ---

# Run the Meson setup and build process.
run_meson_build() {
  log "\n${BLUE}--- Starting Meson Build ---${NC}"
  if [ ! -f "meson.build" ]; then
    log "${RED}[FATAL] meson.build file not found in the current directory. Cannot proceed.${NC}"
    exit 1
  fi

  log "${BLUE}[Info] Running meson setup...${NC}"
  if ! meson setup "${BUILD_DIR}"; then
    log "${RED}[FATAL] Meson setup failed. Check the output above and in ${LOGFILE} for details.${NC}"
    exit 1
  fi
  log "${GREEN}[Success] Meson setup complete.${NC}"

  log "${BLUE}[Info] Running meson compile...${NC}"
  if ! meson compile -C "${BUILD_DIR}"; then
    log "${RED}[FATAL] Meson compile failed. Check the output above and in ${LOGFILE} for details.${NC}"
    exit 1
  fi

  log "\n${GREEN}--- Build Successful! ---${NC}"
  log "${GREEN}The project has been built in the '${BUILD_DIR}' directory.${NC}"
}

# Run the TUI for dependency installation.
run_tui() {
  if ! command -v dialog &>/dev/null; then
    log "${YELLOW}[Warn] The 'dialog' utility is not installed. It is required for TUI mode.${NC}"
    local install_cmd
    install_cmd=$(get_install_cmd "dialog")
    if [ -n "$install_cmd" ]; then
      if prompt_yes_no "Would you like to install it now? (y/n):"; then
        eval "$install_cmd"
        if ! command -v dialog &>/dev/null; then
          log "${RED}[FATAL] Failed to install 'dialog'. Cannot run in TUI mode.${NC}"
          exit 1
        fi
      else
        log "${RED}[FATAL] Cannot run in TUI mode without 'dialog'. Exiting.${NC}"
        exit 1
      fi
    else
      log "${RED}[FATAL] Cannot automatically install 'dialog' on your system. Please install it manually and try again.${NC}"
      exit 1
    fi
  fi

  local choices
  choices=$(dialog --clear --backtitle "Project Dependency Installer" \
    --title "Install Dependencies" \
    --checklist "Select dependencies to install. Dependencies marked with (*) are missing." 20 70 5 \
    "compiler" "C++ Compiler (g++ or clang++)" "$([[ ${DEP_STATUS[compiler]} == "FAIL" ]] && echo "on" || echo "off")" \
    "python-dev" "Python 3 Dev Headers" "$([[ ${DEP_STATUS[python-dev]} == "FAIL" ]] && echo "on" || echo "off")" \
    "cmake" "CMake" "$([[ ${DEP_STATUS[cmake]} == "FAIL" ]] && echo "on" || echo "off")" \
    "meson" "Meson Build System" "$([[ ${DEP_STATUS[meson]} == "FAIL" ]] && echo "on" || echo "off")" \
    "boost" "Boost Libraries" "$([[ ${DEP_STATUS[boost]} == "FAIL" ]] && echo "on" || echo "off")" \
    3>&1 1>&2 2>&3)

  clear

  if [ -z "$choices" ]; then
    log "${YELLOW}[Info] No dependencies selected for installation.${NC}"
    return
  fi

  for choice in $choices; do
    local dep
    dep=$(echo "$choice" | tr -d '"')
    log "\n${BLUE}--- Installing ${dep} ---${NC}"
    local install_cmd
    install_cmd=$(get_install_cmd "$dep")
    if [ -n "$install_cmd" ]; then
      run_install_cmd "$install_cmd"
    else
      log "${RED}[Error] No automatic installation command available for '${dep}' on your system.${NC}"
      log "${RED}Please install it manually.${NC}"
    fi
  done
}

# --- Script Entry Point ---
main() {
  # Handle command-line arguments
  if [[ " $@ " =~ " --help " ]] || [[ " $@ " =~ " -h " ]]; then
    show_help
    exit 0
  fi

  if [[ " $@ " =~ " --clean " ]]; then
    log "${BLUE}[Info] Cleaning up previous build...${NC}"
    rm -rf "$BUILD_DIR" "$LOGFILE"
  fi

  # Start logging
  echo "" > "$LOGFILE"
  log "--- Project Installation Log ---"
  log "Date: $(date)"
  log "OS: ${OS_NAME}, Distro: ${DISTRO_ID}"

  # Perform all dependency checks
  log "\n${BLUE}--- Checking System Dependencies ---${NC}"
  declare -A DEP_STATUS
  check_compiler; DEP_STATUS[compiler]=$?
  check_python_dev; DEP_STATUS[python-dev]=$?
  check_cmake; DEP_STATUS[cmake]=$?
  check_meson; DEP_STATUS[meson]=$?
  check_boost; DEP_STATUS[boost]=$?

  local all_deps_met=0
  for status in "${DEP_STATUS[@]}"; do
    if [ "$status" -ne 0 ]; then
      all_deps_met=1
      break
    fi
  done

  if [ $all_deps_met -eq 0 ]; then
    log "\n${GREEN}--- All dependencies are met! ---${NC}"
    run_meson_build
    exit 0
  fi

  log "\n${YELLOW}--- Some dependencies are missing ---${NC}"

  if [[ " $@ " =~ " --tui " ]]; then
    run_tui
  else
    # Interactive command-line prompts
    for dep in "${!DEP_STATUS[@]}"; do
      if [ "${DEP_STATUS[$dep]}" -ne 0 ]; then
        local install_cmd
        install_cmd=$(get_install_cmd "$dep")
        if [ -n "$install_cmd" ]; then
          if prompt_yes_no "Dependency '${dep}' is missing. Attempt to install it now? (y/n):"; then
            run_install_cmd "$install_cmd"
          fi
        else
          log "${RED}[Error] No automatic installation command available for '${dep}' on your system. Please install it manually.${NC}"
        fi
      fi
    done
  fi

  # Final check before build
  log "\n${BLUE}--- Re-checking dependencies after installation attempts ---${NC}"
  check_compiler; DEP_STATUS[compiler]=$?
  check_python_dev; DEP_STATUS[python-dev]=$?
  check_cmake; DEP_STATUS[cmake]=$?
  check_meson; DEP_STATUS[meson]=$?
  check_boost; DEP_STATUS[boost]=$?

  all_deps_met=0
  local failed_deps=""
  for dep in "${!DEP_STATUS[@]}"; do
    if [ "${DEP_STATUS[$dep]}" -ne 0 ]; then
      all_deps_met=1
      failed_deps+="${dep} "
    fi
  done

  if [ $all_deps_met -eq 0 ]; then
    log "\n${GREEN}--- All dependencies are now met! ---${NC}"
    run_meson_build
  else
    log "\n${RED}[FATAL] The following dependencies are still missing: ${failed_deps}${NC}"
    log "${RED}Please install them manually and re-run this script.${NC}"
    log "${RED}Check ${LOGFILE} for detailed error messages.${NC}"
    exit 1
  fi
}

# Run the main function
main "$@"
