#!/bin/bash
# pip_install_mac_patch.sh - Workaround for meson-python duplicate RPATH bug on macOS

set -e

# Color codes for output
RED='\033[0;31m'
YELLOW='\033[1;33m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

# Returns 0 if the Mach-O binary at $1 has duplicate LC_RPATH entries.
has_duplicate_rpaths() {
    local binary="$1"
    local rpaths dup
    rpaths=$(otool -l "$binary" | awk '/cmd LC_RPATH/{getline; getline; print $2}')
    dup=$(printf '%s\n' "$rpaths" | sort | uniq -d)
    [ -n "$dup" ]
}

echo -e "${YELLOW}"
echo "========================================================================="
echo "          INSTALLATION + DUPLICATE-RPATH SAFETY NET (macOS)"
echo "========================================================================="
echo -e "${NC}"
echo ""
echo "This script installs gridfire with pip and then checks the installed"
echo "extension modules for duplicate LC_RPATH entries (a meson-python bug"
echo "exposed by macOS 26.1, see:"
echo "  https://github.com/mesonbuild/meson-python/issues/813 )."
echo ""
echo "With the current self-contained wheel layout the bug should not"
echo "trigger; binaries are only patched if duplicates are actually found."
echo ""
echo -e "${YELLOW}Continue? [y/N]${NC} "
read -r response

if [[ ! "$response" =~ ^[Yy]$ ]]; then
    echo -e "${RED}Installation cancelled.${NC}"
    exit 1
fi

echo ""
echo -e "${GREEN}Step 1: Finding current Python environment...${NC}"

PYTHON_BIN=$(which python3)
if [ -z "$PYTHON_BIN" ]; then
    echo -e "${RED}Error: python3 not found in PATH${NC}"
    exit 1
fi

echo "Using Python: $PYTHON_BIN"
PYTHON_VERSION=$($PYTHON_BIN --version)
echo "Python version: $PYTHON_VERSION"

SITE_PACKAGES=$($PYTHON_BIN -c "import site; print(site.getsitepackages()[0])")
echo "Site packages: $SITE_PACKAGES"
echo ""

echo -e "${GREEN}Step 2: Installing gridfire with pip...${NC}"
$PYTHON_BIN -m pip install . -v --no-build-isolation

if [ $? -ne 0 ]; then
    echo -e "${RED}Error: pip install failed${NC}"
    exit 1
fi
echo ""

FIX_SCRIPT="build-python/fix_rpaths.py"

check_and_fix() {
    local label="$1" so_file="$2"

    if [ -z "$so_file" ]; then
        echo -e "${YELLOW}Skipping ${label}: extension module not found (package may not be installed).${NC}"
        return 0
    fi

    echo "Found ${label} extension module: $so_file"
    if has_duplicate_rpaths "$so_file"; then
        echo -e "${YELLOW}Duplicate LC_RPATH entries detected in ${label}; applying fix...${NC}"
        if [ ! -f "$FIX_SCRIPT" ]; then
            echo -e "${RED}Error: $FIX_SCRIPT not found${NC}"
            echo "Please run this script from the project root directory."
            exit 1
        fi
        $PYTHON_BIN "$FIX_SCRIPT" "$so_file"
    else
        echo -e "${GREEN}No duplicate LC_RPATH entries in ${label}; no patch needed.${NC}"
    fi
    echo ""
}

echo -e "${GREEN}Step 3: Checking installed extension modules...${NC}"

GRIDFIRE_SO=$(find "$SITE_PACKAGES/gridfire" -name "_gridfire.cpython-*-darwin.so" 2>/dev/null | head -n 1)
check_and_fix "gridfire" "$GRIDFIRE_SO"

FOURDST_SO=$(find "$SITE_PACKAGES/fourdst" -name "_phys.cpython-*-darwin.so" 2>/dev/null | head -n 1)
check_and_fix "fourdst" "$FOURDST_SO"

echo -e "${GREEN}=========================================================================${NC}"
echo -e "${GREEN}                    Installation Complete!${NC}"
echo -e "${GREEN}=========================================================================${NC}"
echo ""
echo "Test the installation with:"
echo "  $PYTHON_BIN -c 'import gridfire; print(gridfire.__version__)'"
echo ""
