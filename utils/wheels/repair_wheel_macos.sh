#!/bin/zsh

set -e
# Color codes for output
RED='\033[0;31m'
YELLOW='\033[1;33m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

function fix_file_rpaths() {
    local file_path="$1"
    echo -e "${YELLOW}Fixing RPATHs in file: $file_path...${NC}"
    python3 "$FIX_RPATH_SCRIPT" "$file_path"
    if [ $? -ne 0 ]; then
        echo -e "${RED}Error: RPATH fix script failed for file: $file_path${NC}"
        exit 1
    fi
    echo -e "${GREEN}RPATHs fixed for file: $file_path${NC}"
}

export -f fix_file_rpaths

echo -e "${YELLOW}"
echo "========================================================================="
echo "                    TEMPORARY WHEEL REPAIR WORKAROUND"
echo "========================================================================="
echo -e "${NC}"
echo ""
echo -e "${YELLOW}WARNING:${NC} This script applies a temporary patch to fix"
echo "a known issue with meson-python that causes duplicate RPATH entries in"
echo "built Python wheels on macOS, preventing module imports."
echo ""
echo "This workaround will:"
echo "  1. Unzip the wheel file"
echo "  2. Locate the extension modules"
echo "  3. Remove duplicate RPATH entries using install_name_tool"
echo "  4. Resign the wheel if necessary"
echo "  5. Repackage the wheel file"
echo ""

FIX_RPATH_SCRIPT="../../build-python/fix_rpaths.py"

# get the wheel directory to scan through
WHEEL_DIR="$1"
if [ -z "$WHEEL_DIR" ]; then
    echo -e "${RED}Error: No wheel directory specified.${NC}"
    echo "Usage: $0 /path/to/wheel_directory"
    exit 1
fi

REPAIRED_WHEELS_DIR="repaired_wheels"
mkdir -p "$REPAIRED_WHEELS_DIR"

REPAIRED_DELOCATED_WHEELS_DIR="${REPAIRED_WHEELS_DIR}/delocated"

# Scal all files ending in .whl and not starting with a dot
for WHEEL_PATH in "$WHEEL_DIR"/*.whl; do
    if [ ! -f "$WHEEL_PATH" ]; then
        echo -e "${YELLOW}No wheel files found in directory: $WHEEL_DIR${NC}"
        exit 0
    fi
    echo ""
    echo -e "${GREEN}Processing wheel: $WHEEL_PATH${NC}"

    WHEEL_NAME=$(basename "$WHEEL_PATH")
    TEMP_DIR=$(mktemp -d)

    echo -e "${GREEN}Step 1: Unzipping wheel...${NC}"
    python -m wheel unpack "$WHEEL_PATH" -d "$TEMP_DIR"

    echo -e "${GREEN}Step 2: Locating extension modules...${NC}"
    while IFS= read -r -d '' so_file; do
        echo "Found library: $so_file"
        fix_file_rpaths "$so_file"
    done < <(find "$TEMP_DIR" -name "*.so" -print0)

    echo -e "${GREEN}Step 4: Repackaging wheel...${NC}"
    python -m wheel pack "$TEMP_DIR/gridfire-0.7.4rc2" -d "$REPAIRED_WHEELS_DIR"

    REPAIRED_WHEEL_PATH="${REPAIRED_WHEELS_DIR}/${WHEEL_NAME}"

    echo -e "${GREEN}Step 5: Delocating wheel...${NC}"
    # Ensure delocate is installed
    pip install delocate
    delocate-wheel -w "$REPAIRED_DELOCATED_WHEELS_DIR" "$REPAIRED_WHEEL_PATH"

    echo -e "${GREEN}Repaired wheel saved to: ${REPAIRED_DELOCATED_WHEELS_DIR}/${WHEEL_NAME}${NC}"
    # Clean up temporary directory
    rm -rf "$TEMP_DIR"
done