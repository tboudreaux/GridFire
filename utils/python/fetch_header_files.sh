#!/bin/bash

# --- Configuration ---
PYTHON_VERSIONS=("3.8.10" "3.9.13" "3.10.11" "3.11.9" "3.12.3" "3.13.0" "3.14.0")

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_OUTPUT_DIR="$SCRIPT_DIR/../../cross/python_includes"

# --- OS Detection ---
OS="$(uname -s)"
echo "Detected OS: $OS"

# --- Dependency Check ---
check_dependencies() {
    if [ "$OS" == "Linux" ]; then
        if ! command -v 7z &> /dev/null; then
            echo "Error: '7z' (p7zip-full) is required on Linux."
            exit 1
        fi
        if ! command -v cpio &> /dev/null; then
            echo "Error: 'cpio' is required."
            exit 1
        fi
    fi
}

# --- Extraction Logic (OS Specific) ---
extract_pkg() {
    local pkg_file="$1"
    local extract_root="$2"
    local major_ver="$3" # e.g., 3.11

    echo "  -> Extracting..."

    if [ "$OS" == "Darwin" ]; then
        pkgutil --expand "$pkg_file" "$extract_root/expanded"
        local payload_path="$extract_root/expanded/Python_Framework.pkg/Payload"
        
        if [ ! -f "$payload_path" ]; then
            echo "  -> Error: Could not find Payload in package."
            return 1
        fi
        mkdir -p "$extract_root/root"
        pushd "$extract_root/root" > /dev/null
        cat "$payload_path" | gunzip | cpio -id "*include/python${major_ver}/*" 2>/dev/null
        popd > /dev/null

    else
        7z x "$pkg_file" -o"$extract_root/expanded" -y > /dev/null
        local payload_path="$extract_root/expanded/Python_Framework.pkg/Payload"
        
        if [ ! -f "$payload_path" ]; then
            echo "  -> Error: Could not find Payload in package."
            return 1
        fi

        mkdir -p "$extract_root/root"
        pushd "$extract_root/root" > /dev/null
        cat "$payload_path" | gunzip | cpio -id "*include/python${major_ver}/*" 2>/dev/null
        popd > /dev/null
    fi
}

check_dependencies

mkdir -p "$BASE_OUTPUT_DIR"

for FULL_VER in "${PYTHON_VERSIONS[@]}"; do
    MAJOR_VER=$(echo "$FULL_VER" | cut -d. -f1,2)
    
    TARGET_DIR="$BASE_OUTPUT_DIR/python-$MAJOR_VER"
    TEMP_DIR="$BASE_OUTPUT_DIR/tmp_$FULL_VER"
    PKG_NAME="python-${FULL_VER}-macos11.pkg"
    
    if [[ "$MAJOR_VER" == "3.8" ]]; then
        PKG_NAME="python-${FULL_VER}-macosx10.9.pkg"
    fi

    DOWNLOAD_URL="https://www.python.org/ftp/python/${FULL_VER}/$PKG_NAME"

    echo "Processing Python $FULL_VER..."

    if [ -d "$TARGET_DIR" ] && [ "$(ls -A $TARGET_DIR)" ]; then
        echo "  -> Headers already exist in $TARGET_DIR. Skipping."
        continue
    fi

    mkdir -p "$TEMP_DIR"

    echo "  -> Downloading from $DOWNLOAD_URL"
    curl -L -s -o "$TEMP_DIR/python.pkg" "$DOWNLOAD_URL"
    
    if [ $? -ne 0 ]; then
        echo "  -> Download failed! Check version number or internet connection."
        rm -rf "$TEMP_DIR"
        continue
    fi

    # 2. Extract
    extract_pkg "$TEMP_DIR/python.pkg" "$TEMP_DIR" "$MAJOR_VER"

    # 3. Move Headers to Final Location
    # The cpio extraction usually results in: ./Versions/X.Y/include/pythonX.Y
    # We want to move that specific include folder to our target dir
    
    FOUND_HEADERS=$(find "$TEMP_DIR/root" -type d -path "*/include/python${MAJOR_VER}" | head -n 1)
    
    if [ -n "$FOUND_HEADERS" ]; then
        echo "  -> Found headers at: $FOUND_HEADERS"
        
        # Move the content to the final destination
        # We want the folder to be .../python-3.11/include/python3.11
        mkdir -p "$TARGET_DIR/include"
        mv "$FOUND_HEADERS" "$TARGET_DIR/include/"
        
        # Verify pyconfig.h exists (sanity check)
        if [ -f "$TARGET_DIR/include/python${MAJOR_VER}/pyconfig.h" ]; then
            echo "  -> Success: Headers installed to $TARGET_DIR"
        else 
            echo "  -> Warning: Header move seemed successful, but pyconfig.h is missing."
        fi
    else
        echo "  -> Error: Could not locate header files after extraction."
    fi

    # 4. Cleanup
    rm -rf "$TEMP_DIR"
    echo "---------------------------------------------------"
done

echo "Done. All headers stored in $BASE_OUTPUT_DIR"

