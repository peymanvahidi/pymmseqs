#!/bin/sh
# scripts/download_mmseqs.sh

# Exit on error and show commands
set -ex

# Configuration
MMSEQS_VERSION="18-8cc5c"
TARGET_DIR="$1"  # Accept target directory as the first argument
BASE_URL="https://github.com/soedinglab/mmseqs2/releases/download/${MMSEQS_VERSION}"

# Check if TARGET_DIR is provided
if [ -z "${TARGET_DIR}" ]; then
    echo "Error: Target directory must be provided as the first argument"
    exit 1
fi

# Create target directory
mkdir -p "${TARGET_DIR}"

# OS-specific download and extraction
case "$(uname -s)" in
    Linux*)
        ARCH=$(uname -m)
        if [ "${ARCH}" = "x86_64" ]; then
            URL="${BASE_URL}/mmseqs-linux-avx2.tar.gz"
        elif [ "${ARCH}" = "aarch64" ]; then
            URL="${BASE_URL}/mmseqs-linux-arm64.tar.gz"
        elif [ "${ARCH}" = "i686" ]; then
            URL="${BASE_URL}/mmseqs-linux-sse2.tar.gz"
        else
            echo "Unsupported Linux architecture: ${ARCH}"
            exit 1
        fi
        
        # Download and extract
        curl -L "${URL}" | tar -zxf - \
            --strip-components=2 \
            -C "${TARGET_DIR}" \
            "mmseqs/bin/mmseqs"
        ;;
    
    Darwin*)
        # macOS universal binary
        curl -L "${BASE_URL}/mmseqs-osx-universal.tar.gz" | tar -zxf - \
            --strip-components=2 \
            -C "${TARGET_DIR}" \
            "mmseqs/bin/mmseqs"
        ;;
    
    # MINGW*|CYGWIN*|MSYS*)
    #     # Windows
    #     curl -L "${BASE_URL}/mmseqs-win64.zip" -o mmseqs.zip
    #     unzip -j mmseqs.zip "mmseqs/bin/mmseqs.exe" -d "${TARGET_DIR}"
    #     rm mmseqs.zip
    #     ;;
    
    *)
        echo "Unsupported operating system"
        exit 1
        ;;
esac

# Verify binary exists
ls -l "${TARGET_DIR}/mmseqs"*

# Set executable permissions (Unix systems)
if [ "$(uname -s)" != "MINGW"* ]; then
    chmod +x "${TARGET_DIR}/mmseqs"
fi