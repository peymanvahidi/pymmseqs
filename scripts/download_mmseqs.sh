#!/bin/sh
# scripts/download_mmseqs.sh
#
# Usage:
#   MMSEQS_VERSION=18-8cc5c [MMSEQS_TARGET=linux-avx2] sh scripts/download_mmseqs.sh <target_dir>
#
# MMSEQS_TARGET names an MMseqs2 release asset (they are all mmseqs-<target>.tar.gz).
# Leave it unset to detect it from the host, which is what build_hook.py does.
# CI sets it explicitly so a single Linux runner can produce the wheels for every
# platform: nothing here compiles, so there is nothing to cross-compile - only a
# different tarball to download and a different platform tag to stamp on.

# Exit on error and show commands
set -ex

TARGET_DIR="$1"
: "${TARGET_DIR:?Target directory must be provided as the first argument}"
: "${MMSEQS_VERSION:?must be set by the caller, see pymmseqs/mmseqs_version.py}"

BASE_URL="https://github.com/soedinglab/mmseqs2/releases/download/${MMSEQS_VERSION}"

# Detect the asset from the host unless the caller pinned one
if [ -z "${MMSEQS_TARGET:-}" ]; then
    case "$(uname -s)" in
        Linux*)
            case "$(uname -m)" in
                x86_64)  MMSEQS_TARGET="linux-avx2" ;;
                aarch64) MMSEQS_TARGET="linux-arm64" ;;
                i686)    MMSEQS_TARGET="linux-sse2" ;;
                *)
                    echo "Unsupported Linux architecture: $(uname -m)"
                    exit 1
                    ;;
            esac
            ;;
        Darwin*)
            # Universal binary, covers both arm64 and x86_64
            MMSEQS_TARGET="osx-universal"
            ;;
        *)
            echo "Unsupported operating system: $(uname -s)"
            exit 1
            ;;
    esac
fi

mkdir -p "${TARGET_DIR}"

curl -fL "${BASE_URL}/mmseqs-${MMSEQS_TARGET}.tar.gz" | tar -zxf - \
    --strip-components=2 \
    -C "${TARGET_DIR}" \
    "mmseqs/bin/mmseqs"

chmod +x "${TARGET_DIR}/mmseqs"
ls -l "${TARGET_DIR}/mmseqs"
