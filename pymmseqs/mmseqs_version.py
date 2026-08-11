# pymmseqs/mmseqs_version.py

# Single source of truth for the pinned MMseqs2 release.
#
# Read at build time by build_hook.py (via runpy, so this module must stay
# import-free) and at runtime by utils/binary.py for the download fallback.
# Bumping this here is enough; nothing else hardcodes the version.
MMSEQS_VERSION = "18-8cc5c"
