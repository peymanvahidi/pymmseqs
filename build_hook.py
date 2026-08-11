# build_hook.py
import os
import runpy
import subprocess

# Read the pin without importing the package - pymmseqs itself is not
# importable yet at build time, but mmseqs_version.py has no imports.
MMSEQS_VERSION = runpy.run_path(
    os.path.join("pymmseqs", "mmseqs_version.py")
)["MMSEQS_VERSION"]

target_dir = os.path.join("pymmseqs", "bin")
if not os.path.exists(os.path.join(target_dir, "mmseqs")):
    os.makedirs(target_dir, exist_ok=True)
    try:
        subprocess.check_call(
            ["sh", "scripts/download_mmseqs.sh", target_dir],
            env={**os.environ, "MMSEQS_VERSION": MMSEQS_VERSION},
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError("Failed to download MMseqs binary") from e
