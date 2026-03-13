# build_hook.py
import os
import subprocess

target_dir = os.path.join("pymmseqs", "bin")
if not os.path.exists(os.path.join(target_dir, "mmseqs")):
    os.makedirs(target_dir, exist_ok=True)
    try:
        subprocess.check_call(["sh", "scripts/download_mmseqs.sh", target_dir])
    except subprocess.CalledProcessError as e:
        raise RuntimeError("Failed to download MMseqs binary") from e
