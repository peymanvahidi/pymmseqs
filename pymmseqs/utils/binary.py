# pymmseqs/utils/binary.py
import os
import platform
import shutil
import tarfile
import urllib.request

from ..mmseqs_version import MMSEQS_VERSION

_RELEASE_URL = (
    "https://github.com/soedinglab/mmseqs2/releases/download/{version}/mmseqs-{target}.tar.gz"
)

# MMseqs2 release asset per platform.
#
# ponytail: the x86_64 build is AVX2, so it needs a Haswell (2013) or newer CPU
# and dies with SIGILL on anything older. A wheel cannot express a CPU-feature
# requirement, so this is a documented ceiling rather than something we detect;
# affected users point MMSEQS2_PATH at an sse41 build. Upgrade path if that ever
# shows up in practice: read /proc/cpuinfo here and fall back to linux-sse41.
_TARGETS = {
    ("Linux", "x86_64"): "linux-avx2",
    ("Linux", "aarch64"): "linux-arm64",
    ("Linux", "arm64"): "linux-arm64",
    ("Darwin", "x86_64"): "osx-universal",
    ("Darwin", "arm64"): "osx-universal",
}


def _cache_dir():
    """
    Directory for binaries downloaded at runtime.

    Honours XDG_CACHE_HOME, which is the standard lever for HPC users whose
    $HOME is small or read-only. Keyed by version so a bump does not silently
    reuse the old binary.
    """
    root = os.environ.get("XDG_CACHE_HOME") or os.path.join(
        os.path.expanduser("~"), ".cache"
    )
    return os.path.join(root, "pymmseqs", MMSEQS_VERSION)


def _download_mmseqs(dest):
    """
    Download the pinned MMseqs2 release and extract the binary to `dest`.

    Returns the path to the extracted binary.
    """
    target = _TARGETS.get((platform.system(), platform.machine()))
    if target is None:
        raise FileNotFoundError(
            f"No MMseqs2 build is available for {platform.system()} "
            f"{platform.machine()}. Install mmseqs2 yourself and point the "
            f"MMSEQS2_PATH environment variable at the binary."
        )

    url = _RELEASE_URL.format(version=MMSEQS_VERSION, target=target)
    os.makedirs(os.path.dirname(dest), exist_ok=True)

    # Stream into <dest>.part and rename only once complete. An interrupted
    # download would otherwise leave a truncated file that still satisfies
    # os.path.exists() and gets reused forever.
    part = dest + ".part"
    try:
        with urllib.request.urlopen(url) as response:
            # Stream mode ("r|gz") so we never buffer the whole tarball, and
            # extractfile() rather than extract() so archive member names can
            # never influence the path we write to.
            with tarfile.open(fileobj=response, mode="r|gz") as tar:
                for member in tar:
                    if not member.isfile() or not member.name.endswith("bin/mmseqs"):
                        continue
                    source = tar.extractfile(member)
                    if source is None:
                        continue
                    with open(part, "wb") as out:
                        shutil.copyfileobj(source, out)
                    break
                else:
                    raise FileNotFoundError(f"No mmseqs binary inside {url}")
    except OSError as e:
        if os.path.exists(part):
            os.remove(part)
        raise RuntimeError(f"Failed to download mmseqs2 from {url}") from e

    os.chmod(part, 0o755)
    os.replace(part, dest)
    return dest


def get_mmseqs_binary():
    """
    Retrieve the path to the mmseqs2 binary.

    Resolution order:
        1. The MMSEQS2_PATH environment variable, if set.
        2. The binary bundled alongside this package by the wheel.
        3. A copy downloaded earlier into the user cache.
        4. A fresh download of the pinned release into the user cache.
    """
    custom_path = os.getenv('MMSEQS2_PATH')
    if custom_path:
        if os.path.exists(custom_path):
            return custom_path
        else:
            raise FileNotFoundError(
                f"mmseqs2 binary specified by MMSEQS2_PATH does not exist: {custom_path}"
            )

    system = platform.system()
    binary_name = 'mmseqs.exe' if system == 'Windows' else 'mmseqs'

    # Resolve relative to this module so the bundled binary is found in every
    # install layout (venv, --user, conda, editable), not just global purelib.
    package_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    bundled_path = os.path.join(package_dir, 'bin', binary_name)
    if os.path.exists(bundled_path):
        return bundled_path

    # No bundled binary: either an sdist install on a platform we ship no wheel
    # for, or a wheel whose binary was stripped. Fall back to the cache.
    cached_path = os.path.join(_cache_dir(), binary_name)
    if os.path.exists(cached_path):
        return cached_path

    return _download_mmseqs(cached_path)
