# pymmseqs/utils/binary.py
import os
import platform

def get_mmseqs_binary():
    """
    Retrieve the path to the mmseqs2 binary.
    Allows overriding via the MMSEQS2_PATH environment variable.
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
    binary_path = os.path.join(package_dir, 'bin', binary_name)
    
    if not os.path.exists(binary_path):
        raise FileNotFoundError(
            f"mmseqs2 binary not found at {binary_path}. Please ensure it is installed correctly."
        )
    
    return binary_path
