# pymmseqs/utils/runner.py

import subprocess
from .binary import get_mmseqs_binary

def run_mmseqs_command(args, capture_output=True):
    """
    Run an mmseqs2 command with the given arguments.
    Raises RuntimeError if the command fails.
    Returns the command's result.
    """
    binary = get_mmseqs_binary()
    cmd = [binary] + args

    print("\n" + " Running pymmseqs command... ")

    result = subprocess.run(cmd, capture_output=capture_output, text=True)
    
    return result
