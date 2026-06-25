# pymmseqs/config/convert2fasta_config.py

from pathlib import Path
from typing import Union

from .base import BaseConfig
from ..defaults import loader
from ..utils import (
    get_caller_dir,
    run_mmseqs_command
)

DEFAULTS = loader.load("convert2fasta")


class Convert2FastaConfig(BaseConfig):
    """
    Convert an MMseqs2 sequence database back to a FASTA file.

    Parameters
    ----------
    `sequence_db` : Union[str, Path]
        Path to the input sequence database (created by createdb, extractorfs, etc.).

    `fasta_file` : Union[str, Path]
        Output path for the FASTA file.

    `use_header_file` : bool, optional
        Use the sequence header DB instead of the body to map the entry keys
        - False (default)
        - True

    `v` : int, optional
        Verbosity level of the output
        - 0: Quiet
        - 1: Errors only
        - 2: Errors and warnings
        - 3: Errors, warnings, and info (default)
    """

    def __init__(
        self,
        sequence_db: Union[str, Path],
        fasta_file: Union[str, Path],
        use_header_file: bool = False,
        v: int = 3,
    ):
        super().__init__()

        self.sequence_db = Path(sequence_db)
        self.fasta_file = Path(fasta_file)
        self.use_header_file = use_header_file
        self.v = v

        self._defaults = DEFAULTS
        self._path_params = [param for param, info in DEFAULTS.items() if info['type'] == 'path']
        self._caller_dir = get_caller_dir()

    def _validate(self) -> None:
        self._check_required_files()
        self._validate_choices()

    def run(self) -> None:
        self._resolve_all_path(self._caller_dir)
        self._validate()

        args = self._get_command_args("convert2fasta")
        mmseqs_output = run_mmseqs_command(args)

        self._handle_command_output(
            mmseqs_output=mmseqs_output,
            output_identifier="Convert to FASTA",
            output_path=str(self.fasta_file)
        )
