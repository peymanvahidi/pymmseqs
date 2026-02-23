# pymmseqs/config/createtaxdb_config.py

from pathlib import Path
from typing import Union

from .base import BaseConfig
from ..defaults import loader
from ..utils import (
    get_caller_dir,
    run_mmseqs_command
)

DEFAULTS = loader.load("createtaxdb")

class CreateTaxDBConfig(BaseConfig):
    """
    Create a taxonomy database for a MMseqs2 sequence database

    Parameters
    ----------
    sequence_db : Union[str, Path]
        Path to the input MMseqs2 sequence database
        - If the path is relative, it is resolved relative to the directory
        of the calling script
        - If the path is absolute, it is used as-is

    tmp_dir : Union[str, Path]
        Path to the temporary directory
        - If the path is relative, it is resolved relative to the directory
        of the calling script
        - If the path is absolute, it is used as-is
        All necessary parent directories will be created automatically

    ncbi_tax_dump : Optional[Union[str, Path]], optional
        Path to the NCBI taxonomy dump directory
        The tax dump can be downloaded from:
        `ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz`

    tax_mapping_file : Optional[Union[str, Path]], optional
        Path to the file mapping sequence identifiers to taxonomical identifiers

    tax_mapping_mode : int, optional
        Mode for mapping taxonomy based on sequence database
        - 0: .lookup file (default)
        - 1: .source file

    tax_db_mode : int, optional
        Mode for creating taxonomy database
        - 0: .dmp flat files (human readable)
        - 1: binary dump (faster reading) (default)

    threads : int, optional
        Number of CPU cores to use
        - Default: 14

    v : int, optional
        Verbosity level of the output
        - 0: quiet
        - 1: + errors
        - 2: + warnings
        - 3: + info (default)
    """

    def __init__(
        self,
        sequence_db: Union[str, Path],
        tmp_dir: Union[str, Path],
        ncbi_tax_dump: Union[str, Path] = "",
        tax_mapping_file: Union[str, Path] = "",
        tax_mapping_mode: int = 0,
        tax_db_mode: int = 1,
        threads: Union[str, int] = 'all',
        v: int = 3
    ):
        super().__init__()
        
        self.sequence_db = Path(sequence_db)
        self.tmp_dir = Path(tmp_dir)
        self.ncbi_tax_dump = Path(ncbi_tax_dump) if ncbi_tax_dump else ""
        self.tax_mapping_file = Path(tax_mapping_file) if tax_mapping_file else ""
        self.tax_mapping_mode = tax_mapping_mode
        self.tax_db_mode = tax_db_mode
        self.threads = threads
        self.v = v

        self._defaults = DEFAULTS
        self._path_params = [param for param, info in DEFAULTS.items() if info['type'] == 'path']
        self._caller_dir = get_caller_dir()

    def _validate(self) -> None:
        self._check_required_files()
        self._validate_choices()

        # Validate numeric constraints
        if self.threads < 1:
            raise ValueError(f"threads must be positive, got {self.threads}")
        
        if not (0 <= self.tax_mapping_mode <= 1):
            raise ValueError(f"tax_mapping_mode must be either 0 or 1, got {self.tax_mapping_mode}")
        
        if not (0 <= self.tax_db_mode <= 1):
            raise ValueError(f"tax_db_mode must be either 0 or 1, got {self.tax_db_mode}")
        
        if not (0 <= self.v <= 3):
            raise ValueError(f"verbosity level (v) must be between 0 and 3, got {self.v}")

    def run(self) -> None:
        self._resolve_all_path(self._caller_dir)

        self._validate()
        
        args = self._get_command_args("createtaxdb")
        mmseqs_output = run_mmseqs_command(args)
        
        self._handle_command_output(
            mmseqs_output=mmseqs_output,
            output_identifier="Create Taxonomy Database",
            output_path=str(self.sequence_db)
        )
