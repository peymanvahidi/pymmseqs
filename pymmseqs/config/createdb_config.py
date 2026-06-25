# pymmseqs/config/createdb_config.py

from pathlib import Path
from typing import Union, List

from .base import BaseConfig
from ..defaults import loader   
from ..utils import (
    get_caller_dir,
    run_mmseqs_command
)

DEFAULTS = loader.load("createdb")

class CreateDBConfig(BaseConfig):
    """
    Create a MMseqs2 database from a FASTA file and save it to the specified path prefix

    Parameters
    ----------
    `fasta_file` : Union[List[Union[str, Path]], Union[str, Path]]
    Path(s) to the input FASTA file(s). This can be:
        - A single string or Path object (e.g., `"input.fasta"` or `Path("input.fasta")`)
        - A list of strings or Path objects (e.g., `["input1.fasta", "input2.fasta"]` or 
          `[Path("input1.fasta"), Path("input2.fasta")]`)
        - A mixed list of strings and Path objects (e.g., `["input1.fasta", Path("input2.fasta")]`)

        If the paths are relative, they are resolved relative to the directory of the calling script
        If the paths are absolute, they are used as-is.

    `sequence_db` : Union[str, Path]
        Database path prefix, including the desired directory structure (e.g., `"output/dbs/mydb"`)
        - If the path is relative, it is resolved relative to the directory
        of the calling script
        - If the path is absolute, it is used as-is
        All necessary parent directories will be created automatically

    `dbtype` : int, optional
        Database type
        - 0: Auto-detect (default)
        - 1: Amino acid sequences
        - 2: Nucleotide sequences

    `shuffle` : bool, optional
        Shuffle the input database entries
        - True (default)
        - False

    `createdb_mode` : int, optional
        Database creation mode
        - 0: Copy data (default)
        - 1: Soft-link data and write a new index (only works with single-line FASTA/Q)
        - 2: GPU compatible db

    `id_offset` : int, optional
        Numeric ID offset in the index file
        - 0 (default)

    `mask` : int, optional
        Mask sequences in prefilter stage with tantan
        - 0: without low complexity masking
        - 1: with low complexity masking (default)

    `mask_prob` : float, optional
        Mask sequences if probability is above threshold
        - 0.9 (default)

    `mask_lower_case` : int, optional
        Whether lowercase letters are excluded from k-mer search
        - 0: include region (default)
        - 1: exclude region

    `mask_n_repeat` : int, optional
        Repeat letters that occur more than this threshold in a row
        - 0 (default)

    `compressed` : bool, optional
        Compress the output files
        - True
        - False (default)

    `gpu` : int, optional
        Use GPU (CUDA) if possible to build a GPU-compatible database
        - 0: off (default)
        - 1: on

    `v` : int, optional
        Verbosity level of the output
        - 0: Quiet
        - 1: Errors only
        - 2: Errors and warnings
        - 3: Errors, warnings, and info (default)

    `write_lookup` : bool, optional
        Create a `.lookup` file mapping internal IDs to FASTA IDs
        - True (default)
        - False

    `threads` : Union[str, int], optional
        Number of CPU cores to use
        - 'all' (default), converted to all available cores
    """

    def __init__(
        self,
        fasta_file: Union[List[Union[str, Path]], Union[str, Path]],
        sequence_db: Union[str, Path],
        dbtype: int = 0,
        shuffle: bool = True,
        createdb_mode: int = 0,
        id_offset: int = 0,
        mask: int = 1,
        mask_prob: float = 0.9,
        mask_lower_case: int = 0,
        mask_n_repeat: int = 0,
        compressed: bool = False,
        gpu: int = 0,
        v: int = 3,
        write_lookup: bool = True,
        threads: Union[str, int] = 'all'
    ):
        super().__init__()

        self.fasta_file = fasta_file if isinstance(fasta_file, list) else [fasta_file]
        self.fasta_file = [Path(f) for f in self.fasta_file]
        self.sequence_db = Path(sequence_db)
        self.dbtype = dbtype
        self.shuffle = shuffle
        self.createdb_mode = createdb_mode
        self.id_offset = id_offset
        self.mask = mask
        self.mask_prob = mask_prob
        self.mask_lower_case = mask_lower_case
        self.mask_n_repeat = mask_n_repeat
        self.compressed = compressed
        self.gpu = gpu
        self.v = v
        self.write_lookup = write_lookup
        self.threads = threads

        self._defaults = DEFAULTS
        self._path_params = [param for param, info in DEFAULTS.items() if info['type'] == 'path']
        self._caller_dir = get_caller_dir()

    def _validate(self) -> None:
        self._check_required_files()
        self._validate_choices()
            
        if self.id_offset < 0:
            raise ValueError(f"id_offset is {self.id_offset} but must be non-negative")

    def run(self) -> None:
        self._resolve_all_path(self._caller_dir)

        self._validate()
        
        args = self._get_command_args("createdb")
        mmseqs_output = run_mmseqs_command(args)
        
        self._handle_command_output(
            mmseqs_output=mmseqs_output,
            output_identifier="Database creation",
            output_path=str(self.sequence_db)
        )
