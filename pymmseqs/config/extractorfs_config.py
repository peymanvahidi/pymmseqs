# pymmseqs/config/extractorfs_config.py

from pathlib import Path
from typing import Union

from .base import BaseConfig
from ..defaults import loader
from ..utils import (
    get_caller_dir,
    run_mmseqs_command
)

DEFAULTS = loader.load("extractorfs")


class ExtractOrfsConfig(BaseConfig):
    """
    Six-frame extraction of open reading frames (ORFs) from a nucleotide
    sequence database using MMseqs2 extractorfs.

    Parameters
    ----------
    `sequence_db` : Union[str, Path]
        Path to the input (nucleotide) sequence database created with createdb.

    `orf_db` : Union[str, Path]
        Output path for the ORF sequence database.

    `min_length` : int, optional
        Minimum codon number in open reading frames
        - 30 (default)

    `max_length` : int, optional
        Maximum codon number in open reading frames
        - 32734 (default)

    `max_gaps` : int, optional
        Maximum number of codons with gaps or unknown residues before an ORF is rejected
        - 2147483647 (default)

    `contig_start_mode` : int, optional
        Contig start handling
        - 0: incomplete
        - 1: complete
        - 2: both (default)

    `contig_end_mode` : int, optional
        Contig end handling
        - 0: incomplete
        - 1: complete
        - 2: both (default)

    `orf_start_mode` : int, optional
        ORF fragment handling
        - 0: from start to stop
        - 1: from any to stop (default)
        - 2: from last encountered start to stop (no start in the middle)

    `forward_frames` : str, optional
        Comma-separated list of frames on the forward strand to be extracted
        - "1,2,3" (default)

    `reverse_frames` : str, optional
        Comma-separated list of frames on the reverse strand to be extracted
        - "1,2,3" (default)

    `translation_table` : int, optional
        Genetic code table to use
        - 1: Canonical (default)
        - 11: Prokaryote, etc. (see MMseqs2 docs for the full list)

    `translate` : bool, optional
        Translate ORF to amino acids
        - False (default)
        - True

    `use_all_table_starts` : bool, optional
        Use all alternative start codons in the genetic table
        - False: only ATG (AUG) (default)
        - True

    `id_offset` : int, optional
        Numeric IDs in index file are offset by this value
        - 0 (default)

    `threads` : Union[str, int], optional
        Number of CPU cores to use
        - 'all' (default), converted to all available cores

    `compressed` : bool, optional
        Write compressed output
        - False (default)

    `v` : int, optional
        Verbosity level of the output
        - 3 (default)

    `create_lookup` : bool, optional
        Create database lookup file (can be very large)
        - False (default)
    """

    def __init__(
        self,
        sequence_db: Union[str, Path],
        orf_db: Union[str, Path],
        min_length: int = 30,
        max_length: int = 32734,
        max_gaps: int = 2147483647,
        contig_start_mode: int = 2,
        contig_end_mode: int = 2,
        orf_start_mode: int = 1,
        forward_frames: str = "1,2,3",
        reverse_frames: str = "1,2,3",
        translation_table: int = 1,
        translate: bool = False,
        use_all_table_starts: bool = False,
        id_offset: int = 0,
        threads: Union[str, int] = 'all',
        compressed: bool = False,
        v: int = 3,
        create_lookup: bool = False,
    ):
        super().__init__()

        self.sequence_db = Path(sequence_db)
        self.orf_db = Path(orf_db)
        self.min_length = min_length
        self.max_length = max_length
        self.max_gaps = max_gaps
        self.contig_start_mode = contig_start_mode
        self.contig_end_mode = contig_end_mode
        self.orf_start_mode = orf_start_mode
        self.forward_frames = forward_frames
        self.reverse_frames = reverse_frames
        self.translation_table = translation_table
        self.translate = translate
        self.use_all_table_starts = use_all_table_starts
        self.id_offset = id_offset
        self.threads = threads
        self.compressed = compressed
        self.v = v
        self.create_lookup = create_lookup

        self._defaults = DEFAULTS
        self._path_params = [param for param, info in DEFAULTS.items() if info['type'] == 'path']
        self._caller_dir = get_caller_dir()

    def _validate(self) -> None:
        self._check_required_files()
        self._validate_choices()

        if not (0 <= self.min_length <= self.max_length):
            raise ValueError("min_length must be >= 0 and <= max_length")
        if self.id_offset < 0:
            raise ValueError(f"id_offset is {self.id_offset} but must be non-negative")

    def run(self) -> None:
        self._resolve_all_path(self._caller_dir)
        self._validate()

        args = self._get_command_args("extractorfs")
        mmseqs_output = run_mmseqs_command(args)

        self._handle_command_output(
            mmseqs_output=mmseqs_output,
            output_identifier="Extract ORFs",
            output_path=str(self.orf_db)
        )
