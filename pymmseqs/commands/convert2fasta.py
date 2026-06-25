# pymmseqs/commands/convert2fasta.py

from pathlib import Path
from typing import Union

from ..config import Convert2FastaConfig
from ..parsers import Convert2FastaParser

def convert2fasta(
    sequence_db: Union[str, Path],
    fasta_file: Union[str, Path],

    # Optional parameters
    use_header_file: bool = False,
    v: int = 3,

) -> Convert2FastaParser:
    """
    Convert an MMseqs2 sequence database back to a FASTA file.

    This is the reverse of `createdb`: it dumps a (possibly filtered or
    transformed) sequence database — e.g. the output of `extractorfs` — to a
    plain FASTA file that can be read back into Python.

    Required parameters
    -------------------
    `sequence_db` : Union[str, Path]
        Path to the input sequence database (from createdb, extractorfs, etc.).

    `fasta_file` : Union[str, Path]
        Output path for the FASTA file.

    Optional parameters
    -------------------
    `use_header_file` : bool, optional
        Use the sequence header DB instead of the body to map the entry keys
        - False (default)
        - True

    `v` : int, optional
        Verbosity level of the output
        - 3 (default)

    Returns
    -------
    Convert2FastaParser object
        - A Convert2FastaParser instance exposing to_path(), to_gen(), to_list(), and to_pandas().
    """

    config = Convert2FastaConfig(
        sequence_db=sequence_db,
        fasta_file=fasta_file,
        use_header_file=use_header_file,
        v=v,
    )

    config.run()

    return Convert2FastaParser(config)
