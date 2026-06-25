# pymmseqs/commands/easy_linsearch.py

from pathlib import Path
from typing import Union

from ..config import EasyLinSearchConfig
from ..parsers import EasySearchParser
from ..utils import tmp_dir_handler

def easy_linsearch(
    query_fasta: Union[str, Path],
    target_fasta_or_db: Union[str, Path],
    alignment_file: Union[str, Path],

    # Optional parameters
    tmp_dir: Union[str, Path, None] = None,
    e: float = 0.001,
    min_seq_id: float = 0.0,
    c: float = 0.0,
    translate: bool = False,
    translation_table: int = 1,
    search_type: int = 0,
    format_output: str = "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits",

) -> EasySearchParser:
    """
    Run a fast, linear-time homology search (MMseqs2 easy-linsearch).

    This is the fast, less sensitive sibling of `easy_search`. It searches a query
    FASTA against a target FASTA/DB and writes a BLAST-tab alignment table that is
    parsed back into Python.

    Required parameters
    -------------------
    `query_fasta` : Union[str, Path]
        Path to a query FASTA file. Can be compressed with .gz or .bz2.

    `target_fasta_or_db` : Union[str, Path]
        Path to a target FASTA file (optionally compressed) or an MMseqs2 target database.

    `alignment_file` : Union[str, Path]
        Path to the output file where alignments will be stored.

    Optional parameters
    -------------------
    `tmp_dir` : Union[str, Path]
        Temporary directory for intermediate files.
        If not provided, a temporary directory will be created in the same directory as the alignment_file.

    `e` : float, optional
        E-value threshold (range 0.0, inf)
        - 0.001 (default)

    `min_seq_id` : float, optional
        Minimum sequence identity (range 0.0, 1.0)
        - 0.0 (default)

    `c` : float, optional
        Coverage threshold for alignments
        - 0.0 (default)
        - Determines the minimum fraction of aligned residues required for a match, based on the selected cov_mode

    `translate` : bool, optional
        Translate ORFs to amino acids before searching
        - False (default)

    `translation_table` : int, optional
        Genetic code table to use when translating
        - 1: Canonical (default)
        - See EasyLinSearchConfig for the full list of supported tables.

    `search_type` : int, optional
        Controls search mode
        - 0: auto-detect (default)
        - 1: amino acid vs amino acid
        - 2: translated
        - 3: nucleotide vs nucleotide
        - 4: translated alignment

    `format_output` : str, optional
        Comma-separated list of output columns to include in results.
        See EasyLinSearchConfig for the full list of available columns.
        - Default: "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"

    Returns
    -------
    EasySearchParser object
        - An EasySearchParser instance that provides methods to access and parse the alignment data.
    """

    tmp_dir = tmp_dir_handler(
        tmp_dir=tmp_dir,
        output_file_path=alignment_file
    )

    config = EasyLinSearchConfig(
        query_fasta=query_fasta,
        target_fasta_or_db=target_fasta_or_db,
        alignment_file=alignment_file,
        tmp_dir=tmp_dir,
        e=e,
        min_seq_id=min_seq_id,
        c=c,
        format_mode=4,
        translate=translate,
        translation_table=translation_table,
        search_type=search_type,
        format_output=format_output
    )

    config.run()

    return EasySearchParser(config)
