# pymmseqs/commands/extractorfs.py

from pathlib import Path
from typing import Union

from ..config import ExtractOrfsConfig
from ..parsers import ExtractOrfsParser

def extractorfs(
    sequence_db: Union[str, Path],
    orf_db: Union[str, Path],

    # Optional parameters
    min_length: int = 30,
    max_length: int = 32734,
    orf_start_mode: int = 1,
    contig_start_mode: int = 2,
    contig_end_mode: int = 2,
    forward_frames: str = "1,2,3",
    reverse_frames: str = "1,2,3",
    translation_table: int = 1,
    translate: bool = False,
    threads: Union[str, int] = 'all',

) -> ExtractOrfsParser:
    """
    Six-frame extraction of open reading frames (ORFs) from a nucleotide database.

    This is the core gene-finding step for metagenomics/genomics pipelines. It reads
    a nucleotide sequence database (from `createdb`) and writes the extracted ORFs to
    a new sequence database. The returned parser exposes the ORF coordinates merged
    with their sequences as a DataFrame; chain `to_path()` into `convert2fasta` or
    `search` for downstream work.

    Required parameters
    -------------------
    `sequence_db` : Union[str, Path]
        Path to the input (nucleotide) sequence database created with createdb.

    `orf_db` : Union[str, Path]
        Output path for the ORF sequence database.

    Optional parameters
    -------------------
    `min_length` : int, optional
        Minimum codon number in ORFs
        - 30 (default)

    `max_length` : int, optional
        Maximum codon number in ORFs
        - 32734 (default)

    `orf_start_mode` : int, optional
        ORF fragment handling
        - 0: from start to stop
        - 1: from any to stop (default)
        - 2: from last encountered start to stop (no start in the middle)

    `contig_start_mode` : int, optional
        Contig start handling
        - 0: incomplete, 1: complete, 2: both (default)

    `contig_end_mode` : int, optional
        Contig end handling
        - 0: incomplete, 1: complete, 2: both (default)

    `forward_frames` : str, optional
        Comma-separated forward-strand frames to extract
        - "1,2,3" (default)

    `reverse_frames` : str, optional
        Comma-separated reverse-strand frames to extract
        - "1,2,3" (default)

    `translation_table` : int, optional
        Genetic code table to use
        - 1: Canonical (default), 11: Prokaryote, ... (see MMseqs2 docs)

    `translate` : bool, optional
        Translate ORFs to amino acids
        - False (default)

    `threads` : Union[str, int], optional
        Number of CPU cores to use
        - 'all' (default), converted to all available cores

    Returns
    -------
    ExtractOrfsParser object
        - An ExtractOrfsParser exposing to_pandas()/to_list()/to_gen()/to_path()/summary().
    """

    config = ExtractOrfsConfig(
        sequence_db=sequence_db,
        orf_db=orf_db,
        min_length=min_length,
        max_length=max_length,
        orf_start_mode=orf_start_mode,
        contig_start_mode=contig_start_mode,
        contig_end_mode=contig_end_mode,
        forward_frames=forward_frames,
        reverse_frames=reverse_frames,
        translation_table=translation_table,
        translate=translate,
        threads=threads,
    )

    config.run()

    return ExtractOrfsParser(config)
