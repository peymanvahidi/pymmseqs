# pymmseqs/commands/convertalis.py

from pathlib import Path
from typing import Union

from ..config import ConvertAlisConfig
from ..parsers import EasySearchParser

def convertalis(
    query_db: Union[str, Path],
    target_db: Union[str, Path],
    alignment_db: Union[str, Path],
    alignment_file: Union[str, Path],

    # Optional parameters
    format_output: str = "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits",
    search_type: int = 0,
    translation_table: int = 1,
    threads: Union[str, int] = 'all',

) -> EasySearchParser:
    """
    Convert an MMseqs2 alignment database into a parsable BLAST-tab table.

    This is the bridge from a database-level `search()` result (an alignment DB,
    produced when `search` is run with `a=True`) to a readable, Python-parsable
    alignment table. The output is written with column headers (format_mode=4) and
    returned via an `EasySearchParser` (to_pandas/to_list/to_gen/to_path).

    For SAM (format_mode=1) or pretty HTML (format_mode=3) output, use
    `ConvertAlisConfig` directly and read the file via `to_path()`.

    Required parameters
    -------------------
    `query_db` : Union[str, Path]
        Path to the query sequence database created with createdb.

    `target_db` : Union[str, Path]
        Path to the target sequence database created with createdb.

    `alignment_db` : Union[str, Path]
        Path to the alignment database to convert (produced by search/another module).

    `alignment_file` : Union[str, Path]
        Output path for the human-readable alignment table.

    Optional parameters
    -------------------
    `format_output` : str, optional
        Comma-separated list of output columns to include in results.
        See ConvertAlisConfig for the full list of available columns.
        - Default: "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"

    `search_type` : int, optional
        Search mode
        - 0: auto (default)
        - 1: amino
        - 2: translated
        - 3: nucleotide
        - 4: translated alignment

    `translation_table` : int, optional
        Genetic code table to use when applicable
        - 1: Canonical (default)
        - See ConvertAlisConfig for the full list of supported tables.

    `threads` : Union[str, int], optional
        Number of CPU cores to use
        - 'all' (default), converted to all available cores

    Returns
    -------
    EasySearchParser object
        - An EasySearchParser instance that provides methods to access and parse the alignment data.
    """

    config = ConvertAlisConfig(
        query_db=query_db,
        target_db=target_db,
        alignment_db=alignment_db,
        alignment_file=alignment_file,
        format_mode=4,
        format_output=format_output,
        search_type=search_type,
        translation_table=translation_table,
        threads=threads,
    )

    config.run()

    return EasySearchParser(config)
