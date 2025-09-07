from pathlib import Path
from typing import Union, List

from ..config import (
    CreateDBConfig,
    CreateIndexConfig,
    TouchDBConfig,
    SearchConfig,
    ConvertAlisConfig,
)
from ..parsers import EasySearchParser
from ..utils import tmp_dir_handler


def fast_easy_search(
    query_fasta: Union[str, Path, List[Union[str, Path]]],
    target_fasta: Union[str, Path],
    alignment_file: Union[str, Path],

    # Optional parameters (mirrored from easy_search where applicable)
    tmp_dir: Union[str, Path, None] = None,
    s: float = 5.7,
    e: float = 0.001,
    min_seq_id: float = 0.0,
    c: float = 0.0,
    max_seqs: int = 300,
    translate: bool = False,
    translation_table: int = 1,
    translation_mode: int = 0,
    search_type: int = 0,
    format_output: str = "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits",
) -> EasySearchParser:
    """
    Perform a fast easy search by creating DBs, indexing, touching the index into memory,
    running search, and converting alignments.

    Parameters are aligned with `easy_search` to keep the same user-facing options.
    """

    tmp_dir = tmp_dir_handler(
        tmp_dir=tmp_dir,
        output_file_path=alignment_file
    )

    # Paths derived from alignment_file for intermediate DBs
    base_dir = Path(tmp_dir)
    query_db = base_dir / "query_db"
    target_db = base_dir / "target_db"
    alignment_db = base_dir / "result_db"

    # 1) createdb for query and target
    CreateDBConfig(
        fasta_file=query_fasta,
        sequence_db=query_db
    ).run()

    CreateDBConfig(
        fasta_file=target_fasta,
        sequence_db=target_db
    ).run()

    # 2) createindex for target
    CreateIndexConfig(
        sequence_db=target_db,
        tmp_dir=tmp_dir
    ).run()

    # 3) touchdb (load index into memory)
    TouchDBConfig(
        sequence_db=target_db,
    ).run()

    # 4) search with db_load_mode=2 (mmap)
    SearchConfig(
        query_db=query_db,
        target_db=target_db,
        alignment_db=alignment_db,
        tmp_dir=tmp_dir,
        s=s,
        e=e,
        min_seq_id=min_seq_id,
        c=c,
        max_seqs=max_seqs,
        translate=translate,
        translation_table=translation_table,
        translation_mode=translation_mode,
        search_type=search_type,
        db_load_mode=2,
    ).run()

    # 5) convertalis to human-readable file with headers
    convert_config = ConvertAlisConfig(
        query_db=query_db,
        target_db=target_db,
        alignment_db=alignment_db,
        alignment_file=alignment_file,
        format_mode=4,
        format_output=format_output,
        translation_table=translation_table,
        search_type=search_type,
    )
    convert_config.run()

    return EasySearchParser(convert_config)


