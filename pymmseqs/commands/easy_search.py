# pymmseqs/commands/easy_search.py

from pathlib import Path
from typing import Union, List

from ..config import EasySearchConfig
from ..parsers import EasySearchParser
from ..utils import tmp_dir_handler

def easy_search(
    query_fasta: Union[str, Path, List[Union[str, Path]]],
    target_fasta_or_db: Union[str, Path],
    alignment_file: Union[str, Path],

    # Optional parameters
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
    Required parameters
    ----------
    `query_fasta` : Union[str, Path]
        Path to one or more query FASTA files. Can be compressed with .gz or .bz2.

    `target_fasta_or_db` : Union[str, Path]
        Path to a target FASTA file (optionally compressed) or an MMseqs2 target database.

    `alignment_file` : Union[str, Path]
        Path to the output file where alignments will be stored.
    
    Optional parameters
    -------------------
    `tmp_dir` : Union[str, Path]
        Temporary directory for intermediate files.
        If not provided, a temporary directory will be created in the same directory as the alignment_file.

    `s` : float, optional
        Sensitivity
        - 1.0: faster
        - 4.0: fast
        - 5.7 (default)
        - 7.5: sensitive
    
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
    
    `max_seqs` : int, optional
        Maximum results per query passing prefilter
        - 300 (default)
        - Higher values increase sensitivity but may slow down the search
    
    `translate` : bool, optional
        Translate nucleotide input to amino acids before searching
        - False (default)
        - Set to True when queries are nucleotides and targets are proteins

    `translation_table` : int, optional
        Specifies the genetic code table to use 
            - 1: Canonical (default)
            - 2: Vert Mitochondrial
            - 3: Yeast Mitochondrial
            - 4: Mold Mitochondrial
            - 5: Invert Mitochondrial
            - 6: Ciliate
            - 9: Flatworm Mitochondrial
            - 10: Euplotid
            - 11: Prokaryote
            - 12: Alt Yeast
            - 13: Ascidian Mitochondrial
            - 14: Alt Flatworm Mitochondrial
            - 15: Blepharisma
            - 16: Chlorophycean Mitochondrial
            - 21: Trematode Mitochondrial
            - 22: Scenedesmus Mitochondrial
            - 23: Thraustochytrium Mitochondrial
            - 24: Pterobranchia Mitochondrial
            - 25: Gracilibacteria
            - 26: Pachysolen
            - 27: Karyorelict
            - 28: Condylostoma
            - 29: Mesodinium
            - 30: Pertrich
            - 31: Blastocrithidia
    
    `translation_mode` : int, optional
        How to translate nucleotide queries when `translate` is True
        - 0: ORFs (default)
        - 1: full reading frames

    `search_type` : int, optional
        Controls search mode
        - 0: auto-detect (default)
        - 1: amino acid vs amino acid
        - 2: translated nucleotide vs amino acid
        - 3: nucleotide vs nucleotide
        - 4: translated nucleotide alignment
    
    `format_output` : str, optional
        Comma-separated list of output columns to include in results.  
        Available columns:  
        - query, target, evalue, gapopen, pident, fident, nident, qstart, qend, qlen  
        - tstart, tend, tlen, alnlen, raw, bits, cigar, qseq, tseq, qheader, theader, qaln, taln  
        - qframe, tframe, mismatch, qcov, tcov, qset, qsetid, tset, tsetid, taxid, taxname, taxlineage  
        - qorfstart, qorfend, torfstart, torfend, ppos  

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

    config = EasySearchConfig(
        query_fasta=query_fasta,
        target_fasta_or_db=target_fasta_or_db,
        alignment_file=alignment_file,
        tmp_dir=tmp_dir,
        s=s,
        e=e,
        min_seq_id=min_seq_id,
        c=c,
        max_seqs=max_seqs,
        format_mode=4,
        translate=translate,
        translation_table=translation_table,
        translation_mode=translation_mode,
        search_type=search_type,
        format_output=format_output
    )

    config.run()

    return EasySearchParser(config)
