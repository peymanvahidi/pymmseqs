# pymmseqs/config/align_config.py

from pathlib import Path
from typing import Union

from .base import BaseConfig
from ..defaults import loader
from ..utils import (
    get_caller_dir,
    run_mmseqs_command
)

DEFAULTS = loader.load("align")

class AlignConfig(BaseConfig):
    def __init__(
        self,
        # Required parameters
        query_db: Union[str, Path],
        target_db: Union[str, Path],
        result_db: Union[str, Path],
        alignment_db: Union[str, Path],
        
        # Alignment parameters
        comp_bias_corr: int = 1,
        comp_bias_corr_scale: float = 1.0,
        add_self_matches: bool = False,
        a: bool = False,
        alignment_mode: int = 0,
        alignment_output_mode: int = 0,
        wrapped_scoring: bool = False,
        e: float = 0.001,
        min_seq_id: float = 0.0,
        min_aln_len: int = 0,
        seq_id_mode: int = 0,
        alt_ali: int = 0,
        c: float = 0.0,
        cov_mode: int = 0,
        max_rejected: int = 2147483647,
        max_accept: int = 2147483647,
        score_bias: float = 0.0,
        realign: bool = False,
        realign_score_bias: float = -0.2,
        realign_max_seqs: int = 2147483647,
        corr_score_weight: float = 0.0,
        gap_open: str = "aa:11,nucl:5",
        gap_extend: str = "aa:1,nucl:2",
        zdrop: int = 40,
        
        # Profile parameters
        pca: float = None,
        pcb: float = None,
        
        # Common parameters
        sub_mat: str = "aa:blosum62.out,nucl:nucleotide.out",
        max_seq_len: int = 65535,
        db_load_mode: int = 0,
        threads: Union[str, int] = 'all',
        compressed: int = 0,
        v: int = 3,
    ):
        """
        Create alignments between sequences using MMseqs2 align module.

        Parameters
        ----------
        `query_db` : Union[str, Path]
            The input query database

        `target_db` : Union[str, Path]
            The input target database

        `result_db` : Union[str, Path]
            The input result database where search results are stored

        `alignment_db` : Union[str, Path]
            Database outputwhere alignments will be stored

        Alignment Parameters
        -------------------
        `comp_bias_corr` : int, optional
            Correct for locally biased amino acid composition
            - 0: disabled
            - 1: enabled (default)

        `comp_bias_corr_scale` : float, optional
            Scale for correcting locally biased amino acid composition
            - Range 0.0-1.0
            - 1.0 (default)

        `add_self_matches` : bool, optional
            Artificially add entries of queries with themselves for clustering
            - True
            - False (default)

        `a` : bool, optional
            Add backtrace string (convert to alignments with mmseqs convertalis module)
            - True
            - False (default)

        `alignment_mode` : int, optional
            How to compute the alignment
            - 0: automatic (default)
            - 1: only score and end_pos
            - 2: also start_pos and cov
            - 3: also seq.id

        `alignment_output_mode` : int, optional
            Output mode for the alignment
            - 0: automatic (default)
            - 1: only score + end_pos
            - 2: also start_pos + cov
            - 3: also seq.id
            - 4: only ungapped alignment
            - 5: score only (output) cluster format

        `wrapped_scoring` : bool, optional
            Double the (nucleotide) query sequence during the scoring process to allow wrapped 
            diagonal scoring around end and start
            - True
            - False (default)

        `e` : float, optional
            List matches below this E-value
            - Range 0.0-inf
            - 0.001 (default)

        `min_seq_id` : float, optional
            List matches above this sequence identity for clustering
            - Range 0.0-1.0
            - 0.0 (default)

        `min_aln_len` : int, optional
            Minimum alignment length
            - Range 0-INT_MAX
            - 0 (default)

        `seq_id_mode` : int, optional
            Sequence identity mode for alignment length
            - 0: alignment length (default)
            - 1: shorter sequence
            - 2: longer sequence

        `alt_ali` : int, optional
            Show up to this many alternative alignments
            - 0 (default)

        `c` : float, optional
            List matches above this fraction of aligned (covered) residues
            - 0.0 (default)

        `cov_mode` : int, optional
            Coverage mode for alignment
            - 0: coverage of query and target (default)
            - 1: coverage of target
            - 2: coverage of query
            - 3: target seq. length has to be at least x% of query length
            - 4: query seq. length has to be at least x% of target length
            - 5: short seq. needs to be at least x% of the other seq. length

        `max_rejected` : int, optional
            Maximum rejected alignments before alignment calculation for a query is stopped
            - 2147483647 (default)

        `max_accept` : int, optional
            Maximum accepted alignments before alignment calculation for a query is stopped
            - 2147483647 (default)

        `score_bias` : float, optional
            Score bias when computing SW alignment (in bits)
            - 0.0 (default)

        `realign` : bool, optional
            Compute more conservative, shorter alignments (scores and E-values not changed)
            - True
            - False (default)

        `realign_score_bias` : float, optional
            Additional bias when computing realignment
            - -0.2 (default)

        `realign_max_seqs` : int, optional
            Maximum number of results to return in realignment
            - 2147483647 (default)

        `corr_score_weight` : float, optional
            Weight of backtrace correlation score that is added to the alignment score
            - 0.0 (default)

        `gap_open` : str, optional
            Gap open cost, different settings for amino acids and nucleotides
            - "aa:11,nucl:5" (default)

        `gap_extend` : str, optional
            Gap extension cost, different settings for amino acids and nucleotides
            - "aa:1,nucl:2" (default)

        `zdrop` : int, optional
            Maximal allowed difference between score values before alignment is truncated (nucleotide alignment only)
            - 40 (default)

        Profile Parameters
        -----------------
        `pca` : float, optional
            Pseudo count admixture strength
            - None (default)

        `pcb` : float, optional
            Pseudo counts: Neff at half of maximum admixture (range 0.0-inf)
            - None (default)

        Common Parameters
        ----------------
        `sub_mat` : str, optional
            Substitution matrix file
            - "aa:blosum62.out,nucl:nucleotide.out" (default)

        `max_seq_len` : int, optional
            Maximum sequence length
            - 65535 (default)

        `db_load_mode` : int, optional
            Database preload mode
            - 0: auto (default)
            - 1: fread
            - 2: mmap
            - 3: mmap+touch

        `threads` : int, optional
            Number of CPU-cores used
            - 14 (default)

        `compressed` : int, optional
            Write compressed output
            - 0: disabled (default)
            - 1: enabled

        `v` : int, optional
            Verbosity level
            - 0: quiet
            - 1: +errors
            - 2: +warnings
            - 3: +info (default)

        Examples
        --------
        >>> config = AlignConfig(
            query_db="query_db",
            target_db="target_db",
            result_db="result_db",
            alignment_db="alignment_db",
            threads=8
        )
        >>> config.run()
        """
        super().__init__()
        
        # Required parameters
        self.query_db = Path(query_db)
        self.target_db = Path(target_db)
        self.result_db = Path(result_db)
        self.alignment_db = Path(alignment_db)
        
        # Alignment parameters
        self.comp_bias_corr = comp_bias_corr
        self.comp_bias_corr_scale = comp_bias_corr_scale
        self.add_self_matches = add_self_matches
        self.a = a
        self.alignment_mode = alignment_mode
        self.alignment_output_mode = alignment_output_mode
        self.wrapped_scoring = wrapped_scoring
        self.e = e
        self.min_seq_id = min_seq_id
        self.min_aln_len = min_aln_len
        self.seq_id_mode = seq_id_mode
        self.alt_ali = alt_ali
        self.c = c
        self.cov_mode = cov_mode
        self.max_rejected = max_rejected
        self.max_accept = max_accept
        self.score_bias = score_bias
        self.realign = realign
        self.realign_score_bias = realign_score_bias
        self.realign_max_seqs = realign_max_seqs
        self.corr_score_weight = corr_score_weight
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        self.zdrop = zdrop
        
        # Profile parameters
        self.pca = pca
        self.pcb = pcb
        
        # Common parameters
        self.sub_mat = sub_mat
        self.max_seq_len = max_seq_len
        self.db_load_mode = db_load_mode
        self.threads = threads
        self.compressed = compressed
        self.v = v
        
        self._defaults = DEFAULTS
        self._path_params = [param for param, info in DEFAULTS.items() if info['type'] == 'path']
        self._caller_dir = get_caller_dir()
    
    def _validate(self) -> None:
        self._check_required_files()
        self._validate_choices()
        
        # Additional validations
        if not (0.0 <= self.comp_bias_corr_scale <= 1.0):
            raise ValueError("comp_bias_corr_scale must be between 0.0 and 1.0")
        if not (0.0 <= self.min_seq_id <= 1.0):
            raise ValueError("min_seq_id must be between 0.0 and 1.0")
        if not (0.0 <= self.c <= 1.0):
            raise ValueError("c must be between 0.0 and 1.0")
        if not (self.threads >= 1):
            raise ValueError("threads must be at least 1")
        if not (self.min_aln_len >= 0):
            raise ValueError("min_aln_len cannot be negative")
    
    def run(self) -> None:
        self._resolve_all_path(self._caller_dir)
        
        self._validate()
        
        args = self._get_command_args("align")
        mmseqs_output = run_mmseqs_command(args)
        
        self._handle_command_output(
            mmseqs_output=mmseqs_output,
            output_identifier="Align",
            output_path=str(self.alignment_db)
        )
