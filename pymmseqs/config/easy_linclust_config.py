from pathlib import Path
from typing import Union, List

from pymmseqs.config.base import BaseConfig
from pymmseqs.defaults import loader
from pymmseqs.utils import (
    get_caller_dir,
    run_mmseqs_command
)

DEFAULTS = loader.load("easy_linclust")

class EasyLinClustConfig(BaseConfig):
    def __init__(
        self,
        # Required parameters
        fasta_files: Union[str, Path],
        cluster_prefix: Union[str, Path],
        tmp_dir: Union[str, Path],

        # Prefilter parameters
        comp_bias_corr: bool = True,
        comp_bias_corr_scale: float = 1.0,
        add_self_matches: bool = False,
        alph_size: str = "aa:21,nucl:5",
        spaced_kmer_mode: int = 1,
        spaced_kmer_pattern: str = "",
        mask: bool = True,
        mask_prob: float = 0.9,
        mask_lower_case: bool = False,
        mask_n_repeat: int = 0,
        k: int = 0,
        split_memory_limit: str = "0",

        # Alignment parameters
        a: bool = False,
        alignment_mode: int = 0,
        alignment_output_mode: int = 0,
        wrapped_scoring: bool = False,
        e: float = 0.001,
        min_seq_id: float = 0.0,
        min_aln_len: int = 0,
        seq_id_mode: int = 0,
        alt_ali: int = 0,
        c: float = 0.8,
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
        
        # Clustering parameters
        cluster_mode: int = 0,
        max_iterations: int = 1000,
        similarity_type: int = 2,

        # K-mer matcher parameters
        weights: str = "",
        cluster_weight_threshold: float = 0.9,
        kmer_per_seq: int = 21,
        kmer_per_seq_scale: str = "aa:0.0,nucl:0.2",
        adjust_kmer_len: bool = False,
        hash_shift: int = 67,
        include_only_extendable: bool = False,
        ignore_multi_kmer: bool = False,

        # Profile parameters
        pca: float = 0.0,
        pcb: float = 0.0,

        # Misc parameters
        rescore_mode: int = 0,
        dbtype: int = 0,
        shuffle: bool = True,
        createdb_mode: int = 1,
        id_offset: int = 0,

        # Common parameters
        threads: Union[str, int] = 'all',
        compressed: bool = False,
        v: int = 3,
        sub_mat: str = "aa:blosum62.out,nucl:nucleotide.out",
        max_seq_len: int = 65535,
        db_load_mode: int = 0,
        remove_tmp_files: bool = True,
        force_reuse: bool = False,
        mpi_runner: str = "",


        # Expert parameters
        filter_hits: bool = False,
        sort_results: int = 0,
        write_lookup: bool = False,
    ):
        """
        Perform sequence clustering using MMseqs2.

        Parameters
        ----------
        `fasta_files` : Union[str, Path]
            Path to MMseqs2 sequence database created with createdb.

        `cluster_prefix` : Union[str, Path]
            Output cluster database path prefix (will create multiple files with this prefix).

        `tmp_dir` : Union[str, Path]
            Temporary directory for intermediate files (will be created if not exists).

        Prefilter Parameters
        --------------------
        `comp_bias_corr` : bool, optional
            Correct for locally biased amino acid composition
            - True (default)
            - False

        `comp_bias_corr_scale` : float, optional
            Scale factor for composition bias correction
            - Range 0, 1
            - 1.0 (default)

        `add_self_matches` : bool, optional
            Add entries of queries with themselves for clustering.
            - Default: False

        `alph_size` : str, optional
            Alphabet sizes for amino acid (protein) and nucleotide sequences (range 2-21)
            - "aa:21,nucl:5" (default)
                - aa:21: 20 amino acids + X for unknown residues
                - nucl:5: 4 nucleotides + N for unknown bases
        
        `spaced_kmer_mode` : int, optional
            Spaced k-mer mode
            - 0: consecutive
            - 1: spaced (default)

        `spaced_kmer_pattern` : str, optional
            Custom pattern for spaced k-mers used during k-mer matching.
            - Define a pattern of 1s (match positions) and 0s (ignore positions)
            - Example: "1101011" means 5 match positions and 2 ignored positions
            - Increases sensitivity by focusing on conserved regions while allowing flexibility in less conserved areas.

        `mask` : bool, optional
            Use low complexity masking
            - True (default)
            - False

        `mask_prob` : float, optional
            Probability threshold for masking low-complexity regions in sequences
            - 0.9 (default)
            - Sequences with low-complexity regions above this threshold are masked during k-mer matching

        `mask_lower_case` : bool, optional
            Mask lowercase letters in k-mer search.
            - True
            - False (default)
        
        `mask_n_repeat` : int, optional
            Repeat letters that occure > threshold in a row
            - 0 (default)

        `k` : int, optional
            k-mer length.
            - 0: automatically set to optimum (default)

        `split_memory_limit` : str, optional
            Maximum memory allocated per split for processing
            - "0":  all available (default)
                - Use suffixes like K, M, or G (e.g., "4G" for 4 gigabytes)

        Alignment Parameters
        --------------------
        `a` : bool, optional
            Add backtrace string (convert to alignments with mmseqs convertalis module)
            - True
            - False (default)

        `alignment_mode` : int, optional
            Alignment detail level
            - 0: auto (default)
            - 1: score + end_po
            - 2: + start_pos + cov
            - 3: + seq.id
            - 4: only ungapped alignment
            - 5: score only (output) cluster format
        
        `alignment_output_mode` : int, optional
            Alignment output mode
            - 0: automatic (default)
            - 1: only score and end_pos
            - 2: also start_pos and cov
            - 3: also seq.id
            - 4: only ungapped alignment
            - 5: score only (output) cluster format
        
        `wrapped_scoring` : bool, optional
            Enable wrapped diagonal scoring for nucleotide sequences by doubling the query sequence
            - True
            - False (default)

        `e` : float, optional
            E-value threshold (range 0.0, inf)
            - 0.001 (default)

        `min_seq_id` : float, optional
            Minimum sequence identity (range 0.0, 1.0)
            - 0.0 (default)

        `min_aln_len` : int, optional
            Minimum alignment length (range 0, inf)
            - 0 (default)

        `seq_id_mode` : int, optional
            Defines how sequence identity calculation is based on
            - 0: Alignment length (default)
            - 1: Shorter sequence
            - 2: Longer sequence

        `alt_ali` : int, optional
            Number of alternative alignments to show
            - 0 (default)

        `c` : float, optional
            Coverage threshold for alignments
            - 0.8 (default)
            - Determines the minimum fraction of aligned residues required for a match, based on the selected cov_mode

        `cov_mode` : int, optional
            Defines how alignment coverage is calculated:
            - 0: query + target (default)
            - 1: target only
            - 2: query only
            - 3: Target length ≥ x% query length
            - 4: Query length ≥ x% target length
            - 5: Short seq length ≥ x% other seq length

        `max_rejected` : int, optional
            Maximum rejected alignments before alignment calculation for a query is stopped
            - 2147483647 (default)

        `max_accept` : int, optional
            Maximum accepted alignments before alignment calculation for a query is stopped
            - 2147483647 (default)

        `score_bias` : float, optional
            Score bias added to alignment scores (in bits)
            - 0.0: no bias (default)
            - Adjusts alignment scores to favor or penalize certain alignments

        `realign` : bool, optional
            Compute more conservative, shorter alignments (scores and E-values not changed)
            - True
            - False (default)

        `realign_score_bias` : float, optional
            Additional score bias applied during realignment to compute more conservative alignments
            - -0.2 (default)
            - A negative value encourages shorter, more precise alignments

        `realign_max_seqs` : int, optional
            Maximum number of results to return in realignment
            - 2147483647 (default)

        `ccorr_score_weight` : float, optional
            Weight of backtrace correlation score that is added to the alignment score
            - 0.0 (default)
            - Higher values increase the influence of the backtrace correlation on the final alignment score

        `gap_open` : str, optional
            Gap open costs for amino acid (protein) and nucleotide alignments
            - "aa:11,nucl:5" (default)
                - aa:x: Gap open cost for protein alignments
                - nucl:x: Gap open cost for nucleotide alignments
            - Higher values penalize gap openings more heavily, favoring fewer gaps in alignments

        `gap_extend` : str, optional
            Gap extension costs for amino acid (protein) and nucleotide alignments
            - "aa:1,nucl:2" (default)
                - aa:x: Cost for extending a gap in protein alignments
                - nucl:x: Cost for extending a gap in nucleotide alignments
            - Lower values allow longer gaps; higher values penalize gaps more heavily

        `zdrop` : int, optional
            Maximum score drop allowed before truncating the alignment (nucleotide alignments only)
            - 40 (default)
                - Terminates alignments early in low-quality regions to improve computational efficiency
        
        Clustering Parameters
        ----------------------
        `cluster_mode` : int, optional
            Clustering method.
            - 0: Set-Cover (greedy) (default)
            - 1: Connected component (BLASTclust)
            - 2, 3: Greedy by sequence length (CDHIT)

        `max_iterations` : int, optional
            Maximum depth of breadth first search for connected component clustering.
            - 1000 (default)

        `similarity_type` : int, optional
            Type of score used for clustering.
            - 1: Alignment score
            - 2: Sequence identity (default)

        K-mer Matcher Parameters
        -------------------------
        `weights` : str, optional
            File containing sequence weights for cluster prioritization.
            Each line in the file corresponds to a sequence in the input
            file, with a floating-point weight determining its influence on clustering (higher values = higher priority)
            - Default: ""
            - Example: A text file like this:\n1.0\n0.5\n0.2\n

        `cluster_weight_threshold` : float, optional
            Weight threshold for cluster prioritization.
            - 0.9 (default)

        `kmer_per_seq` : int, optional
            Number of k-mers per sequence.
            - 21 (default)

        `kmer_per_seq_scale` : str, optional
            Scale k-mer per sequence based on sequence length as (kmer-pers-seq val + scale * seq-len)
            - "aa:0.0,nucl:0.2" (default)

        `adjust_kmer_len` : bool, optional
            Adjust k-mer length based on specificity (only for nucleotide).
            - True
            - False (default)

        `hash_shift` : int, optional
            Shift k-mer hash initialization.
            - 67 (default)

        `include_only_extendable` : bool, optional
            Include only extendable k-mers.
            - True
            - False (default)

        `ignore_multi_kmer` : bool, optional
            Skip k-mers occurring multiple times (>=2)
            - True
            - False (default)

        Profile Parameters
        -----------------
        `pca` : float, optional
            Pseudo count admixture strength for profile construction
            - 0.0 (default)
            - Higher values increase the weight of pseudo counts, making the profile more conservative
            - Lower values reduce their influence, making the profile more specific to the input sequences

        `pcb` : float, optional
            Controls the threshold for pseudo-count admixture based on the effective number of sequences (Neff) (range 0.0, inf)
            - 0.0 (default)
            - Lower values apply pseudo-counts more aggressively

        Misc Parameters
        ---------------
        `rescore_mode` : int, optional
            Rescore diagonals with:
            - 0: Hamming distance (default)
            - 1: local alignment (score only)
            - 2: local alignment
            - 3: global alignment
            - 4: longest alignment fulfilling window quality criterion
        
        `dbtype` : int, optional  
            Database type
            - 0: Auto-detect (default)
            - 1: Amino acid
            - 2: Nucleotide

        `shuffle` : bool, optional
            Shuffle the input database before processing.
            - True (default)
            - False

        `createdb_mode` : int, optional  
            Database creation mode
            - 0: Copy data
            - 1: Soft link data and write a new index (only works with single-line FASTA/Q files) (default)
        
        `id_offset` : int, optional
            Numeric IDs in index file are offset by this value
            - 0 (default)

        Common Parameters
        -----------------
        `threads` : int, optional
            CPU threads
            - 14 (default)
        
        `compressed` : bool, optional
            Compress output
            - True
            - False (default)

        `v` : int, optional
            Output verbosity
            - 0: quiet
            - 1: +errors
            - 2: +warnings
            - 3: +info (default)

        `sub_mat` : str, optional
            Substitution matrix (type:path, type:path)
            type: "aa" or "nucl"
            path: matrix file path
            - "aa:blosum62.out,nucl:nucleotide.out"

            Note: find available matrices in the MMseqs2 data directory: (https://github.com/soedinglab/MMseqs2/tree/master/data)

        `max_seq_len` : int, optional
            Maximum sequence length
            - 65535 (default)

        `db_load_mode` : int, optional
            Database preloading method
            - 0: auto (default)
            - 1: fread
            - 2: mmap
            - 3: mmap+touch

        `remove_tmp_files` : bool, optional
            Delete temporary files
            - True (default)
            - False

        `force_reuse` : bool, optional
            Reuse tmp filse in tmp/latest folder ignoring parameters and version changes
            - True
            - False (default)
        
        `mpi_runner` : str, optional
            Use MPI on compute cluster with this MPI command (e.g., "mpirun -np 42")
            - "" (default)

        Expert Parameters
        -----------------
        `filter_hits` : bool, optional
            Filter hits by sequence ID and coverage
            - True
            - False (default)

        `sort_results` : int, optional
            Result sorting method
            - 0: No sorting (default)
            - 1: E-value (Alignment) or sequence ID (Hamming) 
        
        `write_lookup` : bool, optional
            Create a `.lookup` file mapping internal IDs to FASTA IDs
            - True (default)
            - False

        Examples
        --------
        Linear clustering:
        >>> config = EasyLinClustConfig(
            "input.fasta",
            "results",
            "tmp"
        )
        >>> config.run()
        """
        super().__init__()
        
        # Initialize required parameters
        self.fasta_files = fasta_files if isinstance(fasta_files, list) else [fasta_files]
        self.cluster_prefix = Path(cluster_prefix)
        self.tmp_dir = Path(tmp_dir)

        # Initialize prefilter parameters
        self.comp_bias_corr = comp_bias_corr
        self.comp_bias_corr_scale = comp_bias_corr_scale
        self.add_self_matches = add_self_matches
        self.alph_size = alph_size
        self.spaced_kmer_mode = spaced_kmer_mode
        self.spaced_kmer_pattern = spaced_kmer_pattern
        self.mask = mask
        self.mask_prob = mask_prob
        self.mask_lower_case = mask_lower_case
        self.mask_n_repeat = mask_n_repeat
        self.k = k
        self.split_memory_limit = split_memory_limit

        # Initialize alignment parameters
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

        # Initialize clustering parameters
        self.cluster_mode = cluster_mode
        self.max_iterations = max_iterations
        self.similarity_type = similarity_type

        # Initialize k-mer matcher parameters
        self.weights = weights
        self.cluster_weight_threshold = cluster_weight_threshold
        self.kmer_per_seq = kmer_per_seq
        self.kmer_per_seq_scale = kmer_per_seq_scale
        self.adjust_kmer_len = adjust_kmer_len
        self.hash_shift = hash_shift
        self.include_only_extendable = include_only_extendable
        self.ignore_multi_kmer = ignore_multi_kmer

        # Initialize profile parameters
        self.pca = pca
        self.pcb = pcb

        # Initialize misc parameters
        self.rescore_mode = rescore_mode
        self.dbtype = dbtype
        self.shuffle = shuffle
        self.createdb_mode = createdb_mode
        self.id_offset = id_offset

        # Initialize common parameters
        self.threads = threads
        self.compressed = compressed
        self.v = v
        self.sub_mat = sub_mat
        self.max_seq_len = max_seq_len
        self.db_load_mode = db_load_mode
        self.remove_tmp_files = remove_tmp_files
        self.force_reuse = force_reuse
        self.mpi_runner = mpi_runner

        # Initialize expert parameters
        self.filter_hits = filter_hits
        self.sort_results = sort_results
        self.write_lookup = write_lookup

        self._defaults = DEFAULTS
        self._path_params = [param for param, info in DEFAULTS.items() if info['type'] == 'path']
        self._caller_dir = get_caller_dir()

    def _validate(self) -> None:
        self._check_required_files()
        self._validate_choices()

        # Validate numerical ranges
        if not (0.0 <= self.comp_bias_corr_scale <= 1.0):
            raise ValueError("comp_bias_corr_scale must be between 0.0 and 1.0")
        if not (0.0 <= self.mask_prob <= 1.0):
            raise ValueError("mask_prob must be between 0.0 and 1.0")
        if not (0.0 <= self.e):
            raise ValueError("e-value threshold must be >= 0.0")
        if not (0.0 <= self.min_seq_id <= 1.0):
            raise ValueError("min_seq_id must be between 0.0 and 1.0")
        if not (0.0 <= self.c <= 1.0):
            raise ValueError("coverage threshold (c) must be between 0.0 and 1.0")
        if not (0 <= self.min_aln_len):
            raise ValueError("min_aln_len must be >= 0")
        if not (0 <= self.max_rejected):
            raise ValueError("max_rejected must be >= 0")
        if not (0 <= self.max_accept):
            raise ValueError("max_accept must be >= 0")
        if not (0.0 <= self.cluster_weight_threshold <= 1.0):
            raise ValueError("cluster_weight_threshold must be between 0.0 and 1.0")
        if not (0 <= self.kmer_per_seq):
            raise ValueError("kmer_per_seq must be >= 0")
        if not (0 <= self.hash_shift):
            raise ValueError("hash_shift must be >= 0")
        if not (0 <= self.max_iterations):
            raise ValueError("max_iterations must be >= 0")
        if not (0 <= self.id_offset):
            raise ValueError("id_offset must be >= 0")
        if not (0 <= self.threads):
            raise ValueError("threads must be >= 0")
        if not (0 <= self.max_seq_len):
            raise ValueError("max_seq_len must be >= 0")
        

    def run(self) -> None:
        self._resolve_all_path(self._caller_dir)

        self._validate()

        args = self._get_command_args("easy-linclust")
        mmseqs_output = run_mmseqs_command(args)

        self._handle_command_output(
            mmseqs_output=mmseqs_output,
            output_identifier="Easy Linclust",
            output_path=str(self.cluster_prefix)
        )
