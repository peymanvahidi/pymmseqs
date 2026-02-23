# pymmseqs/config/createindex_config.py

from pathlib import Path
from typing import Union

from pymmseqs.config.base import BaseConfig
from pymmseqs.defaults import loader
from pymmseqs.utils import (
    get_caller_dir,
    run_mmseqs_command
)

DEFAULTS = loader.load("createindex")

class CreateIndexConfig(BaseConfig):
    def __init__(
        self,
        # Required parameters
        sequence_db: Union[str, Path],
        tmp_dir: Union[str, Path],

        # Prefilter parameters
        seed_sub_mat: str = "aa:VTML80.out,nucl:nucleotide.out",
        k: int = 0,
        alph_size: str = "aa:21,nucl:5",
        comp_bias_corr: bool = True,
        comp_bias_corr_scale: float = 1.0,
        max_seqs: int = 300,
        mask: bool = True,
        mask_prob: float = 0.9,
        mask_lower_case: bool = False,
        mask_n_repeat: int = 0,
        spaced_kmer_mode: int = 1,
        spaced_kmer_pattern: str = "",
        s: float = 7.5,
        k_score: str = "seq:0,prof:0",
        split: int = 0,
        split_memory_limit: str = "0",

        # Misc parameters
        check_compatible: int = 0,
        search_type: int = 0,
        min_length: int = 30,
        max_length: int = 32734,
        max_gaps: int = 2147483647,
        contig_start_mode: int = 2,
        contig_end_mode: int = 2,
        orf_start_mode: int = 1,
        forward_frames: str = "1,2,3",
        reverse_frames: str = "1,2,3",
        translation_table: int = 1,
        translate: int = 0,
        use_all_table_starts: bool = False,
        id_offset: int = 0,
        sequence_overlap: int = 0,
        sequence_split_mode: int = 1,
        headers_split_mode: int = 0,
        translation_mode: int = 0,

        # Common parameters
        max_seq_len: int = 65535,
        v: int = 3,
        threads: Union[str, int] = 'all',
        compressed: bool = False,
        remove_tmp_files: bool = False,

        # Expert parameters
        index_subset: int = 0,
        create_lookup: int = 0,
        strand: int = 1,
    ):
        """
        Create an index for a sequence database using MMseqs2.

        Parameters
        ----------
        `sequence_db` : Union[str, Path]
            Path to MMseqs2 sequence database created with createdb.

        `tmp_dir` : Union[str, Path]
            Temporary directory for intermediate files (will be created if not exists).

        Prefilter Parameters
        --------------------
        `seed_sub_mat` : str, optional
            Substitution matrix for k-mer generation as (type:path, type:path)
            - Default: "aa:VTML80.out,nucl:nucleotide.out"

        `k` : int, optional
            k-mer length.
            - 0: automatically set to optimum (default)

        `alph_size` : str, optional
            Alphabet sizes for amino acid (protein) and nucleotide sequences (range 2-21)
            - Default: "aa:21,nucl:5"
                - aa:21: 20 amino acids + X for unknown residues
                - nucl:5: 4 nucleotides + N for unknown bases

        `comp_bias_corr` : bool, optional
            Correct for locally biased amino acid composition
            - True (default)
            - False

        `comp_bias_corr_scale` : float, optional
            Scale factor for composition bias correction
            - Range 0, 1
            - 1.0 (default)

        `max_seqs` : int, optional
            Maximum results per query sequence passing prefilter.
            - Default: 300

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
            Repeat letters that occur > threshold in a row
            - 0 (default)

        `spaced_kmer_mode` : int, optional
            Spaced k-mer mode
            - 0: consecutive
            - 1: spaced (default)

        `spaced_kmer_pattern` : str, optional
            Custom pattern for spaced k-mers used during k-mer matching.
            - Define a pattern of 1s (match positions) and 0s (ignore positions)
            - Example: "1101011" means 5 match positions and 2 ignored positions
            - Increases sensitivity by focusing on conserved regions while allowing flexibility in less conserved areas.

        `s` : float, optional
            Sensitivity.
            - Options: 1.0 (faster), 4.0 (fast), 7.5 (sensitive)
            - Default: 7.5

        `k_score` : str, optional
            k-mer thresholds for sequence and profile searches.
            - Default: "seq:0,prof:0"

        `split` : int, optional
            Split input into N chunks
            - 0: set the best split automatically (default)

        `split_memory_limit` : str, optional
            Maximum memory allocated per split for processing
            - "0": all available (default)
            - Use suffixes like K, M, or G (e.g., "4G" for 4 gigabytes)

        Misc Parameters
        ---------------
        `check_compatible` : int, optional
            Index compatibility check
            - 0: always recreate (default)
            - 1: check if recreating is needed
            - 2: fail if index is incompatible

        `search_type` : int, optional
            Search type
            - 0: auto (default)
            - 1: amino acid
            - 2: translated
            - 3: nucleotide
            - 4: translated nucleotide alignment

        `min_length` : int, optional
            Minimum codon number in open reading frames
            - 30 (default)

        `max_length` : int, optional
            Maximum codon number in open reading frames
            - 32734 (default)

        `max_gaps` : int, optional
            Maximum number of codons with gaps or unknown residues before an open reading frame is rejected
            - 2147483647 (default)

        `contig_start_mode` : int, optional
            Contig start mode
            - 0: incomplete
            - 1: complete
            - 2: both (default)

        `contig_end_mode` : int, optional
            Contig end mode
            - 0: incomplete
            - 1: complete
            - 2: both (default)

        `orf_start_mode` : int, optional
            ORF fragment mode
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
            Translation table for genetic code
            - 1: CANONICAL (default)
            - Other options: 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31

        `translate` : int, optional
            Translate ORF to amino acid
            - 0 (default)

        `use_all_table_starts` : bool, optional
            Use all alternatives for a start codon in the genetic table
            - True
            - False: only ATG (AUG) (default)

        `id_offset` : int, optional
            Numeric IDs in index file are offset by this value
            - 0 (default)

        `sequence_overlap` : int, optional
            Overlap between sequences
            - 0 (default)

        `sequence_split_mode` : int, optional
            Sequence split mode
            - 0: copy data
            - 1: soft link data and write new index (default)

        `headers_split_mode` : int, optional
            Header split mode
            - 0: split position (default)
            - 1: original header

        `translation_mode` : int, optional
            Translation AA seq from nucleotide by
            - 0: ORFs (default)
            - 1: full reading frames

        Common Parameters
        -----------------
        `max_seq_len` : int, optional
            Maximum sequence length
            - 65535 (default)

        `v` : int, optional
            Output verbosity
            - 0: quiet
            - 1: +errors
            - 2: +warnings
            - 3: +info (default)

        `threads` : int, optional
            CPU threads
            - 14 (default)

        `compressed` : bool, optional
            Compress output
            - True
            - False (default)

        `remove_tmp_files` : bool, optional
            Delete temporary files
            - True
            - False (default)

        Expert Parameters
        -----------------
        `index_subset` : int, optional
            Create specialized index with subset of entries
            - 0: normal index (default)
            - 1: index without headers
            - 2: index without prefiltering data
            - 4: index without aln (for cluster db)
            - Flags can be combined bit wise

        `create_lookup` : int, optional
            Create database lookup file (can be very large)
            - 0 (default)

        `strand` : int, optional
            Strand selection for DNA/DNA search
            - 0: reverse
            - 1: forward (default)
            - 2: both

        Returns
        -------
        None
            Creates the index for the specified sequence database.

        Raises
        ------
        FileNotFoundError
            If input database is missing.

        ValueError
            For invalid parameter combinations or values.
        """
        super().__init__()

        # Required parameters
        self.sequence_db = Path(sequence_db)
        self.tmp_dir = Path(tmp_dir)

        # Prefilter parameters
        self.seed_sub_mat = seed_sub_mat
        self.k = k
        self.alph_size = alph_size
        self.comp_bias_corr = comp_bias_corr
        self.comp_bias_corr_scale = comp_bias_corr_scale
        self.max_seqs = max_seqs
        self.mask = mask
        self.mask_prob = mask_prob
        self.mask_lower_case = mask_lower_case
        self.mask_n_repeat = mask_n_repeat
        self.spaced_kmer_mode = spaced_kmer_mode
        self.spaced_kmer_pattern = spaced_kmer_pattern
        self.s = s
        self.k_score = k_score
        self.split = split
        self.split_memory_limit = split_memory_limit

        # Misc parameters
        self.check_compatible = check_compatible
        self.search_type = search_type
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
        self.sequence_overlap = sequence_overlap
        self.sequence_split_mode = sequence_split_mode
        self.headers_split_mode = headers_split_mode
        self.translation_mode = translation_mode

        # Common parameters
        self.max_seq_len = max_seq_len
        self.v = v
        self.threads = threads
        self.compressed = compressed
        self.remove_tmp_files = remove_tmp_files

        # Expert parameters
        self.index_subset = index_subset
        self.create_lookup = create_lookup
        self.strand = strand

        self._defaults = DEFAULTS
        self._path_params = [param for param, info in DEFAULTS.items() if info['type'] == 'path']
        self._caller_dir = get_caller_dir()

    def _validate(self) -> None:
        self._check_required_files()
        self._validate_choices()

        # Validate numerical ranges
        if not (0 <= self.comp_bias_corr_scale <= 1):
            raise ValueError("comp_bias_corr_scale must be between 0 and 1")
        if not (0.0 <= self.mask_prob <= 1.0):
            raise ValueError("mask_prob must be between 0.0 and 1.0")

    def run(self) -> None:
        self._resolve_all_path(self._caller_dir)

        self._validate()

        args = self._get_command_args("createindex")
        mmseqs_output = run_mmseqs_command(args)

        self._handle_command_output(
            mmseqs_output=mmseqs_output,
            output_identifier="Create Index",
            output_path=str(self.sequence_db)
        )
