# pymmseqs/config/easy_taxonomy_config.py

from pathlib import Path
from typing import Union, List, Optional

from .base import BaseConfig
from ..defaults import loader
from ..utils import (
    get_caller_dir,
    run_mmseqs_command
)

DEFAULTS = loader.load("easy_taxonomy")


class EasyTaxonomyConfig(BaseConfig):
    """Assign taxonomic labels to FASTA sequences using MMseqs2.

    Parameters
    ----------
    `fasta_file` : Union[List[Union[str, Path]], Union[str, Path]]
        Path(s) to the input FASTA file(s), optionally compressed as .gz or .bz2.

    `target_db` : Union[str, Path]
        Target sequence database with taxonomy information (e.g., UniProtKB/Swiss-Prot).

    `tax_reports` : Union[str, Path]
        Output prefix for taxonomy report files. Creates:
        - {tax_reports}_tophit_aln: top hits
        - {tax_reports}_tophit_report: coverage profiles per database entry
        - {tax_reports}_report: kraken-style report

    `tmp_dir` : Union[str, Path]
        Temporary directory for intermediate files.

    `lca_mode` : int, optional
        LCA Mode 1: single search LCA, 2/3: approximate 2bLCA, 4: top hit.
        - Default: 3

    `s` : float, optional
        Sensitivity: 1.0 faster; 4.0 fast; 7.5 sensitive.
        - Default: 4.0

    `e` : float, optional
        E-value threshold (range 0.0-inf).
        - Default: 0.001

    `threads` : Union[str, int], optional
        Number of CPU-cores used.
        - Default: 'all'
    """

    def __init__(
        self,
        # Required parameters
        fasta_file: Union[List[Union[str, Path]], Union[str, Path]],
        target_db: Union[str, Path],
        tax_reports: Union[str, Path],
        tmp_dir: Union[str, Path],

        # Prefilter parameters
        comp_bias_corr: bool = True,
        comp_bias_corr_scale: float = 1.0,
        add_self_matches: bool = False,
        seed_sub_mat: str = "aa:VTML80.out,nucl:nucleotide.out",
        s: float = 4.0,
        k: int = 0,
        target_search_mode: int = 0,
        k_score: str = "seq:2147483647,prof:2147483647",
        alph_size: str = "aa:21,nucl:5",
        max_seqs: int = 300,
        split: int = 0,
        split_mode: int = 2,
        split_memory_limit: str = "0",
        diag_score: bool = True,
        exact_kmer_matching: bool = False,
        mask: bool = True,
        mask_prob: float = 0.9,
        mask_lower_case: bool = False,
        mask_n_repeat: int = 0,
        min_ungapped_score: int = 15,
        spaced_kmer_mode: int = 1,
        spaced_kmer_pattern: str = "",
        local_tmp: Union[str, Path] = "",
        disk_space_limit: str = "0",

        # Align parameters
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
        exhaustive_search_filter: int = 0,

        # Profile parameters
        pca: Optional[float] = None,
        pcb: Optional[float] = None,
        mask_profile: int = 1,
        e_profile: float = 0.001,
        wg: bool = False,
        filter_msa: int = 1,
        filter_min_enable: int = 0,
        max_seq_id: float = 0.9,
        qid: str = "0.0",
        qsc: float = -20.0,
        cov: float = 0.0,
        diff: int = 1000,
        pseudo_cnt_mode: int = 0,
        profile_output_mode: int = 0,
        exhaustive_search: bool = False,
        lca_search: bool = False,

        # Misc parameters
        orf_filter: int = 0,
        orf_filter_e: float = 100.0,
        orf_filter_s: float = 2.0,
        lca_mode: int = 3,
        majority: float = 0.5,
        vote_mode: int = 1,
        lca_ranks: str = "",
        tax_lineage: int = 0,
        blacklist: str = "12908:unclassified sequences,28384:other sequences",
        taxon_list: str = "",
        prefilter_mode: int = 0,
        rescore_mode: int = 0,
        allow_deletion: bool = False,
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
        search_type: int = 0,
        translation_mode: int = 0,
        report_mode: int = 0,
        format_mode: int = 0,
        format_output: str = "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits",
        first_seq_as_repr: bool = False,
        target_column: int = 1,
        full_header: bool = False,
        idx_seq_src: int = 0,
        dbtype: int = 0,
        shuffle: bool = True,
        createdb_mode: int = 1,

        # Common parameters
        compressed: int = 0,
        threads: Union[str, int] = 'all',
        v: int = 3,
        sub_mat: str = "aa:blosum62.out,nucl:nucleotide.out",
        max_seq_len: int = 65535,
        db_load_mode: int = 0,
        gpu: int = 0,
        gpu_server: int = 0,
        gpu_server_wait_timeout: int = 600,
        mpi_runner: str = "",
        force_reuse: bool = False,
        remove_tmp_files: bool = True,

        # Expert parameters
        filter_hits: bool = False,
        sort_results: int = 0,
        create_lookup: int = 0,
        chain_alignments: int = 0,
        merge_query: int = 1,
        strand: int = 1,
        db_output: bool = False,
        write_lookup: int = 0,
    ):
        super().__init__()

        # Required parameters
        self.fasta_file = fasta_file if isinstance(fasta_file, list) else [fasta_file]
        self.fasta_file = [Path(f) for f in self.fasta_file]
        self.target_db = Path(target_db)
        self.tax_reports = Path(tax_reports)
        self.tmp_dir = Path(tmp_dir)

        # Prefilter parameters
        self.comp_bias_corr = comp_bias_corr
        self.comp_bias_corr_scale = comp_bias_corr_scale
        self.add_self_matches = add_self_matches
        self.seed_sub_mat = seed_sub_mat
        self.s = s
        self.k = k
        self.target_search_mode = target_search_mode
        self.k_score = k_score
        self.alph_size = alph_size
        self.max_seqs = max_seqs
        self.split = split
        self.split_mode = split_mode
        self.split_memory_limit = split_memory_limit
        self.diag_score = diag_score
        self.exact_kmer_matching = exact_kmer_matching
        self.mask = mask
        self.mask_prob = mask_prob
        self.mask_lower_case = mask_lower_case
        self.mask_n_repeat = mask_n_repeat
        self.min_ungapped_score = min_ungapped_score
        self.spaced_kmer_mode = spaced_kmer_mode
        self.spaced_kmer_pattern = spaced_kmer_pattern
        self.local_tmp = local_tmp
        self.disk_space_limit = disk_space_limit

        # Align parameters
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
        self.exhaustive_search_filter = exhaustive_search_filter

        # Profile parameters
        self.pca = pca
        self.pcb = pcb
        self.mask_profile = mask_profile
        self.e_profile = e_profile
        self.wg = wg
        self.filter_msa = filter_msa
        self.filter_min_enable = filter_min_enable
        self.max_seq_id = max_seq_id
        self.qid = qid
        self.qsc = qsc
        self.cov = cov
        self.diff = diff
        self.pseudo_cnt_mode = pseudo_cnt_mode
        self.profile_output_mode = profile_output_mode
        self.exhaustive_search = exhaustive_search
        self.lca_search = lca_search

        # Misc parameters
        self.orf_filter = orf_filter
        self.orf_filter_e = orf_filter_e
        self.orf_filter_s = orf_filter_s
        self.lca_mode = lca_mode
        self.majority = majority
        self.vote_mode = vote_mode
        self.lca_ranks = lca_ranks
        self.tax_lineage = tax_lineage
        self.blacklist = blacklist
        self.taxon_list = taxon_list
        self.prefilter_mode = prefilter_mode
        self.rescore_mode = rescore_mode
        self.allow_deletion = allow_deletion
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
        self.search_type = search_type
        self.translation_mode = translation_mode
        self.report_mode = report_mode
        self.format_mode = format_mode
        self.format_output = format_output
        self.first_seq_as_repr = first_seq_as_repr
        self.target_column = target_column
        self.full_header = full_header
        self.idx_seq_src = idx_seq_src
        self.dbtype = dbtype
        self.shuffle = shuffle
        self.createdb_mode = createdb_mode

        # Common parameters
        self.compressed = compressed
        self.threads = threads
        self.v = v
        self.sub_mat = sub_mat
        self.max_seq_len = max_seq_len
        self.db_load_mode = db_load_mode
        self.gpu = gpu
        self.gpu_server = gpu_server
        self.gpu_server_wait_timeout = gpu_server_wait_timeout
        self.mpi_runner = mpi_runner
        self.force_reuse = force_reuse
        self.remove_tmp_files = remove_tmp_files

        # Expert parameters
        self.filter_hits = filter_hits
        self.sort_results = sort_results
        self.create_lookup = create_lookup
        self.chain_alignments = chain_alignments
        self.merge_query = merge_query
        self.strand = strand
        self.db_output = db_output
        self.write_lookup = write_lookup

        self._defaults = DEFAULTS
        self._path_params = [param for param, info in DEFAULTS.items() if info['type'] == 'path']
        self._caller_dir = get_caller_dir()

    def _validate(self) -> None:
        self._check_required_files()
        self._validate_choices()

        if not (0 <= self.comp_bias_corr_scale <= 1):
            raise ValueError("comp_bias_corr_scale must be between 0 and 1")
        if not (1.0 <= self.s <= 7.5):
            raise ValueError("Sensitivity (-s) must be between 1.0 and 7.5")
        if not (0.0 <= self.min_seq_id <= 1.0):
            raise ValueError("min_seq_id must be between 0.0 and 1.0")
        if not (0.0 <= self.mask_prob <= 1.0):
            raise ValueError("mask_prob must be between 0.0 and 1.0")
        if not (0.0 <= self.majority <= 1.0):
            raise ValueError("majority must be between 0.0 and 1.0")

    def run(self) -> None:
        self._resolve_all_path(self._caller_dir)
        self._validate()

        args = self._get_command_args("easy_taxonomy")
        mmseqs_output = run_mmseqs_command(args)

        self._handle_command_output(
            mmseqs_output=mmseqs_output,
            output_identifier="Easy Taxonomy",
            output_path=str(self.tax_reports)
        )
