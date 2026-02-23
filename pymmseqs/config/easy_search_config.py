# pymmseqs/config/easy_search_config.py

from pathlib import Path
from typing import Union

from .base import BaseConfig
from ..defaults import loader   
from ..utils import (
    get_caller_dir,
    run_mmseqs_command
)

DEFAULTS = loader.load("easy_search")

class EasySearchConfig(BaseConfig):
    def __init__(
        self,
        # Required parameters
        query_fasta: Union[str, Path],
        target_fasta_or_db: Union[str, Path],
        alignment_file: Union[str, Path],
        tmp_dir: Union[str, Path],
        
        # Prefilter parameters
        comp_bias_corr: bool = True,
        comp_bias_corr_scale: float = 1.0,
        add_self_matches: bool = False,
        seed_sub_mat: str = "aa:VTML80.out,nucl:nucleotide.out",
        s: float = 5.7,
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
        
        # Alignment parameters
        a: bool = False,
        alignment_mode: int = 3,
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
        exhaustive_search_filter: bool = False,
        
        # Profile parameters
        pca: float = None,
        pcb: float = None,
        mask_profile: bool = True,
        e_profile: float = 0.001,
        wg: bool = False,
        filter_msa: bool = True,
        filter_min_enable: int = 0,
        max_seq_id: float = 0.9,
        qid: str = "0.0",
        qsc: float = -20.0,
        cov: float = 0.0,
        diff: int = 1000,
        pseudo_cnt_mode: int = 0,
        profile_output_mode: int = 0,
        num_iterations: int = 1,
        exhaustive_search: bool = False,
        lca_search: bool = False,
        
        # Misc parameters
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
        translate: bool = False,
        use_all_table_starts: bool = False,
        id_offset: int = 0,
        sequence_overlap: int = 0,
        sequence_split_mode: int = 1,
        headers_split_mode: int = 0,
        search_type: int = 0,
        start_sens: float = 4.0,
        sens_steps: int = 1,
        translation_mode: int = 0,
        format_mode: int = 0,
        format_output: str = "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits",
        overlap: float = 0.0,
        dbtype: int = 0,
        shuffle: bool = True,
        createdb_mode: int = 0,
        greedy_best_hits: bool = False,
        
        # Common parameters
        sub_mat: str = "aa:blosum62.out,nucl:nucleotide.out",
        max_seq_len: int = 65535,
        db_load_mode: int = 0,
        threads: Union[str, int] = 'all',
        compressed: bool = False,
        v: int = 3,
        gpu: bool = False,
        gpu_server: bool = False,
        gpu_server_wait_timeout: int = 600,
        mpi_runner: str = "",
        force_reuse: bool = False,
        remove_tmp_files: bool = True,
        
        # Expert parameters
        filter_hits: bool = False,
        sort_results: int = 0,
        create_lookup: bool = False,
        chain_alignments: bool = False,
        merge_query: bool = True,
        strand: int = 1,
        db_output: bool = False,
        write_lookup: bool = False,
    ):
        """
        Perform a fast and sensitive sequence search using MMseqs2 easy-search.

        Parameters
        ----------
        `query_fasta` : Union[str, Path]
            Path to one or more query FASTA files. Can be compressed with .gz or .bz2.

        `target_fasta_or_db` : Union[str, Path]
            Path to a target FASTA file (optionally compressed) or an MMseqs2 target database.

        `alignment_file` : Union[str, Path]
            Path to the output file where alignments will be stored.

        `tmp_dir` : Union[str, Path]
            Temporary directory for intermediate files. Will be created if not existing.

        Prefilter Parameters
        -------------------
        `comp_bias_corr` : bool, optional
            Correct for locally biased amino acid composition
            - True (default)
            - False

        `comp_bias_corr_scale` : float, optional
            Scale factor for composition bias correction
            - Range 0, 1
            - 1.0 (default)

        `add_self_matches` : bool, optional
            Artificially add entries of queries with themselves (for clustering)
            - True
            - False (default)

        `seed_sub_mat` : str, optional
            Substitution matrix for k-mer generation as (type:path, type:path)
            type: "aa" or "nucl"
            path: matrix file path
            - "aa:VTML80.out,nucl:nucleotide.out" (default)

            Note: find available matrices in the MMseqs2 data directory: (https://github.com/soedinglab/MMseqs2/tree/master/data)

        `s` : float, optional
            Sensitivity
            - 1.0: faster
            - 4.0: fast
            - 5.7 (default)
            - 7.5: sensitive

        `k` : int, optional
            k-mer length
            - 0: auto (default)

        `target_search_mode` : int, optional
            Target search mode
            - 0: regular k-mer (default)
            - 1: similar k-mer

        `k_score` : str, optional
            k-mer thresholds for sequence and profile searches
            - "seq:2147483647,prof:2147483647" (default)

        `alph_size` : str, optional
            Alphabet sizes for amino acid (protein) and nucleotide sequences (range 2-21)
            - "aa:21,nucl:5" (default)
                - aa:21: 20 amino acids + X for unknown residues
                - nucl:5: 4 nucleotides + N for unknown bases

        `max_seqs` : int, optional
            Maximum results per query passing prefilter
            - 300 (default)
            - Higher values increase sensitivity but may slow down the search

        `split` : int, optional
            Split input into N chunks
            - 0: set the best split automatically (default)

        `split_mode` : int, optional
            Split strategy
            - 0: split target db
            - 1: split query db
            - 2: auto, depending on main memory (default)

        `split_memory_limit` : str, optional
            Maximum memory allocated per split for processing
            - "0":  all available (default)
                - Use suffixes like K, M, or G (e.g., "4G" for 4 gigabytes)

        `diag_score` : bool, optional
            Use ungapped diagonal scoring during prefilter
            - True (default)
            - False

        `exact_kmer_matching` : bool, optional
            Extract only exact k-mers for matching
            - True
            - False (default)

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
            Maximum number of consecutive Ns allowed in a sequence
            - 0 (default)

        `min_ungapped_score` : int, optional
            Minimum ungapped alignment score
            - 15 (default)
            - Higher values increase specificity but may reduce sensitivity

        `spaced_kmer_mode` : int, optional
            Spaced k-mer mode
            - 0: consecutive
            - 1: spaced (default)

        `spaced_kmer_pattern` : str, optional
            Custom pattern for spaced k-mers used during k-mer matching.
            - Define a pattern of 1s (match positions) and 0s (ignore positions)
            - Example: "1101011" means 5 match positions and 2 ignored positions
            - Increases sensitivity by focusing on conserved regions while allowing flexibility in less conserved areas.

        `local_tmp` : str, optional
            Path to an alternative temporary directory for storing intermediate files
            - Useful for reducing I/O load on shared storage systems (e.g., NFS)
            - Default: Temporary files are stored in the main tmpDir

        `disk_space_limit` : str, optional
            Max disk space usage
            - "0": unlimited (default)
            - Use suffixes like K, M, or G (e.g., "100G" for 100 gigabytes)

        Alignment Parameters
        --------------------
        `a` : bool, optional
            Add backtrace string (convert to alignments with mmseqs convertalis module)
            - True
            - False (default)

        `alignment_mode` : int, optional
            Alignment detail level
            - 0: auto
            - 1: score + end_po
            - 2: + start_pos + cov
            - 3: + seq.id (default)
            - 4: only ungapped alignment

        `alignment_output_mode` : int, optional
            Output detail level
            - 0: auto (default)
            - 1: score + end_po
            - 2: + start_pos + cov
            - 3: + seq.id
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
            - 0.0 (default)
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

        `corr_score_weight` : float, optional
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

        `exhaustive_search_filter` : bool, optional
            Filter result during search
            - True
            - False (default)

        Profile Parameters
        ------------------
        `pca` : float, optional
            Pseudo count admixture strength for profile construction
            - 0.0 (default)
            - Higher values increase the weight of pseudo counts, making the profile more conservative
            - Lower values reduce their influence, making the profile more specific to the input sequences

        `pcb` : float, optional
            Controls the threshold for pseudo-count admixture based on the effective number of sequences (Neff) (range 0.0, inf)
            - 0.0 (default)
            - Lower values apply pseudo-counts more aggressively

        `mask_profile` : bool, optional
            Mask low-complexity regions in the query sequence of a profile using TANTAN
            - True (default)
            - False

        `e_profile` : float, optional
            E-value threshold for including sequence matches in the profile
            - 0.1 (default)

        `wg` : bool, optional
            Use global sequence weighting for profile calculation
            - True
            - False (default)

        `filter_msa` : bool, optional
            Filter MSA
            - True (default)
            - False

        `filter_min_enable` : int, optional
            Minimum number of sequences required to trigger filtering of MSAs
            - 0: Always filter (default)
            - N > 0: Filter only if the MSA contains more than N sequences

        `max_seq_id` : float, optional
            Maximum pairwise sequence identity for redundancy reduction in the output MSA (range 0.0, 1.0)
            - 0.9 (default)
            - Filters sequences to ensure no two sequences in the output share more than the specified identity

        `qid` : str, optional
            Filters output MSAs by minimum sequence identity with the query (range 0.0, 1.0)
            - 0.0: no filtering (default)
            - Can specify multiple thresholds as a comma-separated list (e.g., "0.15,0.30,0.50") to create filter buckets
                - Example: "0.15,0.30,0.50" creates buckets for sequences with identities in ]0.15-0.30] and ]0.30-0.50]

        `qsc` : float, optional
            Filters output MSAs by minimum score per aligned residue with query sequences (range -50.0, 100.0)
            - -20.0 (default)
            - Higher values reduce diversity in the output MSAs by retaining only high-scoring alignments

        `cov` : float, optional
            Minimum fraction of query residues covered by matched sequences to filter output MSAs (range 0.0, 1.0)
            - 0.0 (default)

        `diff` : int, optional
            Filters MSAs by selecting the most diverse sequences, ensuring at least this many sequences are kept in each MSA block of length 50
            - 1000 (default)

        `pseudo_cnt_mode` : int, optional
            Pseudocount method
            - 0: substitution-matrix (default)
            - 1: context-specific pseudocounts
        
        `profile_output_mode` : int, optional
            Output mode for profile search results
            - 0: binary log-odds (default)
            - 1: human readable frequencies
        
        `num_iterations` : int, optional
            Number of iterative profile search iterations
            - 1: (default)

        `exhaustive_search` : bool, optional
            For bigger profile DB, run iteratively the search by greedily swapping the search results
            - True
            - False (default)

        `lca_search` : bool, optional
            Enable LCA candidate search
            - True
            - False (default)

        Misc Parameters
        ---------------
        `taxon_list` : str, optional
            Taxonomy IDs to filter results by. Multiple IDs can be provided, separated by commas (no spaces)
            - "" (default)
            - Example: "9606,10090"

        `prefilter_mode` : int, optional
            Prefilter method
            - 0: kmer/ungapped (default)
            - 1: ungapped
            - 2: nofilter
            - 3: ungapped+gapped

        `rescore_mode` : int, optional
            Rescore diagonals with:
            - 0: Hamming distance (default)
            - 1: local alignment (score only)
            - 2: local alignment
            - 3: global alignment
            - 4: longest alignment fulfilling window quality criterion

        `allow_deletion` : bool, optional
            Allow deletions in MSA
            - True
            - False (default)

        `min_length` : int, optional
            Minimum codon number in open reading frames (ORFs)
            - 30 (default)

        `max_length` : int, optional
            Maximum codon number in open reading frames (ORFs)
            - 32734 (default)

        `max_gaps` : int, optional
            Maximum number of codons with gaps or unknown residues before an open reading frame is rejected
            - 2147483647 (default)

        `contig_start_mode` : int, optional
            Contig start handling
            - 0: incomplete
            - 1: complete
            - 2: both (default)

        `contig_end_mode` : int, optional
            Contig end handling
            - 0: incomplete
            - 1: complete
            - 2: both (default)

        `orf_start_mode` : int, optional
            ORF start handling
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

        `translate` : bool, optional
            Translate open reading frames (ORFs) to amino acid
            - True
            - False (default)

        `use_all_table_starts` : bool, optional
            Use all start codons
            - True
            - False: only ATG (AUG) (default)

        `id_offset` : int, optional
            Numeric IDs in index file are offset by this value
            - 0 (default)

        `sequence_overlap` : int, optional
            Overlap between sequences
            - 0 (default)

        `sequence_split_mode` : int, optional
            Method for splitting sequences during processing
            - 0: Copy data (creates a full copy of the sequence data).
            - 1: Soft link data and write a new index (saves disk space by linking to the original data) (default)

        `headers_split_mode` : int, optional
            Header split method
            - 0: Split positions (Headers are split based on predefined positions) (default)
            - 1: Original header (Headers are preserved as-is without splitting)

        `search_type` : int, optional
            Search mode:
            - 0: auto (default)
            - 1: amino
            - 2: translated
            - 3: nucleotide
            - 4: translated alignment

        `start_sens` : float, optional
            Initial sensitivity
            - 4.0 (default)

        `sens_steps` : int, optional
            Number of search steps performed from `start_sens` argument to `s` argument
            - 1 (default)

        `translation_mode` : int, optional
            Translation AA seq from nucletoide method
            - 0: Open Reading Frames (ORFs) (default)
            - 1: Full Reading Frames
        
        `format_mode` : int, optional  
            Output format type  
            - 0: BLAST-TAB (default)  
            - 1: SAM  
            - 2: BLAST-TAB + query/db length  
            - 3: Pretty HTML  
            - 4: BLAST-TAB + column headers  

            Notes:  
            - BLAST-TAB (0) and BLAST-TAB + column headers (4) support custom output formats via `format_output`.  

        `format_output` : str, optional  
            Comma-separated list of output columns to include in results.  
            Available columns:  
            - query, target, evalue, gapopen, pident, fident, nident, qstart, qend, qlen  
            - tstart, tend, tlen, alnlen, raw, bits, cigar, qseq, tseq, qheader, theader, qaln, taln  
            - qframe, tframe, mismatch, qcov, tcov, qset, qsetid, tset, tsetid, taxid, taxname, taxlineage  
            - qorfstart, qorfend, torfstart, torfend, ppos  

            - Default: "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"

        `overlap` : float, optional  
            Maximum allowed overlap of covered regions between query and target sequences.  
            - 0.0 (default)  
            - Higher values allow more overlap, reducing redundancy.  

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
            - 0: Copy data (default)  
            - 1: Soft link data and write a new index (only works with single-line FASTA/Q files)  

            Notes:  
            - Using soft links reduces disk usage by referencing the original sequence data instead of duplicating it.  

        `greedy_best_hits` : bool, optional  
            Select best hits greedily to maximize query coverage.  
            - True  
            - False (default)
        
        Common Parameters
        ----------------
        `sub_mat` : str, optional
            Substitution matrix (type:path, type:path)
            type: "aa" or "nucl"
            path: matrix file path
            - Default: "aa:blosum62.out,nucl:nucleotide.out"

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

        `gpu` : bool, optional
            Use GPU (CUDA) if possible
            - True
            - False (default)

        `gpu_server` : bool, optional
            Use GPU server
            - True
            - False (default)
        
        `gpu_server_wait_timeout` : int, optional
            - 0: do not wait (default)
            - -1: no wait limit
            - >0: wait for N seconds

        `mpi_runner` : str, optional
            Use MPI on compute cluster with this MPI command (e.g., "mpirun -np 42")
            - "" (default)

        `force_reuse` : bool, optional
            Reuse tmp filse in tmp/latest folder ignoring parameters and version changes
            - True
            - False (default)

        `remove_tmp_files` : bool, optional
            Delete temporary files
            - True (default)
            - False

        Expert Parameters
        ----------------
        `filter_hits` : bool, optional
            Filter hits by sequence ID and coverage
            - True
            - False (default)

        `sort_results` : int, optional
            Result sorting method
            - 0: No sorting (default)
            - 1: E-value (Alignment) or sequence ID (Hamming) 

        `create_lookup` : bool, optional
            Create lookup file (can be very large)
            - True
            - False (default)

        `chain_alignments` : bool, optional
            Chain overlapping alignments
            - True
            - False (default)

        `merge_query` : bool, optional
            Combine ORFs/split sequences to a single entry
            - True (default)
            - False

        `strand` : int, optional
            Strand selection (only works for DNA/DNA search)
            - 0: reverse
            - 1: forward (default)
            - 2: both
        
        `db_output` : bool, optional
            Return the result as DB instead of a text file
            - True: DB output
            - False: Text file output (default)
        
        `write_lookup` : bool, optional
            Create a `.lookup` file mapping internal IDs to FASTA IDs
            - True (default)
            - False

        Returns
        -------
        None
            Writes alignment output to specified file.

        Raises
        ------
        FileNotFoundError
            If input FASTA files or databases are missing.

        ValueError
            For invalid parameter combinations.

        Examples
        --------
        Basic protein sequence search:
        >>> easy_search(
            query_fasta="query.fasta",
            target_fasta_or_db="target.fasta",
            alignment_output="output.m8",
            tmp_dir="tmp_search",
            s=5.7,
            threads=8
        )

        Iterative profile search:
        >>> easy_search(
            query_fasta="query.fasta",
            target_fasta_or_db="database.fasta",
            alignment_output="profile_search.m8",
            tmp_dir="tmp_search",
            num_iterations=3
        )

        Search against a small profile database:
        >>> easy_search(
            query_fasta="query.fasta",
            target_fasta_or_db="pfam_db",
            alignment_output="pfam_results.m8",
            tmp_dir="tmp_search"
        )

        """
        super().__init__()
        
        # Initialize the config manually
        self.query_fasta = query_fasta if isinstance(query_fasta, list) else [query_fasta]
        self.target_fasta_or_db = Path(target_fasta_or_db)
        self.alignment_file = Path(alignment_file)
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
        
        # Alignment parameters
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
        self.num_iterations = num_iterations
        self.exhaustive_search = exhaustive_search
        self.lca_search = lca_search
        
        # Misc parameters
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
        self.start_sens = start_sens
        self.sens_steps = sens_steps
        self.translation_mode = translation_mode
        self.format_mode = format_mode
        self.format_output = format_output
        self.overlap = overlap
        self.dbtype = dbtype
        self.shuffle = shuffle
        self.createdb_mode = createdb_mode
        self.greedy_best_hits = greedy_best_hits
        
        # Common parameters
        self.sub_mat = sub_mat
        self.max_seq_len = max_seq_len
        self.db_load_mode = db_load_mode
        self.threads = threads
        self.compressed = compressed
        self.v = v
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
            
        # Validate numeric constraints
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
        if not (0 <= self.id_offset):
            raise ValueError("id_offset must be non-negative")
        if not (0 <= self.threads):
            raise ValueError("threads must be >= 0")
        if not (0.0 <= self.max_seq_id <= 1.0):
            raise ValueError("max_seq_id must be between 0.0 and 1.0")
        if not (-50.0 <= self.qsc <= 100.0):
            raise ValueError("qsc must be between -50.0 and 100.0")
        if not (0.0 <= self.cov <= 1.0):
            raise ValueError("cov must be between 0.0 and 1.0")
        if not (0 <= self.min_length <= self.max_length):
            raise ValueError("min_length must be >= 0 and <= max_length")
        if not (0 <= self.max_gaps):
            raise ValueError("max_gaps must be >= 0")
        if not (0 <= self.num_iterations):
            raise ValueError("num_iterations must be >= 0")

    def run(self) -> None:
        self._resolve_all_path(self._caller_dir)

        self._validate()
        
        args = self._get_command_args("easy_search")
        mmseqs_output = run_mmseqs_command(args)
        
        self._handle_command_output(
            mmseqs_output=mmseqs_output,
            output_identifier="Easy Search",
            output_path=str(self.alignment_file)
        )
