# pymmseqs/commands/easy_rbh.py

from pathlib import Path
from typing import Optional, Union

from ..config import EasyRbhConfig
from ..parsers import GenericParser
from ..utils import tmp_dir_handler


def easy_rbh(
    query_fasta: Union[str, Path],
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

    # Prefilter parameters
    comp_bias_corr: bool = True,
    comp_bias_corr_scale: float = 1.0,
    add_self_matches: bool = False,
    seed_sub_mat: str = "aa:VTML80.out,nucl:nucleotide.out",
    k: int = 0,
    target_search_mode: int = 0,
    k_score: str = "seq:2147483647,prof:2147483647",
    alph_size: str = "aa:21,nucl:5",
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
    min_aln_len: int = 0,
    seq_id_mode: int = 0,
    alt_ali: int = 0,
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
    pca: Optional[float] = None,
    pcb: Optional[float] = None,
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
    use_all_table_starts: bool = False,
    id_offset: int = 0,
    sequence_overlap: int = 0,
    sequence_split_mode: int = 1,
    headers_split_mode: int = 0,
    start_sens: float = 4.0,
    sens_steps: int = 1,
    format_mode: int = 0,
    overlap: float = 0.0,
    dbtype: int = 0,
    shuffle: bool = True,
    createdb_mode: int = 1,
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

) -> GenericParser:
    """Assign reciprocal best hits between query and target sequences.

    Parameters
    ----------
    `query_fasta` : Union[str, Path]
        Path to the input query FASTA file. Can be compressed with .gz or .bz2.

    `target_fasta_or_db` : Union[str, Path]
        Path to a target FASTA file (optionally compressed) or an MMseqs2 target database.

    `alignment_file` : Union[str, Path]
        Path to the output file where reciprocal best hit alignments will be stored.

    `tmp_dir` : Union[str, Path, None], optional
        Temporary directory for intermediate files.
        If not provided, a temporary directory will be created next to the alignment_file.

    `s` : float, optional
        Sensitivity: 1.0 faster, 4.0 fast, 5.7 default, 7.5 sensitive.

    `e` : float, optional
        E-value threshold (range 0.0-inf). Default: 0.001.

    `min_seq_id` : float, optional
        Minimum sequence identity (range 0.0-1.0). Default: 0.0.

    `c` : float, optional
        Coverage threshold. Default: 0.0.

    `max_seqs` : int, optional
        Maximum results per query passing prefilter. Default: 300.

    `translate` : bool, optional
        Translate ORF to amino acid. Default: False.

    `translation_table` : int, optional
        Genetic code table. Default: 1 (Canonical).

    `translation_mode` : int, optional
        Translation method: 0 ORFs (default), 1 full reading frames.

    `search_type` : int, optional
        Search mode: 0 auto (default), 1 amino acid, 2 translated, 3 nucleotide, 4 translated alignment.

    `format_output` : str, optional
        Comma-separated list of output columns.

    Returns
    -------
    GenericParser
        Parser providing access to command output path.
    """
    tmp_dir = tmp_dir_handler(
        tmp_dir=tmp_dir,
        output_file_path=alignment_file
    )

    config = EasyRbhConfig(
        query_fasta=query_fasta,
        target_fasta_or_db=target_fasta_or_db,
        alignment_file=alignment_file,
        tmp_dir=tmp_dir,
        comp_bias_corr=comp_bias_corr,
        comp_bias_corr_scale=comp_bias_corr_scale,
        add_self_matches=add_self_matches,
        seed_sub_mat=seed_sub_mat,
        s=s,
        k=k,
        target_search_mode=target_search_mode,
        k_score=k_score,
        alph_size=alph_size,
        max_seqs=max_seqs,
        split=split,
        split_mode=split_mode,
        split_memory_limit=split_memory_limit,
        diag_score=diag_score,
        exact_kmer_matching=exact_kmer_matching,
        mask=mask,
        mask_prob=mask_prob,
        mask_lower_case=mask_lower_case,
        mask_n_repeat=mask_n_repeat,
        min_ungapped_score=min_ungapped_score,
        spaced_kmer_mode=spaced_kmer_mode,
        spaced_kmer_pattern=spaced_kmer_pattern,
        local_tmp=local_tmp,
        disk_space_limit=disk_space_limit,
        a=a,
        alignment_mode=alignment_mode,
        alignment_output_mode=alignment_output_mode,
        wrapped_scoring=wrapped_scoring,
        e=e,
        min_seq_id=min_seq_id,
        min_aln_len=min_aln_len,
        seq_id_mode=seq_id_mode,
        alt_ali=alt_ali,
        c=c,
        cov_mode=cov_mode,
        max_rejected=max_rejected,
        max_accept=max_accept,
        score_bias=score_bias,
        realign=realign,
        realign_score_bias=realign_score_bias,
        realign_max_seqs=realign_max_seqs,
        corr_score_weight=corr_score_weight,
        gap_open=gap_open,
        gap_extend=gap_extend,
        zdrop=zdrop,
        exhaustive_search_filter=exhaustive_search_filter,
        pca=pca,
        pcb=pcb,
        mask_profile=mask_profile,
        e_profile=e_profile,
        wg=wg,
        filter_msa=filter_msa,
        filter_min_enable=filter_min_enable,
        max_seq_id=max_seq_id,
        qid=qid,
        qsc=qsc,
        cov=cov,
        diff=diff,
        pseudo_cnt_mode=pseudo_cnt_mode,
        profile_output_mode=profile_output_mode,
        num_iterations=num_iterations,
        exhaustive_search=exhaustive_search,
        lca_search=lca_search,
        taxon_list=taxon_list,
        prefilter_mode=prefilter_mode,
        rescore_mode=rescore_mode,
        allow_deletion=allow_deletion,
        min_length=min_length,
        max_length=max_length,
        max_gaps=max_gaps,
        contig_start_mode=contig_start_mode,
        contig_end_mode=contig_end_mode,
        orf_start_mode=orf_start_mode,
        forward_frames=forward_frames,
        reverse_frames=reverse_frames,
        translation_table=translation_table,
        translate=translate,
        use_all_table_starts=use_all_table_starts,
        id_offset=id_offset,
        sequence_overlap=sequence_overlap,
        sequence_split_mode=sequence_split_mode,
        headers_split_mode=headers_split_mode,
        search_type=search_type,
        start_sens=start_sens,
        sens_steps=sens_steps,
        translation_mode=translation_mode,
        format_mode=format_mode,
        format_output=format_output,
        overlap=overlap,
        dbtype=dbtype,
        shuffle=shuffle,
        createdb_mode=createdb_mode,
        greedy_best_hits=greedy_best_hits,
        sub_mat=sub_mat,
        max_seq_len=max_seq_len,
        db_load_mode=db_load_mode,
        threads=threads,
        compressed=compressed,
        v=v,
        gpu=gpu,
        gpu_server=gpu_server,
        gpu_server_wait_timeout=gpu_server_wait_timeout,
        mpi_runner=mpi_runner,
        force_reuse=force_reuse,
        remove_tmp_files=remove_tmp_files,
        filter_hits=filter_hits,
        sort_results=sort_results,
        create_lookup=create_lookup,
        chain_alignments=chain_alignments,
        merge_query=merge_query,
        strand=strand,
        db_output=db_output,
        write_lookup=write_lookup,
    )

    config.run()

    return GenericParser(config)
