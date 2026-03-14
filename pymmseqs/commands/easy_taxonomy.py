# pymmseqs/commands/easy_taxonomy.py

from pathlib import Path
from typing import List, Optional, Union

from ..config import EasyTaxonomyConfig
from ..parsers import EasyTaxonomyParser
from ..utils import tmp_dir_handler


def easy_taxonomy(
    # Required parameters
    fasta_file: Union[str, Path, List[Union[str, Path]]],
    target_db: Union[str, Path],
    tax_reports: Union[str, Path],

    # Optional parameters
    tmp_dir: Union[str, Path, None] = None,

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
) -> EasyTaxonomyParser:
    """Assign taxonomic labels to FASTA sequences using MMseqs2.

    Parameters
    ----------
    `fasta_file` : Union[str, Path, List[Union[str, Path]]]
        Path(s) to the input FASTA file(s), optionally compressed as .gz or .bz2.

    `target_db` : Union[str, Path]
        Target sequence database with taxonomy information (e.g., UniProtKB/Swiss-Prot).

    `tax_reports` : Union[str, Path]
        Output prefix for taxonomy report files.

    `tmp_dir` : Union[str, Path, None], optional
        Temporary directory for intermediate files.
        If not provided, a temporary directory will be created next to the output.

    `lca_mode` : int, optional
        LCA Mode 1: single search LCA, 2/3: approximate 2bLCA, 4: top hit.
        - Default: 3

    `s` : float, optional
        Sensitivity: 1.0 faster; 4.0 fast; 7.5 sensitive.
        - Default: 4.0

    `e` : float, optional
        E-value threshold (range 0.0-inf).
        - Default: 0.001

    `orf_filter` : int, optional
        Prefilter query ORFs with non-selective search.
        Consider disabling (0) when classifying short reads.
        - Default: 0

    `report_mode` : int, optional
        Taxonomy report mode 0: Kraken, 1: Krona, 2: no report, 3: Kraken per query database.
        - Default: 0

    Returns
    -------
    EasyTaxonomyParser
        Parser providing data access, analysis, statistics, and visualization
        methods for taxonomy results.

    Notes
    -----
    Output files created:
        - {tax_reports}_tophit_aln: top hits alignment
        - {tax_reports}_tophit_report: coverage profiles per database entry
        - {tax_reports}_report: kraken-style report
    """
    tmp_dir = tmp_dir_handler(tmp_dir=tmp_dir, output_file_path=tax_reports)

    config = EasyTaxonomyConfig(
        fasta_file=fasta_file,
        target_db=target_db,
        tax_reports=tax_reports,
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
        exhaustive_search=exhaustive_search,
        lca_search=lca_search,
        orf_filter=orf_filter,
        orf_filter_e=orf_filter_e,
        orf_filter_s=orf_filter_s,
        lca_mode=lca_mode,
        majority=majority,
        vote_mode=vote_mode,
        lca_ranks=lca_ranks,
        tax_lineage=tax_lineage,
        blacklist=blacklist,
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
        translation_mode=translation_mode,
        report_mode=report_mode,
        format_mode=format_mode,
        format_output=format_output,
        first_seq_as_repr=first_seq_as_repr,
        target_column=target_column,
        full_header=full_header,
        idx_seq_src=idx_seq_src,
        dbtype=dbtype,
        shuffle=shuffle,
        createdb_mode=createdb_mode,
        compressed=compressed,
        threads=threads,
        v=v,
        sub_mat=sub_mat,
        max_seq_len=max_seq_len,
        db_load_mode=db_load_mode,
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

    return EasyTaxonomyParser(config)
