# pymmseqs/commands/easy_linclust.py

from pathlib import Path
from typing import List, Union

from ..config import EasyLinClustConfig
from ..parsers import EasyLinClustParser
from ..utils import tmp_dir_handler


def easy_linclust(
    # Required parameters
    fasta_files: Union[str, Path, List[Union[str, Path]]],
    cluster_prefix: Union[str, Path],
    # Optional parameters
    tmp_dir: Union[str, Path, None] = None,
    min_seq_id: float = 0.0,
    c: float = 0.8,
    cov_mode: int = 0,
    e: float = 0.001,
    cluster_mode: int = 0,
    kmer_per_seq: int = 21,
    kmer_per_seq_scale: str = "aa:0.000,nucl:0.200",
    gpu: int = 0,
) -> EasyLinClustParser:
    """
    Perform fast linear-time sequence clustering using MMseqs2.

    This is a faster alternative to easy_cluster that uses linear-time clustering,
    suitable for very large datasets.

    Parameters
    ----------
    `fasta_files` : Union[str, Path, List[Union[str, Path]]]
        Path to a FASTA file (optionally compressed with .gz or .bz2).

    `cluster_prefix` : Union[str, Path]
        Output cluster database path prefix (will create multiple files with this prefix).

    Optional parameters
    --------------------
    `tmp_dir` : Union[str, Path]
        Temporary directory for intermediate files.
        If not provided, a temporary directory will be created in the same directory as the cluster_prefix.

    `min_seq_id` : float, optional
        Minimum sequence identity (range 0.0, 1.0)
        - 0.0 (default)

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

    `e` : float, optional
        E-value threshold (range 0.0, inf)
        - 0.001 (default)

    `cluster_mode` : int, optional
        Clustering method.
        - 0: Set-Cover (greedy) (default)
        - 1: Connected component (BLASTclust)
        - 2, 3: Greedy by sequence length (CDHIT)

    `kmer_per_seq` : int, optional
        Number of k-mers per sequence.
        - 21 (default)

    `kmer_per_seq_scale` : str, optional
        Scale k-mer per sequence based on sequence length as kmer-per-seq val + scale x seqlen.
        - "aa:0.000,nucl:0.200" (default)

    `gpu` : int, optional
        Use GPU (CUDA) if possible
        - 0: off (default)
        - 1: on

    Returns
    -------
    EasyLinClustParser object
        An EasyLinClustParser instance that provides methods to access and parse the clustering results.

    Notes
    -----
    Output files:
        - {cluster_prefix}_rep_seq.fasta: Representative sequences
        - {cluster_prefix}_all_seqs.fasta: FASTA-like per cluster
        - {cluster_prefix}_cluster.tsv: Adjacency list

    For nucleotide sequences, consider using --kmer-per-seq-scale 0.3

    """

    tmp_dir = tmp_dir_handler(tmp_dir=tmp_dir, output_file_path=cluster_prefix)

    config = EasyLinClustConfig(
        fasta_files=fasta_files,
        cluster_prefix=cluster_prefix,
        tmp_dir=tmp_dir,
        min_seq_id=min_seq_id,
        c=c,
        cov_mode=cov_mode,
        e=e,
        cluster_mode=cluster_mode,
        kmer_per_seq=kmer_per_seq,
        kmer_per_seq_scale=kmer_per_seq_scale,
        gpu=gpu,
    )

    config.run()

    return EasyLinClustParser(config)
