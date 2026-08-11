# pymmseqs/parsers/base_cluster_parser.py

import os
import numpy as np
import pandas as pd
from math import ceil
from typing import Generator, Optional, Union

from ..tools.easy_cluster_tools import parse_fasta_clusters
from ..utils import write_fasta


def _train_test_split(items, test_size, shuffle=True, random_state=None):
    """
    Split a list into (train, test), matching sklearn's train_test_split.

    Built on numpy, which is already a dependency, so scikit-learn is not
    pulled in for the two calls below. scikit-learn publishes no musllinux
    wheels, which made `pip install pymmseqs` fail on Alpine.

    Verified against sklearn on 780 combinations of length, test_size, seed
    and shuffle. Two details matter for that equivalence:
      - n_train is n - n_test, not floor(n * (1 - test_size)). The latter
        disagrees at proportions like 0.9, where 1 - 0.9 is 0.09999999999999998.
      - the permutation must come from RandomState, not the newer Generator
        API, since that is what sklearn uses to honour random_state.

    Unlike sklearn this does not raise when a resulting split is empty, or
    when random_state is set alongside shuffle=False; callers here already
    handle empty splits.
    """
    n = len(items)
    n_test = ceil(n * test_size)
    n_train = n - n_test

    if not shuffle:
        return items[:n_train], items[n_train:]

    perm = np.random.RandomState(random_state).permutation(n)
    test_idx, train_idx = perm[:n_test], perm[n_test:n_test + n_train]
    return [items[i] for i in train_idx], [items[i] for i in test_idx]


class BaseClusterParser:
    """
    Base class for parsing the output of MMseqs2 clustering commands
    (easy-cluster and easy-linclust). Provides shared methods for accessing
    cluster results in multiple formats.
    """
    def __init__(
        self,
        config,
        seq_id_separator: str = "|",
        seq_id_index: int = 1
    ):
        """
        Parameters
        ----------
        config
            A clustering configuration object with a `cluster_prefix` attribute.
        seq_id_separator: str, optional
            The separator used in the FASTA headers to separate the sequence ID from other information.
            Default is "|".
        seq_id_index: int, optional
            The index of the sequence ID in the FASTA header.
            Default is 1.
        Note: It tries to extract the seq_id from the header using the separator and index, if it fails, it doesn't add the seq_id to the member.
        """
        self.cluster_prefix = config.cluster_prefix
        self.seq_id_separator = seq_id_separator
        self.seq_id_index = seq_id_index

    def split_rep_as_fasta(
        self,
        train: float,
        val: float,
        test: float,
        shuffle: bool = True,
        seed: int = None
    ) -> tuple[Optional[str], Optional[str], Optional[str]]:
        """
        Splits the clusters into train, validation, and test sets.

        - if train + val + test != 1.0, the proportions are normalized to sum to 1.0

        Parameters
        ----------
        train: float
            The proportion of the data to use for training.
        val: float
            The proportion of the data to use for validation.
        test: float
            The proportion of the data to use for testing.
        shuffle: bool, optional
            Whether to shuffle the data before splitting.
            Default is True.
        seed: int, optional
            The seed for the random number generator.
            Default is None.

        Returns
        -------
        Path
            The path to the train file, or None if the train set is empty.
        Path
            The path to the validation file, or None if the validation set is empty.
        Path
            The path to the test file, or None if the test set is empty.
        """
        train_reps, val_reps, test_reps = self.split_rep_as_list(
            train=train,
            val=val,
            test=test,
            with_seq=True,
            shuffle=shuffle,
            seed=seed
        )

        parent_dir = os.path.dirname(self.cluster_prefix)
        if not parent_dir:
            parent_dir = "."

        base_name = os.path.basename(self.cluster_prefix)

        train_path = None
        val_path = None
        test_path = None

        if train_reps:
            train_path = os.path.join(parent_dir, f"{base_name}_rep_train.fasta")
            write_fasta(train_reps, train_path)

        if val_reps:
            val_path = os.path.join(parent_dir, f"{base_name}_rep_val.fasta")
            write_fasta(val_reps, val_path)

        if test_reps:
            test_path = os.path.join(parent_dir, f"{base_name}_rep_test.fasta")
            write_fasta(test_reps, test_path)

        return train_path, val_path, test_path

    def split_rep_as_list(
        self,
        train: float,
        val: float,
        test: float,
        with_seq: bool = True,
        shuffle: bool = True,
        seed: int = None
    ) -> tuple[list, list, list]:
        """
        Splits the representatives to train, validation, and test sets.

        - if train + val + test != 1.0, the proportions are normalized to sum to 1.0

        Parameters
        ----------
        train: float
            The proportion of the data to use for training.
        val: float
            The proportion of the data to use for validation.
        test: float
            The proportion of the data to use for testing.
        with_seq: bool, optional
            If True, returns a list of tuples with the representative sequence.
            Default is True.
        shuffle: bool, optional
            Whether to shuffle the data before splitting.
            Default is True.
        seed: int, optional
            The seed for the random number generator.
            Default is None.

        Returns
        -------
        tuple[list, list, list]
            A tuple of lists containing the train, validation, and test sets.
        """
        rep_seqs = self.to_rep_list(
            with_seq=with_seq
        )

        if train + val + test != 1.0:
            total = train + val + test
            train = train / total
            val = val / total
            test = test / total

        if val == 0 and test == 0:
            return rep_seqs, [], []

        train_rep_seqs, temp_rep_seqs = _train_test_split(
            rep_seqs,
            test_size=(val + test),
            shuffle=shuffle,
            random_state=seed
        )

        if val == 0:
            return train_rep_seqs, [], temp_rep_seqs
        elif test == 0:
            return train_rep_seqs, temp_rep_seqs, []

        test_proportion_in_temp = test / (val + test)
        val_rep_seqs, test_rep_seqs = _train_test_split(
            temp_rep_seqs,
            test_size=test_proportion_in_temp,
            shuffle=shuffle,
            random_state=seed
        )

        return train_rep_seqs, val_rep_seqs, test_rep_seqs

    def to_rep_list(self, with_seq: bool = True) -> list[Union[tuple[str, str], str]]:
        """
        Returns a list of representatives.

        Parameters
        ----------
        with_seq: bool, optional
            If True, returns a list of tuples with the representative sequence.
            Default is True.

        Returns
        -------
        list[Union[tuple[str, str], str]]
            A list of representatives.
        """
        gen = self.to_gen()
        rep_seqs = []

        if with_seq:
            for cluster in gen:
                rep = cluster["rep"]
                rep_seq = cluster["members"][0]["sequence"]
                rep_seqs.append((rep, rep_seq))
        else:
            for cluster in gen:
                rep_seqs.append(cluster["rep"])

        return rep_seqs

    def to_rep_gen(self, with_seq: bool = True) -> Generator:
        """
        Returns a generator of representatives.

        Parameters
        ----------
        with_seq: bool, optional
            If True, returns a generator of tuples with the representative sequence.
            Default is True.

        Returns
        -------
        Generator
            A generator of representatives.
        """
        gen = self.to_gen()

        if with_seq:
            for cluster in gen:
                rep = cluster["rep"]
                rep_seq = cluster["members"][0]["sequence"]
                yield (rep, rep_seq)
        else:
            for cluster in gen:
                yield cluster["rep"]

    def to_list(self) -> list:
        """
        Parses a FASTA file containing clustered sequences and returns a list of dictionaries,
        where each dictionary represents a cluster.

        Returns:
        --------
        list of dict
            A list of dictionaries where each dictionary represents a single cluster with the following keys:
            - "rep": The representative sequence ID.
            - "members": List of member dictionaries in the cluster with the following keys:
                - "seq_id": Unique sequence identifier extracted from the header.
                    - If the header has format like ">seq_id|header", the seq_id is extracted from the header.
                - "header": Full FASTA header for the sequence.
                - "sequence": Nucleotide or protein sequence.

        When to Use:
        ------------
        - When you need to preserve the order of clusters as they appear in the file.
        - When you need to process all clusters at once and memory usage is not a concern.
        """
        return [
            {
            "rep": rep,
            "members": members
            }
            for rep, members in parse_fasta_clusters(f"{self.cluster_prefix}_all_seqs.fasta", self.seq_id_separator, self.seq_id_index)
        ]

    def to_pandas(self) -> pd.DataFrame:
        clusters = self.to_list()
        rows = []
        for cluster in clusters:
            rep = cluster["rep"]
            for member in cluster["members"]:
                # If the member has a seq_id, add it to the rows
                if "seq_id" in member:
                    rows.append({
                        "rep": rep,
                        "seq_id": member["seq_id"],
                        "header": member["header"],
                        "sequence": member["sequence"]
                    })
                else:
                    rows.append({
                        "rep": rep,
                        "header": member["header"],
                        "sequence": member["sequence"]
                    })
        return pd.DataFrame(rows).set_index('rep')

    def to_gen(self) -> Generator:
        """
        Generator that yields clusters one at a time from a FASTA file as dictionaries.

        Yields:
        -------
        dict
            A dictionary which represents a single cluster with the following keys:
            - "rep": The representative sequence ID.
            - "members": List of member dictionaries in the cluster with the following keys:
                - "seq_id": Unique sequence identifier extracted from the header.
                - "header": Full FASTA header for the sequence.
                - "sequence": Nucleotide or protein sequence.

        When to Use:
        ------------
        - When processing very large files where loading all clusters at once would consume too much memory.
        - When implementing streaming pipelines that process one cluster at a time.
        - When you need a dictionary format but want to avoid loading the entire dataset into memory.
        """
        for rep, members in parse_fasta_clusters(f"{self.cluster_prefix}_all_seqs.fasta", self.seq_id_separator, self.seq_id_index):
            yield {
                "rep": rep,
                "members": members
            }

    def to_path(self) -> list[str]:
        """
        Returns a list of file paths for the output files.

        Returns:
        --------
        list of str
        """
        return [
            f"{self.cluster_prefix}_all_seqs.fasta",
            f"{self.cluster_prefix}_cluster.tsv",
            f"{self.cluster_prefix}_rep_seq.fasta",
        ]
