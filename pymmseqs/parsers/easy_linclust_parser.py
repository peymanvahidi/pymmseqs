# pymmseqs/parsers/easy_linclust_parser.py

from ..config import EasyLinClustConfig
from .base_cluster_parser import BaseClusterParser

class EasyLinClustParser(BaseClusterParser):
    """
    A class for parsing the output of the EasyLinClustConfig.
    """
    def __init__(
        self,
        config: EasyLinClustConfig,
        seq_id_separator: str = "|",
        seq_id_index: int = 1
    ):
        """
        Parameters
        ----------
        config: EasyLinClustConfig
            The configuration object for the EasyLinClust command.
        seq_id_separator: str, optional
            The separator used in the FASTA headers to separate the sequence ID from other information.
            Default is "|".
        seq_id_index: int, optional
            The index of the sequence ID in the FASTA header.
            Default is 1.
        Note: It tries to extract the seq_id from the header using the separator and index, if it fails, it doesn't add the seq_id to the member.
        """
        super().__init__(config, seq_id_separator, seq_id_index)
