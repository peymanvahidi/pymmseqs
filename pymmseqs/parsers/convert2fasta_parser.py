# pymmseqs/parsers/convert2fasta_parser.py

import pandas as pd
from typing import Generator

from ..config import Convert2FastaConfig


def _iter_fasta(path) -> Generator[dict, None, None]:
    """Yield {'header': ..., 'sequence': ...} records from a FASTA file."""
    header, seq_lines = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield {"header": header, "sequence": "".join(seq_lines)}
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            yield {"header": header, "sequence": "".join(seq_lines)}


class Convert2FastaParser:
    """
    Parser for convert2fasta output.

    Provides access to the generated FASTA file as a path, a streaming generator,
    a list of records, or a pandas DataFrame.

    Parameters
    ----------
    config : Convert2FastaConfig
        Configuration object after run() has been called.
    """

    def __init__(self, config: Convert2FastaConfig):
        self.fasta_file = config.fasta_file

    def to_path(self) -> str:
        """Return the path to the output FASTA file."""
        return str(self.fasta_file)

    def to_gen(self) -> Generator[dict, None, None]:
        """Yield one {'header', 'sequence'} dict per FASTA record."""
        return _iter_fasta(self.fasta_file)

    def to_list(self) -> list[dict]:
        """Return all FASTA records as a list of {'header', 'sequence'} dicts."""
        return list(_iter_fasta(self.fasta_file))

    def to_pandas(self) -> pd.DataFrame:
        """Return the FASTA records as a DataFrame with 'header' and 'sequence' columns."""
        return pd.DataFrame(self.to_list(), columns=["header", "sequence"])

    def __len__(self) -> int:
        """Return the number of sequences in the FASTA file."""
        return sum(1 for _ in _iter_fasta(self.fasta_file))
