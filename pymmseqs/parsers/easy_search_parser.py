# pymmseqs/parsers/easy_search_parser.py

import pandas as pd
import csv
from typing import Generator, Union
import numpy as np

from ..config import EasySearchConfig, ConvertAlisConfig

class EasySearchParser:
    """
    A class for parsing the output of the EasySearchConfig.
    """
    def __init__(self, config: Union[EasySearchConfig, ConvertAlisConfig]):
        if not config.format_mode == 4:
            raise ValueError(f"Using EasySearchParser with format_mode={config.format_mode} is not supported. Please use format_mode=4.")
        
        self.alignment_file = config.alignment_file
    
    def to_pandas(self) -> pd.DataFrame:
        """
        Returns a pandas DataFrame containing the alignment data.
        """
        return pd.read_csv(self.alignment_file, sep="\t")
    
    def to_list(self) -> list[dict]:
        """
        Returns a list of dictionaries containing the alignment data.
        """
        return self.to_pandas().to_dict(orient="records")
    
    def to_gen(self) -> Generator[dict, None, None]:
        """
        Returns a generator that yields dictionaries for each row in the alignment file.

        Each dictionary represents a row in the TSV file, with keys corresponding to 
        the column names in the header.
        """
        # Note that here the alignment_file is not a prefix, but a tsv file
        df_sample = pd.read_csv(self.alignment_file, sep="\t", nrows=1)
        column_types = df_sample.dtypes.to_dict()
        
        with open(self.alignment_file, 'r') as file:
            reader = csv.DictReader(file, delimiter='\t')
            for row in reader:
                for field, dtype in column_types.items():
                    if field in row:
                        try:
                            if dtype in (np.float64, float):
                                row[field] = float(row[field].replace('E', 'e'))
                            elif dtype in (np.int64, int):
                                row[field] = int(row[field])
                        except ValueError:
                            pass
                yield row
    
    def to_path(self) -> str:
        """
        Returns a list of file paths for the output files.

        Returns:
        --------
        list of str
        """
        return self.alignment_file
