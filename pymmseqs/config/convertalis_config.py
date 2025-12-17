from pathlib import Path
from typing import Union

from .base import BaseConfig
from ..defaults import loader
from ..utils import (
    get_caller_dir,
    run_mmseqs_command
)

DEFAULTS = loader.load("convertalis")

class ConvertAlisConfig(BaseConfig):
    """
    Convert alignment database to human-readable tabular format using MMseqs2 convertalis module.
    
    This class provides a Python interface to the MMseqs2 convertalis module,
    which converts alignment databases to readable formats like BLAST tabular format.
    """
    
    def __init__(
        self,
        # Required parameters
        query_db: Union[str, Path],
        target_db: Union[str, Path],
        alignment_db: Union[str, Path],
        alignment_file: Union[str, Path],
        
        # alignment parameters
        gap_open: str = "aa:11,nucl:5",
        gap_extend: str = "aa:1,nucl:2",

        # Misc parameters
        format_mode: int = 0,
        format_output: str = "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits",
        translation_table: int = 1,
        search_type: int = 0,

        # Common parameters
        sub_mat: str = "aa:blosum62.out,nucl:nucleotide.out",
        db_load_mode: int = 0,
        threads: Union[str, int] = 'all',
        compressed: bool = False,
        v: int = 3,

        # Expert parameters
        db_output: bool = False,
    ):
        """
        Convert alignment database to human-readable tabular format using MMseqs2 convertalis module.

        Parameters
        ----------
        `query_db` : Union[str, Path]
            Path to the query sequence database created with createdb
            
        `target_db` : Union[str, Path]
            Path to the target sequence database created with createdb
            
        `alignment_db` : Union[str, Path]
            Path to the alignment database to be converted (created by search or another module)
            
        `alignment_file` : Union[str, Path]
            Output path for the alignment file in human-readable format

        Alignment Parameters
        -------------------
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
        
        Misc Parameters
        --------------
        `format_mode` : int, optional  
            Output format type  
            - 0: BLAST-TAB (default)  
            - 1: SAM  
            - 2: BLAST-TAB + query/db length  
            - 3: Pretty HTML  
            - 4: BLAST-TAB + column headers  
            
        `format_output` : str, optional  
            Comma-separated list of output columns to include in results.  
            Available columns:  
            - query, target, evalue, gapopen, pident, fident, nident, qstart, qend, qlen  
            - tstart, tend, tlen, alnlen, raw, bits, cigar, qseq, tseq, qheader, theader, qaln, taln  
            - qframe, tframe, mismatch, qcov, tcov, qset, qsetid, tset, tsetid, taxid, taxname, taxlineage  
            - qorfstart, qorfend, torfstart, torfend, ppos  

            - Default: "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"
        
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
        
        `search_type` : int, optional
            Search mode:
            - 0: auto (default)
            - 1: amino
            - 2: translated
            - 3: nucleotide
            - 4: translated alignment

        Common Parameters
        ----------------
        `sub_mat` : str, optional
            Substitution matrix (type:path, type:path)
            type: "aa" or "nucl"
            path: matrix file path
            - "aa:blosum62.out,nucl:nucleotide.out" (default)

            Note: find available matrices in the MMseqs2 data directory: (https://github.com/soedinglab/MMseqs2/tree/master/data)
            
        `db_load_mode` : int, optional
            Database preloading method
            - 0: auto (default)
            - 1: fread
            - 2: mmap
            - 3: mmap+touch
            
        `threads` : int, optional
            Number of CPU-cores used
            - 14 (default)
            
        `compressed` : bool, optional
            Write compressed output
            - True
            - False (default)
            
        `v` : int, optional
            Output verbosity
            - 0: quiet
            - 1: +errors
            - 2: +warnings
            - 3: +info (default)

        Expert Parameters
        ----------------
        `db_output` : bool, optional
            Return a database instead of a flat file
            - True
            - False (default)

        Returns
        -------
        None
            Creates alignment file at specified output path

        Raises
        ------
        FileNotFoundError
            If input databases are missing

        ValueError
            For invalid parameter combinations/values

        Examples
        --------
        Basic alignment conversion to BLAST tabular format:
        >>> config = ConvertAlisConfig(
            query_db="query.db",
            target_db="target.db",
            alignment_db="results.db",
            alignment_file="results.m8",
            threads=8
        )
        >>> config.run()
        """
        super().__init__()
        
        # Required parameters
        self.query_db = Path(query_db)
        self.target_db = Path(target_db)
        self.alignment_db = Path(alignment_db)
        self.alignment_file = Path(alignment_file)
        
        # Alignment parameters
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        
        # Misc parameters
        self.format_mode = format_mode
        self.format_output = format_output
        self.translation_table = translation_table
        self.search_type = search_type
        
        # Common parameters
        self.sub_mat = sub_mat
        self.db_load_mode = db_load_mode
        self.threads = threads
        self.compressed = compressed
        self.v = v
        
        # Expert parameters
        self.db_output = db_output
        
        self._defaults = DEFAULTS
        self._path_params = [param for param, info in DEFAULTS.items() if info['type'] == 'path']
        self._caller_dir = get_caller_dir()
        
    def _validate(self) -> None:
        self._check_required_files()
        self._validate_choices()
        
        # Additional validations
        if not self.threads >= 0:
            raise ValueError(f"Number of threads must be >= 0, got {self.threads}")
        
    
    def run(self) -> None:
        self._resolve_all_path(self._caller_dir)
        
        self._validate()
        
        args = self._get_command_args("convertalis")
        mmseqs_output = run_mmseqs_command(args)
        
        self._handle_command_output(
            mmseqs_output=mmseqs_output,
            output_identifier="ConvertAlis",
            output_path=str(self.alignment_file)
        )
