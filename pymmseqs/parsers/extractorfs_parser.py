# pymmseqs/parsers/extractorfs_parser.py

import re
from typing import Generator, Optional

import pandas as pd

from ..config import ExtractOrfsConfig

# MMseqs2 ORF headers look like:  "<source_key>\t<start><strand><end>\t<frame>"
# e.g. "0\t0+149\t1"  (forward)  or  "0\t152-137\t1" (reverse)
_COORD_RE = re.compile(r"^(\d+)([+-])(\d+)$")


def _read_db_entries(path) -> list[str]:
    """Read an MMseqs2 DB data file into an ordered list of decoded entries.

    Entries are separated by NUL bytes; trailing whitespace/newlines are stripped.
    """
    with open(path, "rb") as f:
        raw = f.read()
    return [chunk.decode("utf-8", errors="replace").strip()
            for chunk in raw.split(b"\x00") if chunk.strip()]


class ExtractOrfsParser:
    """
    Parser for extractorfs output.

    extractorfs writes ORFs to a sequence database whose header DB encodes the
    coordinates of each ORF on its source contig. This parser merges those
    coordinates with the ORF sequences into a single table, and (best effort)
    resolves the numeric source key back to the original contig name.

    Parameters
    ----------
    config : ExtractOrfsConfig
        Configuration object after run() has been called.
    """

    def __init__(self, config: ExtractOrfsConfig):
        self.orf_db = str(config.orf_db)
        self.sequence_db = str(config.sequence_db)

    # ------------------------------------------------------------------ #
    # internal helpers
    # ------------------------------------------------------------------ #
    def _source_names(self) -> dict:
        """Map source DB key -> source contig accession (best effort).

        Uses the source header DB's `.index` (key->offset) so it stays correct
        even when source keys are non-contiguous (e.g. id_offset, merged DBs),
        rather than assuming keys run 0..N-1 in file order.
        """
        try:
            headers = _read_db_entries(f"{self.sequence_db}_h")
            with open(f"{self.sequence_db}_h.index") as f:
                rows = [line.split("\t") for line in f if line.strip()]
            # Header entries are stored in ascending-offset order; align keys the same way.
            keys_in_order = [int(r[0]) for r in sorted(rows, key=lambda r: int(r[1]))]
            return {
                key: (headers[i].split()[0] if headers[i] else "")
                for i, key in enumerate(keys_in_order)
                if i < len(headers)
            }
        except (FileNotFoundError, OSError, IndexError, ValueError):
            return {}

    def _records(self) -> list[dict]:
        headers = _read_db_entries(f"{self.orf_db}_h")
        sequences = _read_db_entries(self.orf_db)
        source_names = self._source_names()

        records = []
        for orf_id, header in enumerate(headers):
            fields = header.split("\t")
            source_id = int(fields[0]) if fields and fields[0].isdigit() else None
            start = end = strand = frame = None
            if len(fields) >= 2:
                m = _COORD_RE.match(fields[1])
                if m:
                    # MMseqs encodes the coordinate as start{strand}offset, where
                    # offset == length-1. Convert to an absolute end position on
                    # the source contig (end > start for '+', end < start for '-').
                    start, strand, offset = int(m.group(1)), m.group(2), int(m.group(3))
                    end = start + offset if strand == "+" else start - offset
            if len(fields) >= 3 and fields[2].lstrip("-").isdigit():
                frame = int(fields[2])

            sequence = sequences[orf_id] if orf_id < len(sequences) else ""

            source_name = source_names.get(source_id) if source_id is not None else None

            records.append({
                "orf_id": orf_id,
                "source_id": source_id,
                "source_name": source_name,
                "start": start,
                "end": end,
                "strand": strand,
                "frame": frame,
                "length": len(sequence),
                "sequence": sequence,
            })
        return records

    # ------------------------------------------------------------------ #
    # data access
    # ------------------------------------------------------------------ #
    def to_pandas(self) -> pd.DataFrame:
        """Return one row per ORF with coordinates, source contig, and sequence.

        `start` and `end` are absolute positions on the source contig. For a
        forward ORF (`strand == '+'`) `end > start`; for a reverse ORF
        (`strand == '-'`) `end < start` (the ORF reads from `start` down to `end`).
        In all cases `abs(end - start) + 1 == length`.
        """
        columns = ["orf_id", "source_id", "source_name", "start", "end",
                   "strand", "frame", "length", "sequence"]
        return pd.DataFrame(self._records(), columns=columns)

    def to_list(self) -> list[dict]:
        """Return all ORF records as a list of dictionaries."""
        return self._records()

    def to_gen(self) -> Generator[dict, None, None]:
        """Yield one ORF record dict at a time."""
        for record in self._records():
            yield record

    def to_path(self) -> str:
        """Return the path to the ORF sequence database (for chaining, e.g. convert2fasta/search)."""
        return self.orf_db

    def __len__(self) -> int:
        """Return the number of extracted ORFs."""
        return len(_read_db_entries(f"{self.orf_db}_h"))

    # ------------------------------------------------------------------ #
    # analysis
    # ------------------------------------------------------------------ #
    def summary(self) -> dict:
        """Compute summary statistics over the extracted ORFs.

        Returns
        -------
        dict
            total_orfs, source_sequences, forward_orfs, reverse_orfs,
            mean_length, median_length, min_length, max_length.
        """
        df = self.to_pandas()
        lengths = df["length"]
        return {
            "total_orfs": len(df),
            "source_sequences": int(df["source_id"].nunique()),
            "forward_orfs": int((df["strand"] == "+").sum()),
            "reverse_orfs": int((df["strand"] == "-").sum()),
            "mean_length": round(float(lengths.mean()), 2) if len(df) else 0,
            "median_length": int(lengths.median()) if len(df) else 0,
            "min_length": int(lengths.min()) if len(df) else 0,
            "max_length": int(lengths.max()) if len(df) else 0,
        }

    def __repr__(self) -> str:
        try:
            stats = self.summary()
        except Exception:
            return f"{self.__class__.__name__}(orf_db={self.orf_db!r})"
        lines = [f"{self.__class__.__name__}:"]
        for key, value in stats.items():
            lines.append(f"  {key}: {value}")
        return "\n".join(lines)
