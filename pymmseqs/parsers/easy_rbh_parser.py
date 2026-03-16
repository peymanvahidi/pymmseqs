# pymmseqs/parsers/easy_rbh_parser.py

import csv
from pathlib import Path
from typing import Generator, Union

import numpy as np
import pandas as pd

from ..config import EasyRbhConfig


class EasyRbhParser:
    """Parser for easy-rbh reciprocal best hit output.

    Provides methods for accessing ortholog pairs as Python objects,
    computing summary statistics, and generating visualizations.

    Supports both format_mode=0 (headerless BLAST-TAB, uses format_output
    for column names) and format_mode=4 (BLAST-TAB with column headers).

    Parameters
    ----------
    config : EasyRbhConfig
        Configuration object after run() has been called.
    """

    def __init__(self, config: EasyRbhConfig):
        self._alignment_file = config.alignment_file
        self._format_mode = config.format_mode
        self._format_output = config.format_output

    def _get_columns(self) -> list[str]:
        """Derive column names from format_output string."""
        return [col.strip() for col in self._format_output.split(",")]

    def to_pandas(self) -> pd.DataFrame:
        """Return reciprocal best hits as a pandas DataFrame.

        Returns
        -------
        pd.DataFrame
            DataFrame with one row per ortholog pair.
        """
        if self._format_mode == 4:
            return pd.read_csv(self._alignment_file, sep="\t")
        return pd.read_csv(
            self._alignment_file, sep="\t",
            header=None, names=self._get_columns()
        )

    def to_list(self) -> list[dict]:
        """Return reciprocal best hits as a list of dictionaries.

        Returns
        -------
        list[dict]
            Each dict represents one ortholog pair.
        """
        return self.to_pandas().to_dict(orient="records")

    def to_gen(self) -> Generator[dict, None, None]:
        """Yield reciprocal best hits one row at a time with type conversion.

        Yields
        ------
        dict
            One ortholog pair per iteration.
        """
        if self._format_mode == 4:
            df_sample = pd.read_csv(self._alignment_file, sep="\t", nrows=1)
            column_types = df_sample.dtypes.to_dict()

            with open(self._alignment_file, "r") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    for field, dtype in column_types.items():
                        if field in row:
                            try:
                                if dtype in (np.float64, float):
                                    row[field] = float(row[field].replace("E", "e"))
                                elif dtype in (np.int64, int):
                                    row[field] = int(row[field])
                            except (ValueError, AttributeError):
                                pass
                    yield row
        else:
            columns = self._get_columns()
            df_sample = pd.read_csv(
                self._alignment_file, sep="\t",
                header=None, names=columns, nrows=1
            )
            column_types = df_sample.dtypes.to_dict()

            with open(self._alignment_file, "r") as f:
                reader = csv.reader(f, delimiter="\t")
                for fields in reader:
                    row = dict(zip(columns, fields))
                    for col, dtype in column_types.items():
                        if col in row:
                            try:
                                if dtype in (np.float64, float):
                                    row[col] = float(row[col].replace("E", "e"))
                                elif dtype in (np.int64, int):
                                    row[col] = int(row[col])
                            except (ValueError, AttributeError):
                                pass
                    yield row

    def to_path(self) -> str:
        """Return the output alignment file path.

        Returns
        -------
        str
            Path to the reciprocal best hit alignment file.
        """
        return str(self._alignment_file)

    def to_json(self, path: Union[str, Path, None] = None) -> Union[str, None]:
        """Export reciprocal best hits as JSON.

        Parameters
        ----------
        path : str or Path, optional
            File path to write JSON. If None, returns JSON string.

        Returns
        -------
        str or None
            JSON string if path is None, otherwise None.
        """
        import json

        data = self.to_list()
        if path is not None:
            with open(path, "w") as f:
                json.dump(data, f, indent=2)
            return None
        return json.dumps(data, indent=2)

    def to_csv(self, path: Union[str, Path]) -> None:
        """Export reciprocal best hits as CSV.

        Parameters
        ----------
        path : str or Path
            Output CSV file path.
        """
        self.to_pandas().to_csv(path, index=False)

    def summary(self) -> dict:
        """Compute summary statistics for reciprocal best hits.

        Returns
        -------
        dict
            Always contains 'total_pairs'. Additionally includes
            'mean_identity' and 'median_identity' (if fident column present),
            'median_evalue' (if evalue column present),
            'mean_alignment_length' (if alnlen column present), and
            'median_bit_score' (if bits column present).
            Available keys depend on format_output columns.
        """
        df = self.to_pandas()
        cols = df.columns.tolist()
        result = {"total_pairs": len(df)}

        if "fident" in cols:
            result["mean_identity"] = round(df["fident"].mean(), 4)
            result["median_identity"] = round(df["fident"].median(), 4)
        if "evalue" in cols:
            result["median_evalue"] = df["evalue"].median()
        if "alnlen" in cols:
            result["mean_alignment_length"] = round(df["alnlen"].mean(), 2)
        if "bits" in cols:
            result["median_bit_score"] = round(df["bits"].median(), 2)

        return result

    def plot_identity_distribution(self, bins=50, figsize=(10, 6), ax=None):
        """Plot histogram of sequence identity distribution across ortholog pairs.

        Parameters
        ----------
        bins : int, optional
            Number of histogram bins. Default is 50.
        figsize : tuple, optional
            Figure size as (width, height). Default is (10, 6).
        ax : matplotlib.axes.Axes, optional
            Axes to plot on. If None, creates a new figure.

        Returns
        -------
        tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
            The figure and axes objects for further customization.
        """
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.get_figure()

        df = self.to_pandas()
        if "fident" not in df.columns and "pident" not in df.columns:
            raise ValueError(
                f"Neither 'fident' nor 'pident' column found. Available: {df.columns.tolist()}"
            )
        col = "fident" if "fident" in df.columns else "pident"
        ax.hist(df[col], bins=bins, edgecolor="black", alpha=0.7, color="steelblue")
        ax.set_xlabel("Sequence Identity")
        ax.set_ylabel("Count")
        ax.set_title("RBH Sequence Identity Distribution")
        ax.axvline(
            df[col].mean(), color="red", linestyle="--",
            label=f"Mean: {df[col].mean():.3f}"
        )
        ax.legend()

        fig.tight_layout()
        return fig, ax

    def plot_evalue_distribution(self, bins=50, figsize=(10, 6), ax=None):
        """Plot log-scale histogram of e-value distribution.

        Parameters
        ----------
        bins : int, optional
            Number of histogram bins. Default is 50.
        figsize : tuple, optional
            Figure size as (width, height). Default is (10, 6).
        ax : matplotlib.axes.Axes, optional
            Axes to plot on. If None, creates a new figure.

        Returns
        -------
        tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
            The figure and axes objects for further customization.
        """
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.get_figure()

        df = self.to_pandas()
        if "evalue" not in df.columns:
            raise ValueError(f"Column 'evalue' not found. Available: {df.columns.tolist()}")
        evalues = df["evalue"][df["evalue"] > 0]
        ax.hist(np.log10(evalues), bins=bins, edgecolor="black", alpha=0.7, color="coral")
        ax.set_xlabel("log10(E-value)")
        ax.set_ylabel("Count")
        ax.set_title("RBH E-value Distribution")

        fig.tight_layout()
        return fig, ax

    def plot_alignment_length(self, bins=50, figsize=(10, 6), ax=None):
        """Plot histogram of alignment lengths.

        Parameters
        ----------
        bins : int, optional
            Number of histogram bins. Default is 50.
        figsize : tuple, optional
            Figure size as (width, height). Default is (10, 6).
        ax : matplotlib.axes.Axes, optional
            Axes to plot on. If None, creates a new figure.

        Returns
        -------
        tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
            The figure and axes objects for further customization.
        """
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.get_figure()

        df = self.to_pandas()
        if "alnlen" not in df.columns:
            raise ValueError(f"Column 'alnlen' not found. Available: {df.columns.tolist()}")
        ax.hist(df["alnlen"], bins=bins, edgecolor="black", alpha=0.7, color="mediumpurple")
        ax.set_xlabel("Alignment Length")
        ax.set_ylabel("Count")
        ax.set_title("RBH Alignment Length Distribution")
        ax.axvline(
            df["alnlen"].mean(), color="red", linestyle="--",
            label=f"Mean: {df['alnlen'].mean():.1f}"
        )
        ax.legend()

        fig.tight_layout()
        return fig, ax

    def plot_score_distribution(self, score="bits", bins=50, figsize=(10, 6), ax=None):
        """Plot histogram of any numeric score column.

        Parameters
        ----------
        score : str, optional
            Column name to plot. Default is 'bits'.
        bins : int, optional
            Number of histogram bins. Default is 50.
        figsize : tuple, optional
            Figure size as (width, height). Default is (10, 6).
        ax : matplotlib.axes.Axes, optional
            Axes to plot on. If None, creates a new figure.

        Returns
        -------
        tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
            The figure and axes objects for further customization.
        """
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.get_figure()

        df = self.to_pandas()
        if score not in df.columns:
            raise ValueError(f"Column '{score}' not found. Available: {df.columns.tolist()}")

        ax.hist(df[score], bins=bins, edgecolor="black", alpha=0.7, color="seagreen")
        ax.set_xlabel(score)
        ax.set_ylabel("Count")
        ax.set_title(f"RBH {score} Distribution")

        fig.tight_layout()
        return fig, ax
