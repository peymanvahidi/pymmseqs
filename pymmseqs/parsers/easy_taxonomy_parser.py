# pymmseqs/parsers/easy_taxonomy_parser.py

import pandas as pd
from pathlib import Path
from typing import Generator, Union

from ..config import EasyTaxonomyConfig


# MMseqs2 easy-taxonomy report uses full rank names (e.g., 'phylum', 'class')
# as rank_code values, not Kraken-style single-letter codes.
# Normalize common synonyms and alternative names to canonical forms.
RANK_ALIASES = {
    'superkingdom': 'domain',
}

RANK_ORDER = [
    'unclassified', 'root', 'domain', 'superkingdom', 'kingdom',
    'phylum', 'subphylum', 'superclass', 'class', 'subclass',
    'superorder', 'order', 'suborder', 'infraorder', 'parvorder',
    'superfamily', 'family', 'subfamily', 'genus', 'subgenus',
    'species group', 'species', 'subspecies', 'strain',
]


class EasyTaxonomyParser:
    """Parser for easy-taxonomy output with data access, analysis, and visualization.

    Provides methods for accessing taxonomy results as Python objects,
    computing taxonomic composition and diversity metrics, and generating
    visualizations of classification results.

    Parameters
    ----------
    config : EasyTaxonomyConfig
        Configuration object after run() has been called.
    """

    def __init__(self, config: EasyTaxonomyConfig):
        self._prefix = str(config.tax_reports)
        self._format_output = config.format_output
        self._format_mode = config.format_mode
        self._tax_lineage = config.tax_lineage
        self._lca_ranks = config.lca_ranks

    # ------------------------------------------------------------------ #
    #  Internal helpers                                                    #
    # ------------------------------------------------------------------ #

    def _parse_report_file(self, filepath: str) -> pd.DataFrame:
        """Parse a Kraken-style report file preserving hierarchy depth.

        Kraken report columns (tab-separated, no header):
        percentage, num_reads_clade, num_reads_direct, rank_code, taxid, name
        """
        rows = []
        with open(filepath) as f:
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 6:
                    name_raw = parts[5]
                    depth = (len(name_raw) - len(name_raw.lstrip(' '))) // 2
                    rows.append({
                        'percentage': float(parts[0]),
                        'num_reads_clade': int(parts[1]),
                        'num_reads_direct': int(parts[2]),
                        'rank_code': parts[3].strip(),
                        'taxid': int(parts[4]),
                        'name': name_raw.strip(),
                        'depth': depth,
                    })
        return pd.DataFrame(rows)

    def _get_lca_columns(self, num_cols: int) -> list[str]:
        """Determine LCA file column names based on actual column count and config."""
        # Default LCA format: query, taxid, rank, taxon_name (4 columns)
        # With --tax-lineage: query, taxid, rank, taxon_name, lineage
        # With --lca-ranks: additional rank columns appended
        base = ['query', 'taxid', 'rank', 'taxon_name']
        columns = base[:min(num_cols, len(base))]
        if num_cols > len(base):
            if self._tax_lineage in (1, 2):
                columns.append('lineage')
            if self._lca_ranks:
                for rank in self._lca_ranks.split(','):
                    rank = rank.strip()
                    if rank:
                        columns.append(rank)
        # Pad with generic names if still short
        while len(columns) < num_cols:
            columns.append(f'col_{len(columns)}')
        return columns

    def _get_classified_counts(self) -> tuple[int, int, int]:
        """Return (total, classified, unclassified) from report."""
        report_df = self.report()

        # MMseqs2 report: unclassified has taxid=0, root has taxid=1
        unclassified_row = report_df[report_df['taxid'] == 0]
        root_row = report_df[report_df['taxid'] == 1]

        unclassified = int(unclassified_row['num_reads_clade'].iloc[0]) if not unclassified_row.empty else 0
        classified = int(root_row['num_reads_clade'].iloc[0]) if not root_row.empty else 0
        total = classified + unclassified
        return total, classified, unclassified

    # ------------------------------------------------------------------ #
    #  Data Access Methods (S1-S8)                                         #
    # ------------------------------------------------------------------ #

    def to_pandas(self, output: str = 'lca') -> pd.DataFrame:
        """Return results as a pandas DataFrame.

        Parameters
        ----------
        output : str
            Which output to parse:
            - 'lca': LCA assignments (query -> taxid)
            - 'report': Kraken-style taxonomic report
            - 'tophit_aln': Top-hit alignments
            - 'tophit_report': Top-hit coverage report

        Returns
        -------
        pd.DataFrame
        """
        if output == 'lca':
            return self.lca_assignments()
        elif output == 'report':
            return self.report()
        elif output == 'tophit_aln':
            return self.top_hits()
        elif output == 'tophit_report':
            return self._parse_report_file(f"{self._prefix}_tophit_report")
        else:
            valid = "'lca', 'report', 'tophit_aln', 'tophit_report'"
            raise ValueError(f"Unknown output '{output}'. Choose from: {valid}")

    def to_list(self, output: str = 'lca') -> list[dict]:
        """Return results as a list of dictionaries.

        Parameters
        ----------
        output : str
            Which output to parse (see to_pandas for options).
        """
        return self.to_pandas(output=output).to_dict(orient='records')

    def to_gen(self) -> Generator[dict, None, None]:
        """Yield LCA assignments one row at a time.

        Yields
        ------
        dict
            Dictionary with keys matching LCA file columns.
        """
        import csv

        lca_file = f"{self._prefix}_lca.tsv"
        with open(lca_file) as f:
            reader = csv.reader(f, delimiter='\t')
            first_row = next(reader, None)
            if first_row is None:
                return
            columns = self._get_lca_columns(len(first_row))

            def _make_entry(row):
                entry = {}
                for i, val in enumerate(row):
                    key = columns[i] if i < len(columns) else f'col_{i}'
                    if key == 'taxid':
                        try:
                            val = int(val)
                        except ValueError:
                            pass
                    entry[key] = val
                return entry

            yield _make_entry(first_row)
            for row in reader:
                yield _make_entry(row)

    def to_path(self) -> dict[str, str]:
        """Return output file paths.

        Returns
        -------
        dict[str, str]
            Mapping of output name to file path.
        """
        return {
            'lca': f"{self._prefix}_lca.tsv",
            'report': f"{self._prefix}_report",
            'tophit_aln': f"{self._prefix}_tophit_aln",
            'tophit_report': f"{self._prefix}_tophit_report",
        }

    def to_json(self, path: Union[str, Path, None] = None) -> Union[str, None]:
        """Export LCA results as JSON.

        Parameters
        ----------
        path : str or Path, optional
            File path to write JSON. If None, returns JSON string.
        """
        import json

        data = self.to_list()
        if path is not None:
            with open(path, 'w') as f:
                json.dump(data, f, indent=2)
            return None
        return json.dumps(data, indent=2)

    def to_csv(self, path: Union[str, Path], output: str = 'lca') -> None:
        """Export results as CSV file.

        Parameters
        ----------
        path : str or Path
            Output CSV file path.
        output : str
            Which output to export (see to_pandas for options).
        """
        self.to_pandas(output=output).to_csv(path, index=False)

    def __len__(self) -> int:
        """Return number of query sequences in LCA results."""
        lca_file = f"{self._prefix}_lca.tsv"
        with open(lca_file) as f:
            return sum(1 for _ in f)

    def __repr__(self) -> str:
        """Pretty summary for notebooks and REPL."""
        try:
            stats = self.summary()
            lines = ['EasyTaxonomyParser:']
            for key, value in stats.items():
                if isinstance(value, list):
                    value = ', '.join(str(v) for v in value)
                lines.append(f'  {key}: {value}')
            return '\n'.join(lines)
        except Exception:
            return f'EasyTaxonomyParser(prefix="{self._prefix}")'

    # ------------------------------------------------------------------ #
    #  Core Analysis Methods (C1, C3-C9, X1, X3)                          #
    # ------------------------------------------------------------------ #

    def report(self) -> pd.DataFrame:
        """Parse the Kraken-style taxonomic report.

        Returns
        -------
        pd.DataFrame
            Columns: percentage, num_reads_clade, num_reads_direct,
            rank_code, taxid, name, depth.
        """
        return self._parse_report_file(f"{self._prefix}_report")

    def lca_assignments(self) -> pd.DataFrame:
        """Parse LCA assignment file.

        Returns
        -------
        pd.DataFrame
            Columns depend on config (always has query, taxid; may
            include lineage or per-rank columns).
        """
        lca_file = f"{self._prefix}_lca.tsv"
        df = pd.read_csv(lca_file, sep='\t', header=None)
        columns = self._get_lca_columns(len(df.columns))
        df.columns = columns[:len(df.columns)]
        if 'taxid' in df.columns:
            df['taxid'] = pd.to_numeric(df['taxid'], errors='coerce').fillna(0).astype(int)
        return df

    def top_hits(self) -> pd.DataFrame:
        """Parse top-hit alignment file.

        Returns
        -------
        pd.DataFrame
            Alignment data with columns matching format_output setting.
        """
        tophit_file = f"{self._prefix}_tophit_aln"
        if self._format_mode == 4:
            return pd.read_csv(tophit_file, sep='\t')
        else:
            columns = [c.strip() for c in self._format_output.split(',')]
            return pd.read_csv(tophit_file, sep='\t', header=None, names=columns)

    def summary(self) -> dict:
        """Compute overall taxonomy classification summary.

        Returns
        -------
        dict
            Keys: total_queries, classified, unclassified,
            classified_pct, num_taxa, top_phyla.
        """
        total, classified, unclassified = self._get_classified_counts()
        report_df = self.report()

        try:
            comp = self.composition('phylum')
            top_phyla = comp.head(5)['taxon'].tolist()
        except (ValueError, KeyError):
            top_phyla = []

        taxa_count = len(report_df[~report_df['rank_code'].isin(['U', 'R'])])

        return {
            'total_queries': total,
            'classified': classified,
            'unclassified': unclassified,
            'classified_pct': round(classified / total * 100, 2) if total > 0 else 0.0,
            'num_taxa': taxa_count,
            'top_phyla': top_phyla,
        }

    def composition(self, rank: str = 'phylum') -> pd.DataFrame:
        """Get taxonomic composition at a given rank.

        Parameters
        ----------
        rank : str
            Taxonomic rank: 'domain', 'superkingdom', 'kingdom', 'phylum',
            'class', 'order', 'family', 'genus', 'species', etc.

        Returns
        -------
        pd.DataFrame
            Columns: taxon, count, proportion. Sorted by count descending.
        """
        report_df = self.report()
        rank_lower = rank.lower()

        # Resolve aliases (e.g., 'superkingdom' -> also match 'domain')
        search_ranks = {rank_lower}
        for alias, canonical in RANK_ALIASES.items():
            if rank_lower == alias:
                search_ranks.add(canonical)
            elif rank_lower == canonical:
                search_ranks.add(alias)

        # MMseqs2 report uses full rank names as rank_code
        filtered = report_df[report_df['rank_code'].str.lower().isin(search_ranks)].copy()

        if filtered.empty:
            return pd.DataFrame(columns=['taxon', 'count', 'proportion'])

        total = filtered['num_reads_clade'].sum()
        result = filtered[['name', 'num_reads_clade']].copy()
        result.columns = ['taxon', 'count']
        result['proportion'] = result['count'] / total if total > 0 else 0.0
        return result.sort_values('count', ascending=False).reset_index(drop=True)

    def diversity(self, rank: str = 'phylum') -> dict:
        """Compute alpha diversity metrics at a given rank.

        Parameters
        ----------
        rank : str
            Taxonomic rank to compute diversity at.

        Returns
        -------
        dict
            Keys: richness, shannon_entropy, simpson_index, evenness.
        """
        import numpy as np

        comp = self.composition(rank=rank)
        if comp.empty:
            return {'richness': 0, 'shannon_entropy': 0.0, 'simpson_index': 0.0, 'evenness': 0.0}

        proportions = comp['proportion'].values
        proportions = proportions[proportions > 0]

        richness = len(proportions)
        shannon = float(-np.sum(proportions * np.log(proportions)))
        simpson = float(1 - np.sum(proportions ** 2))
        evenness = shannon / np.log(richness) if richness > 1 else 0.0

        return {
            'richness': richness,
            'shannon_entropy': round(shannon, 4),
            'simpson_index': round(simpson, 4),
            'evenness': round(evenness, 4),
        }

    def rank_summary(self) -> pd.DataFrame:
        """Get classification count at each taxonomic level.

        Returns
        -------
        pd.DataFrame
            Columns: rank, rank_code, num_taxa, num_reads.
            Sorted by canonical rank order.
        """
        report_df = self.report()
        summary = report_df.groupby('rank_code').agg(
            num_taxa=('taxid', 'count'),
            num_reads=('num_reads_direct', 'sum')
        ).reset_index()

        # rank_code already contains the rank name in MMseqs2 reports
        summary['rank'] = summary['rank_code']

        def rank_sort_key(r):
            try:
                return RANK_ORDER.index(r.lower())
            except ValueError:
                return 99

        summary['_order'] = summary['rank'].apply(rank_sort_key)
        return summary.sort_values('_order').drop(columns='_order').reset_index(drop=True)[
            ['rank', 'rank_code', 'num_taxa', 'num_reads']
        ]

    def unclassified_report(self) -> dict:
        """Get classification rate statistics.

        Returns
        -------
        dict
            Keys: total_queries, classified, unclassified,
            classified_pct, unclassified_pct.
        """
        total, classified, unclassified = self._get_classified_counts()
        return {
            'total_queries': total,
            'classified': classified,
            'unclassified': unclassified,
            'classified_pct': round(classified / total * 100, 2) if total > 0 else 0.0,
            'unclassified_pct': round(unclassified / total * 100, 2) if total > 0 else 0.0,
        }

    def filter_by_taxon(self, name: str, rank: str = 'phylum') -> pd.DataFrame:
        """Filter LCA assignments to a specific taxonomic lineage.

        Returns all queries classified at or below the specified taxon.

        Parameters
        ----------
        name : str
            Taxon name (e.g., 'Proteobacteria').
        rank : str
            Taxonomic rank of the taxon.

        Returns
        -------
        pd.DataFrame
            Filtered LCA assignments.
        """
        report_df = self.report()
        rank_lower = rank.lower()

        # Resolve aliases
        search_ranks = {rank_lower}
        for alias, canonical in RANK_ALIASES.items():
            if rank_lower == alias:
                search_ranks.add(canonical)
            elif rank_lower == canonical:
                search_ranks.add(alias)

        target = report_df[
            (report_df['name'] == name) & (report_df['rank_code'].str.lower().isin(search_ranks))
        ]
        if target.empty:
            raise ValueError(f"Taxon '{name}' at rank '{rank}' not found in report.")

        target_idx = target.index[0]
        target_depth = target.iloc[0]['depth']

        descendant_taxids = set()
        for i in range(target_idx, len(report_df)):
            row = report_df.iloc[i]
            if i > target_idx and row['depth'] <= target_depth:
                break
            descendant_taxids.add(row['taxid'])

        lca_df = self.lca_assignments()
        return lca_df[lca_df['taxid'].isin(descendant_taxids)].reset_index(drop=True)

    # ------------------------------------------------------------------ #
    #  Visualization Methods (C2, C7, X2, X4, X5)                         #
    # ------------------------------------------------------------------ #

    def plot_composition(self, rank: str = 'phylum', top_n: int = 10,
                         figsize: tuple = (10, 6), ax=None):
        """Plot horizontal bar chart of taxonomic composition.

        Parameters
        ----------
        rank : str
            Taxonomic rank to plot.
        top_n : int
            Number of top taxa to show.
        figsize : tuple
            Figure size as (width, height).
        ax : matplotlib.axes.Axes, optional
            Axes to plot on. If None, creates a new figure.

        Returns
        -------
        tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
        """
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.get_figure()

        comp = self.composition(rank=rank)
        if comp.empty:
            ax.text(0.5, 0.5, f'No {rank}-level classifications',
                    ha='center', va='center', transform=ax.transAxes)
            fig.tight_layout()
            return fig, ax

        plot_data = comp.head(top_n)
        ax.barh(range(len(plot_data)), plot_data['count'],
                color='steelblue', edgecolor='black', alpha=0.7)
        ax.set_yticks(range(len(plot_data)))
        ax.set_yticklabels(plot_data['taxon'])
        ax.invert_yaxis()
        ax.set_xlabel('Number of Reads')
        ax.set_title(f'Top {min(top_n, len(plot_data))} Taxa at {rank.capitalize()} Level')

        fig.tight_layout()
        return fig, ax

    def plot_rank_resolution(self, figsize: tuple = (10, 6), ax=None):
        """Plot bar chart showing reads assigned at each taxonomic rank.

        Parameters
        ----------
        figsize : tuple
            Figure size as (width, height).
        ax : matplotlib.axes.Axes, optional
            Axes to plot on.

        Returns
        -------
        tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
        """
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.get_figure()

        rank_df = self.rank_summary()
        # Exclude non-taxonomic ranks (unclassified/root use 'no rank' or 'cellular root')
        rank_df = rank_df[~rank_df['rank'].isin(['no rank', 'cellular root'])]

        if rank_df.empty:
            ax.text(0.5, 0.5, 'No rank data available',
                    ha='center', va='center', transform=ax.transAxes)
            fig.tight_layout()
            return fig, ax

        ax.bar(range(len(rank_df)), rank_df['num_reads'],
               color='teal', edgecolor='black', alpha=0.7)
        ax.set_xticks(range(len(rank_df)))
        ax.set_xticklabels(rank_df['rank'].str.capitalize(), rotation=45, ha='right')
        ax.set_ylabel('Number of Reads')
        ax.set_title('Classification Depth: Reads Assigned per Rank')

        fig.tight_layout()
        return fig, ax

    def plot_composition_pie(self, rank: str = 'phylum', top_n: int = 8,
                             figsize: tuple = (8, 8), ax=None):
        """Plot pie chart of taxonomic composition with 'Other' slice.

        Parameters
        ----------
        rank : str
            Taxonomic rank to plot.
        top_n : int
            Number of top taxa to show individually.
        figsize : tuple
            Figure size as (width, height).
        ax : matplotlib.axes.Axes, optional
            Axes to plot on.

        Returns
        -------
        tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
        """
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.get_figure()

        comp = self.composition(rank=rank)
        if comp.empty:
            ax.text(0.5, 0.5, f'No {rank}-level classifications',
                    ha='center', va='center', transform=ax.transAxes)
            fig.tight_layout()
            return fig, ax

        top = comp.head(top_n)
        other_count = comp.iloc[top_n:]['count'].sum() if len(comp) > top_n else 0

        labels = list(top['taxon'])
        sizes = list(top['count'])
        if other_count > 0:
            labels.append('Other')
            sizes.append(other_count)

        ax.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90)
        ax.set_title(f'Taxonomic Composition at {rank.capitalize()} Level')

        fig.tight_layout()
        return fig, ax

    def plot_classified_vs_unclassified(self, figsize: tuple = (6, 6), ax=None):
        """Plot pie chart of classified vs unclassified queries.

        Parameters
        ----------
        figsize : tuple
            Figure size as (width, height).
        ax : matplotlib.axes.Axes, optional
            Axes to plot on.

        Returns
        -------
        tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
        """
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.get_figure()

        stats = self.unclassified_report()
        sizes = [stats['classified'], stats['unclassified']]
        labels = [
            f"Classified ({stats['classified_pct']}%)",
            f"Unclassified ({stats['unclassified_pct']}%)",
        ]
        colors = ['steelblue', 'lightcoral']

        ax.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
        ax.set_title('Classification Rate')

        fig.tight_layout()
        return fig, ax

    def plot_diversity_comparison(self, ranks: list[str] = None,
                                 figsize: tuple = (10, 6), ax=None):
        """Plot diversity metrics across taxonomic ranks.

        Parameters
        ----------
        ranks : list[str], optional
            Ranks to compare. Defaults to phylum through species.
        figsize : tuple
            Figure size as (width, height).
        ax : matplotlib.axes.Axes, optional
            Axes to plot on.

        Returns
        -------
        tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
        """
        import matplotlib.pyplot as plt
        import numpy as np

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.get_figure()

        if ranks is None:
            ranks = ['phylum', 'class', 'order', 'family', 'genus', 'species']

        rank_labels = []
        shannon_values = []
        simpson_values = []
        for rank in ranks:
            try:
                div = self.diversity(rank=rank)
                rank_labels.append(rank.capitalize())
                shannon_values.append(div['shannon_entropy'])
                simpson_values.append(div['simpson_index'])
            except (ValueError, KeyError):
                continue

        if not rank_labels:
            ax.text(0.5, 0.5, 'No diversity data available',
                    ha='center', va='center', transform=ax.transAxes)
            fig.tight_layout()
            return fig, ax

        x = np.arange(len(rank_labels))
        width = 0.35

        ax.bar(x - width / 2, shannon_values, width,
               label='Shannon Entropy', color='steelblue', alpha=0.7, edgecolor='black')
        ax.bar(x + width / 2, simpson_values, width,
               label="Simpson's Index", color='coral', alpha=0.7, edgecolor='black')

        ax.set_xticks(x)
        ax.set_xticklabels(rank_labels, rotation=45, ha='right')
        ax.set_ylabel('Diversity Index')
        ax.set_title('Diversity Metrics Across Taxonomic Ranks')
        ax.legend()

        fig.tight_layout()
        return fig, ax
