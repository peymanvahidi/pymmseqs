# pymmseqs.parsers
This module is for parsing the output of MMseqs2 commands as Python objects, making it easier to work with the results in your Python scripts.

# Important Things to Know

- You obtain a parser object when using the `pymmseqs.commands` module to execute a command.
- Advanced users can pass the `pymmseqs.config` object to the parser after execution to leverage additional features of `pymmseqs.parsers`.
- Refer to the sections below for detailed instructions on how to use parser objects effectively.

---

# [CreateDBParser](https://github.com/heispv/pymmseqs/blob/master/pymmseqs/parsers/createdb_parser.py)
This parser processes the output of the MMseqs2 `createdb` command.

## For Basic Users
When using the `createdb` from `pymmseqs.commands`, you receive a `CreateDBParser` object.

Example:
```python
from pymmseqs.commands import createdb

query_database = createdb(
    fasta_file="data/query.fasta",
    sequence_db="output/query_db",
)

# Retrieve the path to the database
print(query_database.get_path())
```

## For Advanced Users
Advanced users can utilize the `pymmseqs.config.CreateDBConfig` object for additional flexibility.

Example:
```python
from pymmseqs.config import CreateDBConfig
from pymmseqs.parsers import CreateDBParser

# Create the configuration
createdb_config = CreateDBConfig(
    fasta_file="data/query.fasta",
    sequence_db="output/query_db",
    dbtype=1,
)

# Execute the configuration
createdb_config.run()

# Obtain the parser object
createdb_parser = CreateDBParser(createdb_config)

# Retrieve the path to the database
print(createdb_parser.get_path())
```

**Want to know what happens when you run `config.run()`?**

Check out the [`pymmseqs.config`](../pymmseqs.config) module for details.

---

# [EasyClusterParser](https://github.com/heispv/pymmseqs/blob/master/pymmseqs/parsers/easy_cluster_parser.py)
The `EasyClusterParser` provides several methods for handling clustering results:

## Methods:

### `to_list()`
- Returns a list of dictionaries, each representing a cluster.
- Each dictionary contains:
    - **`rep`**: The representative sequence ID.
    - **`members`**: A list of dictionaries for cluster members with:
        - `seq_id`: Unique sequence identifier.
        - `header`: Full FASTA header.
        - `sequence`: Sequence data (nucleotide or protein).

**Example Output:**
```python
[
    {
        "rep": "seq1",
        "members": [
            {"seq_id": "seq1", "header": "sp|seq1", "sequence": "SEQWENCE"},
            {"seq_id": "seq2", "header": "sp|seq2", "sequence": "PRTEINSEQWENCE"},
        ]
    },
    {
        "rep": "seq3",
        "members": [
            {"seq_id": "seq3", "header": "sp|seq3", "sequence": "PRTEIN"}
        ]
    }
]
```

### `to_pandas()`
- Converts cluster data into a **pandas DataFrame** for easier analysis.
- Columns include:
    - `rep`: Representative sequence ID (index)
    - `seq_id`: Sequence identifier
    - `header`: FASTA header
    - `sequence`: Sequence data

### `to_gen()`
- Returns a **generator** of clusters, allowing memory-efficient iteration.

### `to_rep_list()`
- Returns a list of representative sequences from all clusters.
- **Parameters:**
    - `with_seq` (bool): If True, returns tuples of (header, sequence). If False, returns only headers.
- **Returns:** List of representative sequences or tuples.

**Example:**
```python
# Get representatives with sequences
reps_with_seq = cluster_result.to_rep_list(with_seq=True)
# Returns: [("seq1", "SEQWENCE"), ("seq3", "PRTEIN")]

# Get only representative headers
rep_headers = cluster_result.to_rep_list(with_seq=False)
# Returns: ["seq1", "seq3"]
```

### `to_rep_gen()`
- Returns a **generator** of representative sequences for memory-efficient processing.
- **Parameters:**
    - `with_seq` (bool): If True, yields tuples of (header, sequence). If False, yields only headers.
- **Yields:** Representative sequences or tuples one at a time.

**Example:**
```python
# Memory-efficient processing of representatives
for rep_header, rep_seq in cluster_result.to_rep_gen(with_seq=True):
    print(f"Representative: {rep_header}, Length: {len(rep_seq)}")
```

### `split_rep_as_list()`
- Splits cluster representatives into train, validation, and test sets as lists.
- **Parameters:**
    - `train` (float): Proportion for training set
    - `val` (float): Proportion for validation set  
    - `test` (float): Proportion for test set
    - `with_seq` (bool): Include sequences in output
    - `shuffle` (bool): Shuffle data before splitting
    - `seed` (int): Random seed for reproducibility
- **Returns:** Tuple of three lists (train, validation, test)
- **Note:** Proportions are automatically normalized if they don't sum to 1.0

**Example:**
```python
# Split representatives for machine learning
train_reps, val_reps, test_reps = cluster_result.split_rep_as_list(
    train=0.7, 
    val=0.15, 
    test=0.15,
    with_seq=True,
    shuffle=True,
    seed=42
)
```

### `split_rep_as_fasta()`
- Splits cluster representatives into train, validation, and test sets and saves as FASTA files.
- **Parameters:**
    - `train` (float): Proportion for training set
    - `val` (float): Proportion for validation set
    - `test` (float): Proportion for test set  
    - `shuffle` (bool): Shuffle data before splitting
    - `seed` (int): Random seed for reproducibility
- **Returns:** Tuple of file paths (train_path, val_path, test_path)
- **Note:** Files are saved with suffixes `_rep_train.fasta`, `_rep_val.fasta`, `_rep_test.fasta` and in the parent directory of the clustering outputs.

**Example:**
```python
# Create train/val/test FASTA files for ML workflows
train_file, val_file, test_file = cluster_result.split_rep_as_fasta(
    train=0.8,
    val=0.1, 
    test=0.1,
    shuffle=True,
    seed=42
)
print(f"Training set saved to: {train_file}")
print(f"Validation set saved to: {val_file}")
print(f"Test set saved to: {test_file}")
```

### `to_path()`
- Returns a list of relevant output file paths:
    - `cluster_prefix_all_seqs.fasta`: All sequences in clusters
    - `cluster_prefix_cluster.tsv`: Cluster information
    - `cluster_prefix_rep_seqs.fasta`: Representative sequences

**Note:** `to_list()`, `to_pandas()`, `to_gen()`, `to_rep_list()`, and `to_rep_gen()` rely on `cluster_prefix_all_seqs.fasta`.

## For Basic Users
Using the `pymmseqs.commands.easy_cluster` command provides an `EasyClusterParser` object.

Example:
```python
from pymmseqs.commands import easy_cluster

my_cluster = easy_cluster(
    fasta_files="query.fasta",
    cluster_prefix="output/cluster",
    tmp_dir="output/tmp",
    min_seq_id=0.3,
)

# If you have a small dataset...
cluster_list = my_cluster.to_list()
print(cluster_list[:2])

# Get just the representative sequences
representatives = my_cluster.to_rep_list(with_seq=True)
print(f"Found {len(representatives)} clusters")

# Split representatives for machine learning workflows
train_reps, val_reps, test_reps = my_cluster.split_rep_as_list(
    train=0.7, val=0.15, test=0.15, 
    with_seq=True, shuffle=True, seed=42
)
print(f"Train: {len(train_reps)}, Val: {len(val_reps)}, Test: {len(test_reps)}")

# Or save directly as FASTA files
train_file, val_file, test_file = my_cluster.split_rep_as_fasta(
    train=0.8, val=0.1, test=0.1, shuffle=True, seed=42
)
```

## For Advanced Users
Advanced users can utilize the `pymmseqs.config.EasyClusterConfig` object for additional control.

Example:
```python
from pymmseqs.config import EasyClusterConfig
from pymmseqs.parsers import EasyClusterParser

# Create the configuration
easy_cluster_config = EasyClusterConfig(
    fasta_files="data/query.fasta",
    cluster_prefix="output/cluster",
    tmp_dir="output/tmp",
    min_seq_id=0.3,
)

# Execute the configuration
easy_cluster_config.run()

# Obtain the parser object
easy_cluster_parser = EasyClusterParser(easy_cluster_config)

# Create a generator for clusters
cluster_gen = easy_cluster_parser.to_gen()

# Retrieve the representative sequence of a cluster with more than 10 members
for cluster in cluster_gen:
    if len(cluster["members"]) > 10:
        print(cluster["rep"])
        break

# Work with representative sequences efficiently
for rep_header, rep_seq in easy_cluster_parser.to_rep_gen(with_seq=True):
    if len(rep_seq) > 500:  # Process long sequences
        print(f"Long representative: {rep_header}")

# Split data for machine learning with custom configuration
train_reps, val_reps, test_reps = easy_cluster_parser.split_rep_as_list(
    train=0.6, val=0.2, test=0.2,
    with_seq=True, shuffle=True, seed=123
)

# Save splits as FASTA files for downstream analysis
train_file, val_file, test_file = easy_cluster_parser.split_rep_as_fasta(
    train=0.7, val=0.15, test=0.15, shuffle=True, seed=456
)
```

---

# [EasyTaxonomyParser](https://github.com/heispv/pymmseqs/blob/master/pymmseqs/parsers/easy_taxonomy_parser.py)
The `EasyTaxonomyParser` provides methods for accessing taxonomy classification results as Python objects, computing taxonomic composition and diversity metrics, and generating visualizations.

## Methods:

### `to_pandas(output='lca')`
- Returns a pandas DataFrame for any of the four output files.
- **Parameters:**
    - `output` (str): Which file to parse — `'lca'`, `'report'`, `'tophit_aln'`, or `'tophit_report'`.

### `to_list(output='lca')`
- Returns a list of dictionaries from the specified output file.

### `to_gen()`
- Returns a generator that yields LCA assignments one row at a time with automatic type conversion.

### `to_path()`
- Returns a dict mapping output names to file paths:
    - `'lca'`: LCA assignments TSV
    - `'report'`: Kraken-style taxonomic report
    - `'tophit_aln'`: Top-hit alignment file
    - `'tophit_report'`: Top-hit coverage report

### `to_json(path=None)`
- Exports LCA results as JSON. Returns a string if `path` is None, otherwise writes to file.

### `to_csv(path, output='lca')`
- Exports the specified output as a CSV file.

### `summary()`
- Returns a dict with overall classification statistics:
    - `total_queries`, `classified`, `unclassified`, `classified_pct`, `num_taxa`, `top_phyla`

### `report()`
- Parses the Kraken-style taxonomic report as a DataFrame with columns: `percentage`, `num_reads_clade`, `num_reads_direct`, `rank_code`, `taxid`, `name`, `depth`.

### `lca_assignments()`
- Returns per-query LCA assignments as a DataFrame with columns: `query`, `taxid`, `rank`, `taxon_name`.

### `top_hits()`
- Returns the top-hit alignment data as a DataFrame. Column format depends on the `format_output` config setting.

### `composition(rank='phylum')`
- Returns taxonomic composition at a given rank as a DataFrame with columns: `taxon`, `count`, `proportion`.
- Supports all standard ranks: `domain`, `kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species`, etc.

**Example:**
```python
comp = tax_result.composition('phylum')
# Returns:
#             taxon  count  proportion
# 0  Pseudomonadota    151    0.310062
# 1        Chordata     70    0.143737
# 2       Bacillota     58    0.119097
```

### `diversity(rank='phylum')`
- Computes alpha diversity metrics at a given rank.
- Returns a dict with: `richness`, `shannon_entropy`, `simpson_index`, `evenness`.

**Example:**
```python
div = tax_result.diversity('genus')
# Returns: {'richness': 202, 'shannon_entropy': 4.924, 'simpson_index': 0.989, 'evenness': 0.928}
```

### `rank_summary()`
- Returns a DataFrame showing the number of taxa and reads at each taxonomic rank, sorted by canonical rank order.

### `unclassified_report()`
- Returns a dict with classified vs unclassified statistics: `total_queries`, `classified`, `unclassified`, `classified_pct`, `unclassified_pct`.

### `filter_by_taxon(name, rank='phylum')`
- Filters LCA assignments to queries classified at or below the specified taxon.
- Uses the report hierarchy to identify all descendant taxids.

**Example:**
```python
chordata_queries = tax_result.filter_by_taxon('Chordata', 'phylum')
# Returns DataFrame with only queries classified within Chordata
```

### `plot_composition(rank='phylum', top_n=10)`
- Horizontal bar chart of the top taxa at a given rank.

### `plot_composition_pie(rank='phylum', top_n=8)`
- Pie chart with an "Other" slice for remaining taxa.

### `plot_classified_vs_unclassified()`
- Pie chart showing classified vs unclassified query proportions.

### `plot_rank_resolution()`
- Bar chart showing how many reads were assigned at each taxonomic rank.

### `plot_diversity_comparison(ranks=None)`
- Grouped bar chart comparing Shannon entropy and Simpson's index across multiple ranks.

**Visualization notes:** All plot methods accept an optional `ax` parameter for subplot composition, import matplotlib lazily, and return `(fig, ax)` tuples for further customization.

## For Basic Users
When using the `pymmseqs.commands.easy_taxonomy` command, you receive an `EasyTaxonomyParser` object.

Example:
```python
from pymmseqs.commands import easy_taxonomy

tax_result = easy_taxonomy(
    fasta_file="query.fasta",
    target_db="swissprotDB",
    tax_reports="output/tax_result",
)

# Quick overview
print(tax_result)
# EasyTaxonomyParser:
#   total_queries: 500
#   classified: 497
#   classified_pct: 99.4
#   num_taxa: 1356
#   top_phyla: Pseudomonadota, Chordata, Bacillota, Ascomycota, Streptophyta

# Taxonomic composition at any rank
comp = tax_result.composition('class')
print(comp.head(5))

# Diversity metrics
for rank in ['phylum', 'class', 'order', 'family', 'genus', 'species']:
    d = tax_result.diversity(rank)
    print(f"{rank}: Shannon={d['shannon_entropy']:.2f}, richness={d['richness']}")

# Filter to a lineage
mammals = tax_result.filter_by_taxon('Mammalia', 'class')
print(f"Mammalian queries: {len(mammals)}")

# Export
tax_result.to_csv('results.csv')

# Visualize
fig, ax = tax_result.plot_composition('phylum', top_n=12)
fig.savefig('phylum_composition.png', dpi=150)

fig, ax = tax_result.plot_diversity_comparison()
fig.savefig('diversity.png', dpi=150)
```

## For Advanced Users
Advanced users can utilize the `pymmseqs.config.EasyTaxonomyConfig` object for full control over all 80+ parameters.

Example:
```python
from pymmseqs.config import EasyTaxonomyConfig
from pymmseqs.parsers import EasyTaxonomyParser

config = EasyTaxonomyConfig(
    fasta_file="query.fasta",
    target_db="swissprotDB",
    tax_reports="output/tax_result",
    tmp_dir="output/tmp",
    lca_mode=4,           # top-hit mode
    s=7.5,                # high sensitivity
    e=1e-5,               # strict e-value
    tax_lineage=1,        # include full lineage in LCA output
    report_mode=0,        # Kraken-style report
)

config.run()

parser = EasyTaxonomyParser(config)

# Access different output files
lca_df = parser.to_pandas('lca')
report_df = parser.to_pandas('report')
tophit_df = parser.to_pandas('tophit_aln')

# Create a multi-panel figure
import matplotlib.pyplot as plt
fig, axes = plt.subplots(1, 3, figsize=(24, 7))
parser.plot_composition('phylum', top_n=10, ax=axes[0])
parser.plot_classified_vs_unclassified(ax=axes[1])
parser.plot_diversity_comparison(ax=axes[2])
fig.savefig('taxonomy_overview.png', dpi=150, bbox_inches='tight')
```

---

# [SearchParser](https://github.com/heispv/pymmseqs/blob/master/pymmseqs/parsers/search_parser.py)
The `SearchParser` processes the output of the MMseqs2 `search` command, which performs sequence similarity searches between a query database and a target database.

## Methods:

### `to_list()`
- Returns a list of dictionaries, each representing a row in the alignment file.
- Each dictionary contains keys corresponding to the alignment data, such as:
  - `query`
  - `target`
  - `fident`
  - `alnlen`
  - `mismatch`
  - `gapopen`
  - `qstart`
  - `qend`
  - `tstart`
  - `tend`
  - `evalue`
  - `bits`

### `to_tsv()`
- Returns nothing, just creates a new .tsv file.
- This is like running a `mmseqs convertalis` function with `format_mode=4` under the hood.

### `to_pandas()`
- Returns a pandas DataFrame containing the alignment data.
- The columns are the same as the keys in the dictionaries returned by `to_list()`.

### `to_gen()`
- Returns a generator of dictionaries, each representing a row in the alignment file.
- The keys are the same as the ones in the dictionaries returned by `to_list()`.

### `to_path()`
- Returns the path to the alignment database.

**Note:** When using `to_list()`, `to_pandas()`, or `to_gen()` methods, the `SearchParser` automatically runs the `to_tsv` command which is equal to `mmseqs convertalis` command to convert the binary alignment database to a readable TSV format if needed.

## For Basic Users
When using the `pymmseqs.commands.search` command, you receive a `SearchParser` object.

Example:
```python
from pymmseqs.commands import search
search_result = search(
query_db="output/query_db",
target_db="output/target_db",
    alignment_db="output/search_results",
    tmp_dir="output/tmp",
    min_seq_id=0.3,
    max_seqs=1000,
)

# Get all alignments with sequence identity > 50%
filtered_alignments = []
for alignment in search_result.to_gen():
    if float(alignment["fident"]) > 0.5:
        filtered_alignments.append(alignment)
print(f"Found {len(filtered_alignments)} alignments with >50% identity")
```

## For Advanced Users
Advanced users can utilize the `pymmseqs.config.SearchConfig` object for additional control.

Example:
```python
from pymmseqs.config import SearchConfig
from pymmseqs.parsers import SearchParser

# Create the configuration with advanced parameters
search_config = SearchConfig(
    query_db="output/query_db",
    target_db="output/target_db",
    alignment_db="output/search_results",
    tmp_dir="output/tmp",
    e=1e-5,
    min_seq_id=0.3,
    c=0.8,
    a=True,
)

# Turn of the logging file saving
search_config.set_config_options(has_log=False)

# Execute the configuration
search_config.run()

# Obtain the parser object
search_parser = SearchParser(search_config)

# Convert results to pandas DataFrame for analysis
results_df = search_parser.to_pandas()

# Filter and analyze results
high_quality_hits = results_df[(results_df['fident'] > 0.7) & (results_df['evalue'] < 1e-10)]
print(f"Found {len(high_quality_hits)} high-quality alignments")
```

---

# [EasySearchParser](https://github.com/heispv/pymmseqs/blob/master/pymmseqs/parsers/easy_search_parser.py)
One of the main differences between the `EasySearchParser` and the `SearchParser` is that the `EasySearchParser` also accepts the inputs as a fasta files, but in the case of the `SearchParser` you need to pass the database paths.

* When we are running the `easy_search` command it will run a `EasySearchConfig` under the hood, with the `format_mode` set to 4. So, we get a .tsv file as output with headers which can be later parsed by the `EasySearchParser`.

## Methods:

### `to_list()`
- Returns a list of dictionaries, each representing a row in the alignment file.
- Each dictionary contains the following keys (these are the defined in the `format_output` parameter of the `EasySearchConfig` object):
    - `query`
    - `target`
    - `fident`
    - `alnlen`
    - `mismatch`
    - `gapopen`
    - `qstart`
    - `qend`
    - `tstart`
    - `tend`
    - `evalue`
    - `bits`


### `to_pandas()`
- Returns a pandas DataFrame containing the alignment data.
- The columns are the same as the keys in the dictionaries returned by `to_list()`.

### `to_gen()`
- Returns a generator of dictionaries, each representing a row in the alignment file.
- The keys are the same as the ones in the dictionaries returned by `to_list()`.

### `to_path()`
- Returns a list of file paths for the output files.

## Basic Users
When using the `pymmseqs.commands.easy_search` command, you get a `EasySearchParser` object.
```python
from pymmseqs.commands import easy_search

my_search = easy_search(
    query_fasta="query.fasta",
    target_fasta_or_db="target.fasta",
    alignment_file="output/results",
    tmp_dir="output/tmp",
    min_seq_id=0.5,
    max_seqs=10000,
)

```

Since you get a parser object...

Let's say you want to get all the alignments with a sequence identity greater than 50%.
```python
filtered_alignments = []
for alignment in my_search.to_gen():
    if alignment["fident"] > 0.5:
        filtered_alignments.append(alignment)
```

Maybe your dataset is too large to fit into memory...

```python
def get_filtered_alignments(parser, threshold=0.5):
    for alignment in parser.to_gen():
        if alignment["fident"] > threshold:
            yield alignment
```

Then just use the generator to get the filtered alignments, and you don't need to worry about memory issues.

## Advanced Users
Advanced users can utilize the `pymmseqs.config.EasySearchConfig` object for additional control.

Example:
```python
from pymmseqs.config import EasySearchConfig
from pymmseqs.parsers import EasySearchParser

# Create the configuration
easy_search_config = EasySearchConfig(
    query_fasta="query.fasta",
    target_fasta_or_db="target.fasta",
    alignment_file="output/results",
    tmp_dir="output/tmp",
    realign=True,
    format_mode=4,
)

# Note that if you don't set the format_mode to 4, you will get an error.

# Execute the configuration
easy_search_config.run()

# Obtain the parser object
easy_search_parser = EasySearchParser(easy_search_config)

# Get the alignment list
alignment_list = easy_search_parser.to_list()

```
