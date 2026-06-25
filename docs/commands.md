# pymmseqs.commands
This module provides a Python wrapper around the MMseqs2 command-line tools with only the most commonly used parameters.


## Important Things to Know
- The paths are relative to the parent directory of the python script it is called from.
- You have the ability to adjust commonly used parameters for each command, while others default to their corresponding MMseqs2 command settings. For more extensive customization, refer to the [pymmseqs.config](./config.md) module to access and modify the full range of parameters available.
- After executing a command, you will see the following output in your terminal, which indicates that the command has been executed successfully.
```
-------------------- Running a mmseqs2 command --------------------
✓ Detailed execution log has been saved
✓ <Command> completed successfully
  Results saved to: <path_to_Command_results>
```
- Detailed MMseqs2 execution logs are stored in the `logs` directory located in the parent directory of your output file path.

---

## [createdb](https://github.com/heispv/pymmseqs/blob/master/pymmseqs/commands/createdb.py)
Wrapper around the `mmseqs createdb` command.

```python
from pymmseqs.commands import createdb

my_db = createdb(
  fasta_file="data/example.fasta",
  sequence_db="output/example_db"
)
```

Optional parameters:
- `shuffle`: bool = True,
- `compressed`: bool = False,
- `createdb_mode`: int = 0 (0: copy data, 1: soft-link data, 2: GPU compatible db),
- `dbtype`: int = 0,
- `gpu`: int = 0 (0: off, 1: build a GPU-compatible database using CUDA if possible),
- `threads`: Union[str, int] = 'all'

Output of `createdb` is an `CreateDBParser` object.

Methods:
- `to_path()`: Get the path prefix of the database.

[See an example](../examples/commands_module/createdb_ex.py)

---

## [createindex](https://github.com/heispv/pymmseqs/blob/master/pymmseqs/commands/createindex.py)
Wrapper around the `mmseqs createindex` command.

```python
from pymmseqs.commands import createindex

index_result = createindex(
  sequence_db="output/example_db"
)
```

Optional parameters:
- `tmp_dir`: Path = None,
- `s`: float = 7.5,
- `k`: int = 0,
- `v`: int = 3,
- `threads`: int = 14,
- `compressed`: bool = False,
- `create_lookup`: int = 0,
- `search_type`: int = 0,
- `headers_split_mode`: int = 0,
- `max_seqs`: int = 300,
- `max_seq_len`: int = 65535

- Note: if `tmp_dir` is None, tmp folder would be created in the parent dir of `sequence_db`
Output of `createindex` is a `CreateIndexParser` object.

Methods:
- `to_path()`: Get the path prefix of the indexed database.


---

## [easy_cluster](https://github.com/heispv/pymmseqs/blob/master/pymmseqs/commands/easy_cluster.py)
Wrapper around the mmseqs easy-cluster command.

```python
from pymmseqs.commands import easy_cluster

cluster_result = easy_cluster(
  fasta_files="data/example.fasta",
  cluster_prefix="output/example_clusters"
)
```

Optional parameters:
- `tmp_dir`: Path = None,
- `min_seq_id`: float = 0.0,
- `s`: float = 4.0,
- `c`: float = 0.8,
- `cov_mode`: int = 0,
- `e`: float = 0.001,
- `cluster_mode`: int = 0,
- `gpu`: int = 0 (0: off, 1: use GPU (CUDA) if possible),

- Note: if `tmp_dir` is None, tmp folder would be created in the parent dir of `cluster_prefix`
Output of `easy_cluster` is an `EasyClusterParser` object.

Methods:
- `to_path()`: Get the paths to the cluster output files.
- `to_list()`: Returns a list of clusters with their representative and member sequences.
- `to_pandas()`: Converts cluster data into a pandas DataFrame for easy manipulation.
- `to_gen()`: Generator that yields clusters one at a time for efficient processing of large files.
- `to_rep_list()`: Returns a list of representative sequences (with or without sequences).
- `to_rep_gen()`: Generator that yields representative sequences one at a time.
- `split_rep_as_list()`: Splits cluster representatives into train/validation/test sets as lists.
- `split_rep_as_fasta()`: Splits cluster representatives into train/validation/test sets and saves as FASTA files.

[See an example](../examples/commands_module/easy_cluster_ex.py)

---

## [easy_linclust](https://github.com/heispv/pymmseqs/blob/master/pymmseqs/commands/easy_linclust.py)
Wrapper around the mmseqs easy-linclust command. A faster, less sensitive linear-time alternative to `easy_cluster`, suited to very large datasets.

```python
from pymmseqs.commands import easy_linclust

cluster_result = easy_linclust(
  fasta_files="data/example.fasta",
  cluster_prefix="output/example_clusters"
)
```

Optional parameters:
- `tmp_dir`: Path = None,
- `min_seq_id`: float = 0.0,
- `c`: float = 0.8,
- `cov_mode`: int = 0,
- `e`: float = 0.001,
- `cluster_mode`: int = 0,
- `kmer_per_seq`: int = 21,
- `kmer_per_seq_scale`: str = "aa:0.000,nucl:0.200",
- `gpu`: int = 0 (0: off, 1: use GPU (CUDA) if possible),

- Note: if `tmp_dir` is None, tmp folder would be created in the parent dir of `cluster_prefix`
Output of `easy_linclust` is an `EasyLinClustParser` object.

Methods:
- `to_path()`: Get the paths to the cluster output files.
- `to_list()`: Returns a list of clusters with their representative and member sequences.
- `to_pandas()`: Converts cluster data into a pandas DataFrame for easy manipulation.
- `to_gen()`: Generator that yields clusters one at a time for efficient processing of large files.
- `to_rep_list()`: Returns a list of representative sequences (with or without sequences).
- `to_rep_gen()`: Generator that yields representative sequences one at a time.
- `split_rep_as_list()`: Splits cluster representatives into train/validation/test sets as lists.
- `split_rep_as_fasta()`: Splits cluster representatives into train/validation/test sets and saves as FASTA files.

---

## [search](https://github.com/heispv/pymmseqs/blob/master/pymmseqs/commands/search.py)
Wrapper around the mmseqs search command.

```python
from pymmseqs.commands import search

search_result = search(
  query_db="query_db",
  target_db="target_db",
  alignment_db="output/search_db"
)
```

Optional parameters:
- `tmp_dir`: Path = None,
- `s`: float = 5.7,
- `e`: float = 0.001,
- `min_seq_id`: float = 0.0,
- `c`: float = 0.0,
- `cov_mode`: int = 0,
- `a`: bool = False,
- `max_seqs`: int = 300,
- `threads`: int = 14,
- `compressed`: bool = False,

- Note: if `tmp_dir` is None, tmp folder would be created in the parent dir of `alignment_db`
Output of `search` is an `SearchParser` object.

Methods:
- `to_list()`: Returns a list of search results with their query, target, and alignment information.
- `to_pandas()`: Converts search results into a pandas DataFrame for easy manipulation.
- `to_gen()`: Generator that yields search results one at a time for efficient processing of large files.
- `to_path()`: Get the paths to the search output files.

* Note: When using the `to_list()`, `to_pandas()`, or `to_gen()` methods, internally the `SearchParser` object will run a `mmseqs convertalis` command to turn the results into a .tsv file, and then performs the necessary parsing of that file. Check out the [pymmseqs.parsers.SearchParser](./parsers.md#SearchParser) class for more details.

[See an example](../examples/commands_module/search_ex.py)

---

## [easy_search](https://github.com/heispv/pymmseqs/blob/master/pymmseqs/commands/easy_search.py)
Wrapper around the mmseqs easy-search command.

```python
from pymmseqs.commands import easy_search

search_result = easy_search(
  query_fasta="query.fasta",
  target_fasta_or_db="target.fasta",
  alignment_file="output/search_results.m8"
)
```

Optional parameters:
- `tmp_dir`: Path = None,
- `s`: float = 5.7,
- `e`: float = 0.001,
- `min_seq_id`: float = 0.0,
- `c`: float = 0.0,
- `max_seqs`: int = 300,
- `translate`: bool = False,
- `translation_table`: int = 1,
- `translation_mode`: int = 0,
- `search_type`: int = 0,

- Note: if `tmp_dir` is None, tmp folder would be created in the parent dir of `alignment_file`
Output of `easy_search` is an `EasySearchParser` object.

Methods:
- `to_list()`: Returns a list of search results with their query, target, and alignment information.
- `to_pandas()`: Converts search results into a pandas DataFrame for easy manipulation.
- `to_gen()`: Generator that yields search results one at a time for efficient processing of large files.
- `to_path()`: Get the paths to the search output files.

[See an example](../examples/commands_module/easy_search_ex.py)

---

## [easy_linsearch](https://github.com/heispv/pymmseqs/blob/master/pymmseqs/commands/easy_linsearch.py)
Fast, less-sensitive linear-time homology search (the fast sibling of `easy_search`).

```python
from pymmseqs.commands import easy_linsearch

search_result = easy_linsearch(
  query_fasta="query.fasta",
  target_fasta_or_db="target.fasta",
  alignment_file="output/search_results.m8"
)
```

Optional parameters:
- `tmp_dir`: Path = None,
- `e`: float = 0.001,
- `min_seq_id`: float = 0.0,
- `c`: float = 0.0,
- `translate`: bool = False,
- `translation_table`: int = 1,
- `search_type`: int = 0,
- `format_output`: str = "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"

- Note: if `tmp_dir` is None, tmp folder would be created in the parent dir of `alignment_file`
Output of `easy_linsearch` is an `EasySearchParser` object.

Methods:
- `to_list()`: Returns a list of search results with their query, target, and alignment information.
- `to_pandas()`: Converts search results into a pandas DataFrame for easy manipulation.
- `to_gen()`: Generator that yields search results one at a time for efficient processing of large files.
- `to_path()`: Get the paths to the search output files.

---

## [easy_taxonomy](https://github.com/heispv/pymmseqs/blob/master/pymmseqs/commands/easy_taxonomy.py)
Wrapper around the `mmseqs easy-taxonomy` command. Assigns taxonomic labels to query sequences using LCA (Lowest Common Ancestor) against a taxonomy-enabled database.

```python
from pymmseqs.commands import easy_taxonomy

tax_result = easy_taxonomy(
    fasta_file="data/query.fasta",
    target_db="swissprotDB",
    tax_reports="output/tax_result",
)
```

Optional parameters:
- `tmp_dir`: Path = None,
- `s`: float = 4.0,
- `e`: float = 0.001,
- `lca_mode`: int = 3 (1: single search LCA, 2/3: approximate 2bLCA, 4: top hit),
- `orf_filter`: int = 0 (set 0 for short reads),
- `report_mode`: int = 0 (0: Kraken, 1: Krona),
- `search_type`: int = 0,
- `min_seq_id`: float = 0.0,
- `max_seqs`: int = 300

- Note: if `tmp_dir` is None, tmp folder would be created in the parent dir of `tax_reports`

Output of `easy_taxonomy` is an `EasyTaxonomyParser` object.

Methods:

**Data Access:**
- `to_pandas(output='lca')`: Parse any output file as DataFrame ('lca', 'report', 'tophit_aln', 'tophit_report').
- `to_list()`: Returns LCA assignments as a list of dictionaries.
- `to_gen()`: Generator that yields LCA assignments one at a time.
- `to_path()`: Returns a dict of all output file paths.
- `to_json(path=None)`: Export LCA results as JSON string or to file.
- `to_csv(path)`: Export results as CSV file.
- `len(tax_result)`: Number of query sequences.
- `print(tax_result)`: Pretty summary with classification stats.

**Analysis:**
- `summary()`: Overall stats — total queries, classified count/%, top phyla.
- `report()`: Parse the Kraken-style taxonomic report as DataFrame.
- `lca_assignments()`: Per-query LCA assignments as DataFrame.
- `top_hits()`: Top-hit alignment data as DataFrame.
- `composition(rank='phylum')`: Taxonomic composition at any rank (taxon, count, proportion).
- `diversity(rank='phylum')`: Alpha diversity metrics — Shannon entropy, Simpson index, richness, evenness.
- `rank_summary()`: Classification counts at each taxonomic level.
- `unclassified_report()`: Classified vs unclassified statistics.
- `filter_by_taxon(name, rank='phylum')`: Filter LCA assignments to a specific lineage.

**Visualization (matplotlib):**
- `plot_composition(rank, top_n)`: Horizontal bar chart of top taxa.
- `plot_composition_pie(rank, top_n)`: Pie chart with "Other" slice.
- `plot_classified_vs_unclassified()`: Pie chart of classification rate.
- `plot_rank_resolution()`: Bar chart showing reads assigned per rank.
- `plot_diversity_comparison(ranks)`: Grouped bar chart of Shannon/Simpson across ranks.

All plot methods accept optional `ax` parameter for subplot composition and return `(fig, ax)` tuples.

```python
# Quick overview
print(tax_result.summary())

# Taxonomic breakdown
comp = tax_result.composition('phylum')
print(comp)

# Diversity analysis
div = tax_result.diversity('genus')
print(f"Shannon: {div['shannon_entropy']}, Richness: {div['richness']}")

# Filter to a specific lineage
chordata = tax_result.filter_by_taxon('Chordata', 'phylum')

# Generate plots
fig, ax = tax_result.plot_composition('phylum', top_n=10)
fig.savefig('phylum_composition.png')

fig, ax = tax_result.plot_diversity_comparison()
fig.savefig('diversity.png')
```

---

## [fast_easy_search](https://github.com/heispv/pymmseqs/blob/master/pymmseqs/commands/fast_easy_search.py)
Run an easy-style search optimized for single or few queries by preloading the target index into memory.

```python
from pymmseqs.commands import fast_easy_search

search_result = fast_easy_search(
  query_fasta="query.fasta",
  target_fasta="target.fasta",
  alignment_file="output/search_results.m8"
)
```

Optional parameters:
- `tmp_dir`: Path = None,
- `s`: float = 5.7,
- `e`: float = 0.001,
- `min_seq_id`: float = 0.0,
- `c`: float = 0.0,
- `max_seqs`: int = 300,
- `translate`: bool = False,
- `translation_table`: int = 1,
- `translation_mode`: int = 0,
- `search_type`: int = 0,
- `format_output`: str = "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"

- Note: if `tmp_dir` is None, a `tmp` folder is created next to `alignment_file`.

What it does under the hood (workflow):
- Create query DB using `mmseqs createdb`.
- Create target DB using `mmseqs createdb`.
- Index the target DB using `mmseqs createindex`.
- Preload the target index into memory using `mmseqs touchdb`.
- Search using `mmseqs search` with `db_load_mode=2` (mmap).
- Convert results to a readable TSV with headers using `mmseqs convertalis` (format_mode=4).

How this differs from `easy_search`:
- `easy_search` is a single MMseqs2 command that builds what it needs and runs the search, usually loading the index from disk for each run.
- `fast_easy_search` explicitly adds a `touchdb` step to keep the target index in memory and forces mmap loading (`db_load_mode=2`) during `search`.
- This typically speeds up single or few-query workflows when the index fits into RAM; for very large query sets, copying into RAM might be less beneficial.

Output of `fast_easy_search` is an `EasySearchParser` object.

Methods:
- `to_list()`: Returns a list of search results with their query, target, and alignment information.
- `to_pandas()`: Converts search results into a pandas DataFrame for easy manipulation.
- `to_gen()`: Generator that yields search results one at a time for efficient processing of large files.
- `to_path()`: Get the path to the alignment file.

---

## [convertalis](https://github.com/heispv/pymmseqs/blob/master/pymmseqs/commands/convertalis.py)
Convert an alignment database (e.g. from `search`) into a parsable BLAST-tab table.

```python
from pymmseqs.commands import convertalis

alis = convertalis(
  query_db="query_db",
  target_db="target_db",
  alignment_db="output/search_db",
  alignment_file="output/search_results.m8"
)
```

Optional parameters:
- `format_output`: str = "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits",
- `search_type`: int = 0,
- `translation_table`: int = 1,
- `threads`: Union[str, int] = 'all'

Output of `convertalis` is an `EasySearchParser` object.

Methods:
- `to_list()`: Returns a list of search results with their query, target, and alignment information.
- `to_pandas()`: Converts search results into a pandas DataFrame for easy manipulation.
- `to_gen()`: Generator that yields search results one at a time for efficient processing of large files.
- `to_path()`: Get the path to the alignment file.

---

## [convert2fasta](https://github.com/heispv/pymmseqs/blob/master/pymmseqs/commands/convert2fasta.py)
Convert a sequence database back to a FASTA file (reverse of `createdb`).

```python
from pymmseqs.commands import convert2fasta

fasta = convert2fasta(
  sequence_db="output/example_db",
  fasta_file="output/example.fasta"
)
```

Optional parameters:
- `use_header_file`: bool = False,
- `v`: int = 3

Output of `convert2fasta` is a `Convert2FastaParser` object.

Methods:
- `to_path()`: Get the path to the FASTA file.
- `to_gen()`: Generator that yields sequences one at a time for efficient processing of large files.
- `to_list()`: Returns a list of sequences with their headers.
- `to_pandas()`: Converts the FASTA records into a pandas DataFrame for easy manipulation.

---

## [extractorfs](https://github.com/heispv/pymmseqs/blob/master/pymmseqs/commands/extractorfs.py)
Six-frame extraction of open reading frames from a nucleotide database.

```python
from pymmseqs.commands import extractorfs

orfs = extractorfs(
  sequence_db="output/contigs_db",
  orf_db="output/orfs_db"
)
```

Optional parameters:
- `min_length`: int = 30,
- `max_length`: int = 32734,
- `orf_start_mode`: int = 1,
- `contig_start_mode`: int = 2,
- `contig_end_mode`: int = 2,
- `forward_frames`: str = "1,2,3",
- `reverse_frames`: str = "1,2,3",
- `translation_table`: int = 1,
- `translate`: bool = False,
- `threads`: Union[str, int] = 'all'

Output of `extractorfs` is an `ExtractOrfsParser` object.

Methods:
- `to_pandas()`: Converts the extracted ORFs into a pandas DataFrame, yielding per-ORF rows with absolute start/end coordinates on the source contig.
- `to_list()`: Returns the extracted ORFs as a list.
- `to_gen()`: Generator that yields ORFs one at a time for efficient processing of large files.
- `to_path()`: Get the path prefix of the ORF database.
- `summary()`: Summary statistics of the extracted ORFs.
