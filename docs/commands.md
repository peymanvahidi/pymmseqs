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
- `createdb_mode`: int = 0,
- `dbtype`: int = 0

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
