# pymmseqs.config
This module provides configuration classes for MMseqs2 commands, allowing advanced customization of parameters and execution options.

---

# Important Things to Know
By using the `Config` class, you can use all the parameters of the MMseqs2 commands.

## Parameter Naming Convention
When using the configuration classes, MMseqs2 command parameters are converted from hyphenated format to underscore format:
- MMseqs2 parameter: `--min-seq-id`
- Python parameter: `min_seq_id`

## Path Resolution
All file paths are automatically resolved relative to the directory of the calling script, unless absolute paths are provided.

If you want to run the following command:
```
mmseqs easy-cluster input.fasta cluster_res tmp --min-seq-id 0.9 -c 0.8
```
You can do it by using the following code:
```python
from pymmseqs.config import EasyClusterConfig

config = EasyClusterConfig(
    fasta_files="input.fasta",
    cluster_prefix="cluster_res",
    tmp_dir="tmp",
    min_seq_id=0.9,
    c=0.8,
)

config.run()
```
* Note that the required parameters can be passed as positional arguments or as keyword arguments.

## Execution Process

When you call `config.run()`, the following happens:

1. **Path Resolution**: All paths are resolved relative to the directory of the calling script
2. **Parameter Validation**: The configuration validates all parameters to ensure they meet MMseqs2 requirements
3. **Command Execution**: The MMseqs2 command is executed with your specified parameters
4. **Logging**: Detailed execution logs are saved to the `logs` directory (located in the parent directory of your output path)
5. **Terminal Output**: You'll see a concise execution summary in your terminal:
   ```
   -------------------- Running a mmseqs2 command --------------------
   ✓ Detailed execution log has been saved
   ✓ <Command> completed successfully
     Results saved to: <path_to_results>
   ```
6. **Result Storage**: All output files are saved to your specified output paths

You can pass the config object to their parser class to parse the output to python objects.

```python
from pymmseqs.parser import EasyClusterParser

parser = EasyClusterParser(config)

parser.to_gen()
```

## Available Configurations

### CreateDBConfig

### CreateTaxDBConfig

### EasyClusterConfig

### EasyLinClustConfig

### SearchConfig

### EasySearchConfig

### EasyLinSearchConfig

### EasyTaxonomyConfig

### ConvertalisConfig

### AlignConfig