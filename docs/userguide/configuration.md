# Configuration

Parameters can be freely set to fit specific analytical needs.
Different sets of parameters should be saved to configuration files in YAML format.

A template of such configuration is stored in the repository under `FooDME/config/config.yaml`.

The `config` folder will be populated with configuration sets for specific experiments.
Feel free to submit yours!

!!! warning 

    Path to reference files are system dependent and will still need to be changed
    even for preset configurations.

## How to use to configuration file

Modify the values of each parameters as you need (see table below).
Then save your own configuration locally (for example unde `FooDme/config/`).
This configuration can be reused for successive analysis.

## Sample sheet

The input files must be linked using a sample sheet. A template for such a file 
is available under `FooDME/config/samples.tsv`. and takes the following form:

| sample | fq1 | fq2 |
| --- | --- | --- |
| nameA | A_R1.fastq.gz | A_R2.fastq.gz |
| nameB | B_R1.fastq.gz | B_R2.fastq.gz |

Simply modify the template with your own files 
or use the provided script to automatically generate a sample sheet form FASTQ files in a folder:

```bash
bash ~/FooDMe/ressources/create_sampleSheet.sh --mode illumina --fastxDir ~/raw_data
```

This will create a file called `samples.tsv` in the `raw_data` folder.

!!! info 

    This script was developed by the Federal institute of risk assessment (BfR) in Berlin.
    More information is available in their [repository](https://gitlab.com/bfr_bioinformatics/AQUAMIS).

## List of parameters


| Category | Parameter                 | Expected values           | Description |
| --- | ---                       | ---                       | --- |
| | `workdir`                 | Path                      | Path to the output directory, will be created if <br>it doesn´t exist |
| | `samples`                 | Path                      | Path to the sample sheet                           |
| | `threads_sample`          | Number                    | Number of threads assigned to each job             |
| | `threads`                 | Number                    | Number of threads assigned to the workflow         |
| `trimming` | `length_required`         | Number                    | Minimal length of the reads after primer trimming  |
| `trimming` | `qualifier-quality-phred` | Number                    | Minimal quality value per nucleotide               |
| `trimming` | `window-size`             | Number                    | Size of the sliding window for 3´ quality trimming |
| `trimming` | `mean_quality`            | Number                    | Minimal quality thresold for sliding average       |
| `trimming` | `primers_fasta`           | Path                      | Path to the fasta file containing primer sequences |
| `trimming` | `primers_end`             | True/False                | Should primers be trimmed on the 3´ end of <br>the reads? Only relevant if the sequencing length is <br>larger than the amplicon length |
| `trimming` | `primer_error_rate`       | Number [0, 1]             | Maximum error-rate allowed for primer matching     |
| `read_filter` | `min_length`              | Number                    | Minimal length of the merge reads or <br>ASV sequences |
| `read_filter` | `max_length`              | Number                    | Maximal lenght of the merge reads or <br>ASV sequences |
| `read_filter` | `max_expected_errors`     | Number                    | Maximum number of expected errors allowed in <br>merged reads (OTUs) or trimmed reads (ASVs) |
| `read_filter` | `max_ns`                  | Number                    | Maximal number of undetermined `N` nucleotide <br>per sequence. This will automatically be set <br>to 0 for ASVs |
| `cluster` | `method`                  | `otu` or `asv`            | Clustering method |
| `cluster` | `cluster_identity`        | Number [0, 1]             | OTU identity threshold. Only for OTU |
| `cluster` | `cluster_minsize`         | Number                    | Minimal size of clusters to keep |
| `cluster` | `max_mismatch`            | Number                    | Maximum number of mismatch allowed in the<br> overlap between reads. Only for ASV |
| | `chimera`                 | True/False                | Should predicted chimeric sequences be removed? |
| `taxonomy` | `rankedlineage_dmp`       | Path                      | Path to the `rankedlineage.dmp` file from the <br>`taxdump` archive |
| `taxonomy` | `nodes_dmp`               | Path                      | Path to the `nodes.dmp` file from <br>the `taxdump` archive |
| `taxonomy` | `min_consensus`           | Number [0.51, 1]          | Minimal agreement for taxonomic consensus <br>determination. 0.51 is a majority vote and 1.0 is<br> a last common ancestor determination |
| `blast` | `blast_DB`                | Path                      | Path to the BLAST database in the form <br>`path/to/folder/db-name` |
| `blast` | `taxdb`                   | Path                      | Path to the folder containing the `taxdb`files |
| `blast` | `taxid_filter`            | Taxonomic identifier      | Node under which to perform the BLAST search. <br>Equivalent to pruning the taxonomy above <br>this node |
| `blast` | `blocklist`               | `extinct` or custom path  | Path to a list of taxonomic identifier to exclude <br>from the BLAST search |
| `blast` | `e_value`                 | Number (scientific)       | Minimal E-value threshold for the BLAST search |
| `blast` | `perc_identity`           | Number [0, 100]           | Minimal identity (in percent) between the query and <br>hit sequence for the BLAST search |
| `blast` | `qcov`                    | Number [0, 100]           | Percent of the query to be covered by the hit <br>sequence for the BLAST search |
| `blast` | `bit_score_diff`          | Number                    | Maximum bit-score difference between the best<br> and last hit to keep for each query after the <br>BLAST search |
