# Benchmark mode

FooDMe offers a benchmarking module that can be used to
automatically compare the anaylsis results to a theoritical 
sample composition and compute performance indicators.

This can be used directly for method validation or verification
or to compare the performance of the pipeline with different 
parameter sets.

## Benchmarking parameters

The benchmarking module requires three additional paramters:

* `reference`: Path to a table containing the sample composition information containing the followuing fields:
    * `sample`: Sample name
    * `taxid`: Taxonomic identifier
    * `proportion`: Expected fraction of the sample made up from this component, in the interval [0, 1]
* `threshold`: A minimal quantity for a component to be considered a 'true' result, given in percent.
* `target_rank`: A maximum taxonomic rank for a component to be considered a 'true' result

!!! warning
    
    Because of limitaitons of the NCBI taxonomiy classification, only
    Linnaean ranks are supported. This means you can choose only from 
    the following (case-sensitive) options:
    * species
    * genus
    * family
    * order
    * class
    * phylum
    * kingdom

## Usage

The basic usage only differs to the basic analysis by specifying the `benchmark` 
target directly after the `snakemake` call:

```bash
conda activate snakemake
snakemake benchmark \
  --use-conda --conda-prefix ~/conda-envs --cores 1 \
  --configfile ~/FooDMe/config/myconfig.yaml
```

This will either run the normal analysis followed by the benchmarking module,
or just the benchmarking module if the analysis already exists.

## Benchmarking results

### Yield

### Metrics

### Confusion matrix

