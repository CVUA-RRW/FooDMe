# Parameter space exploration

Parameter space exploration allows to easily try different sets of parameters,
measure their impact on analysis performance and guide method optimization for new
matrices, targets, or databases.

Provided a grid of parameter combination, foodme will be run several time with a 
different paramter set. Performance metrics will be extracted for each run using using
a expected sample composition and the metrics for each run will be aggregated.

Importantly you will need to provide a Foodme configuration file (see [configuration](../userguide/configuration)) 
containing default values and a parameter grid containing a set of parameters to update on each row.

## Parameter grid

The paraneter grid should be organized as a tab-delimited text file where
the first row contains the name of the parameters to vary. Each row then contains 
a set of parameter that will be used to update the default configuration.

For example with the file below, three independent foodme runs will be triggered with ASV clustering,
dereplication, and 97% identity clustering. The other parameters will be taken for the default configuration file.

| cluster_method | cluster_identity |
| --- | --- |
| asv | 1 |
| otu | 1 |
| otu | 0.97 |

## Configuration

The parameter space exploration mode only requires a few argument, that you can either pass to snakemake
as a YAML configuration file (see the template in `config/` or directly through the command line with 
the `--config` argument.

| Parameter                 | Expected values           | Description |
| ---                       | ---                       | --- |
| `workdir`                 | Path                      | Path to the output directory, will be created if <br>it doesn´t exist |
| `foodme_config` | Path | Path to the foodme configuration file. |
| `paramspace` | Path | Path to the parameter grid |
| `force_rerun` | True/False | Whether to force a recalculation of the foodme results <br> for all parameter combinations. <br> Equivalent to the `--forceall` directive.|

!!! note
    
    As parameter space exploration only makes sense if you can compare the results 
    to known sample compositions, it is required to provide values for the `benchmark_*`
    paramters. These values can be part of the parameter space exploration too!


!!! warning
    
    In the parameter space mode, the foodme configuration file only supports
    **absolute** file paths. Relative paths will result in errors, you've been warned.

## Running a parameter space analyis

The parameter space analysis is organized in a separate workflow from foodme, it is
therefore nescessary to point snakemake to the correct workflow definition:

```
conda activate snakemake
snakemake -s ~/FooDMe/workflow/paramspace \
  --use-conda --conda-prefix ~/conda-envs --cores 1 \
  --configfile ~/FooDMe/config/myconfig.yaml
```

!!! warning
    
    The command `--forceall` will not force a re-run of foodme.
    To trigger a rerun you will have to set the `force_rerun` parameter
    to `True`in the configuration or manually delete the `foodme_runs` folder.

## Results

The parameter space exploration will produce an aggregate of foodme `benchmark` 
results: yields, metrics, PR-curves, and confusion matrices, all with indication
of the specific parameter set that was used for each run.

Additonally the ressource usage (CPU, memory, I/O) will be measured for each foodme run
and saved as a table. See also [snakemake's documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html?benchmark-rules#benchmark-rules).

The output is organized as follows:

```
wordir/
 |- aggregated/             \\ Perfomance indicators
 |- benchmark/              \\ Ressource usage benchmarking
 |- foodme_runs/            \\ Individual foodme runs results 
 |   |- param~value/... 
 |- logs/                   \\ Log files
```



