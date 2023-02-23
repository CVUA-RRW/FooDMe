# Starting a run 

Now that the configuration and sample sheets are ready, we can start an 
analysis!

## Basic usage

Starting and controlling the workflow execution is done through the `snakemake` command.

```bash
conda activate snakemake
snakemake --use-conda --conda-prefix ~/conda-envs --cores 1 \
  --configfile ~/FooDMe/config/myconfig.yaml
```

!!! info

    Depending on your system you may want to assign more cores to the analysis. 
    Each additional core will allow running another job in parallel, considerably speeding up execution time.
    You can adjust the number of assiged cores with `--cores N` whereby N is the number of cores.

!!! info

    Specifying a prefix for the conda environment is not nescessary but will allow you to reuse created environments 
    between runs. Doing so will save up, time, ressources, and memory on each execution.

## Run-specific parameters

As creating or modifying the configuration file to change the sample sheet and the output 
directory for each run can be a bit cumbersome, it is possible to dynamically modify 
parameters at execution time. 

```bash
WORKDIR=$(date +%F)_foodme
SAMPLES=~/raw_data/samples.tsv

snakemake --use-conda --conda-prefix ~/conda-envs --cores 1 \
  --configfile ~/FooDMe/config/myconfig.yaml \
  --config workdir=${WORKDIR} samples=${SAMPLES}
```

The above example will use the sample sheet located under `~/raw_data/samples.tsv`
and outputs the analysis results in the `[YYYY-MM-DD]_foodme` folder, regardless 
of the value for the `workdir` and `samples` arguments in the config file. The rest of
the parameters will be taken form the configuration as usual.

## Advanced Snakemake options

A large number of options are available when running snakemake that be can 
seen using the `--help` argument or in the [online documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html).

A few notables commands are:

- `--forceall` will force the recalculation of all the workflow steps. Usefull to repeat an analysis when something went wrong.
- `--notemp` will disable temporary file removal. This will considerably increase the size of the output but can be useful to troubleshoot an analysis or puzzling results.
- `--report` will output run statistics.
