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

Whereby the number of dedicated cores can be adapted to your specific system.

## Run-specific parameters

As creating or modifying the configuration file to change the sample sheet and the output 
directory for each run can be a bit cumbersome, it is possible to dynamically modify 
parameters at execution time. 

```bash
WORKDIR=$(date +%F)_foodme
SAMPLES=~/raw_data/samples.tsv

snakemake --use-conda --conda-prefix ~/conda-envs --cores 1 \
  --configfile ~/FooDMe/config/myconfig.yaml \
  --config workdir=$WORKDIR samples=$SAMPLES
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
- `--notemp` will disable temporary file removal. THis will considerably increase the size of the output be can be usefull to troubleshoot an analysis or puzzling results.
- `--report` will output run statistics.
