# Parameter space exploration

Tutorial for Snakemake Paramspace coming soon

## Parameter grid


## Configuration

| Parameter                 | Expected values           | Description |
| ---                       | ---                       | --- |
| `workdir`                 | Path                      | Path to the output directory, will be created if <br>it doesn´t exist |
| `foodme_config` | Path | Path to the foodme configuration file. |
| `paramspace` | Path | Path to the parameter grid |
| `force_rerun` | True/False | Whether to force a recalculation of the foodme results <br> for all parameter combinations. <br> Equivalent to the `--forceall` directive.|

!!! warning
    
    In the parameter space mode, the foodme configuration file only supports
    **absolute** file paths. Relative paths will result in errors, you've been warned.

## Running a parameter space analyis

!!! warning
    
    The command `--forceall` will not force a re-run of foodme.
    To trigger a rerun you will have to set the `force_rerun` parameter
    to `True`in the configuration or manually delete the `foodme_runs` folder.

## Results


