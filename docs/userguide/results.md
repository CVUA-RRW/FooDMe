# Viewing the results

The analysis produces a number of files containing quality reports, tables and graphs, 
as well as some of the processed data at different steps of the analysis.

> **_NOTE:_** By default, large files are removed during the analysis
> use snakemake´s `--notemp` to keep all files.

## Folder structure

The output folder follows the following structure:

```
workdir/
|- common/             \\ Taxonomy processing files
|- logs/
|   |- all/            \\ Logs for steps processing aggregated results
|   |- common/         \\ Logs for taxonomy produces
|   |- sample_name/    \\ Sample-wise log files
|- reports/            \\ All analysis reports
|- sample_name/        \\ Sample-wise processed data are located
    |- denoising/      \\    here, splitted in the relevant folders
    |- krona/
    |- taxonomy
    |- trimmmed
    |- reports         \\ Analysis reports for individual samples
```

## Report

All the reports are stored as tab-delimited text files in the `reports` folder.
Additionaly an HTML report is produced that contains the information of the individual files.
Tables in this report are interactive and can be filtered, columns moved or removed, and 
the tables can be exported in excel, csv, and pdf formats.

### Panel descriptions

#### Quality summary

#### Trimming statistics

#### Read filtering statistics (OTU only)

#### Clustering statistics (OTU only)

#### Denoising statistics (ASV only)

#### Taxonomic assignment

#### Taxonomic assignment statistics

#### Metabarcoding results

#### Graphical results

#### Versions
