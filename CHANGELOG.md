### 1.7.2

Fixing fix

### 1.7.1

#### Bugfix

Fixing some issues with blast report generation when using newer BLAST databases. Source of the problemis still unclear but appears to be linked to new realeses of the BLAST db and concerns only specific taxids or sequences.

### 1.7.0

#### Breaking changes

Older configuration files are not compatible anymore. Update the files by adding the following line
should yield the same results as before:

```{yaml}
blast_filter_low_complexity: True
```

#### New features

It is now possible to desactivate the defualt low-complexity filter of the BLAST search.
This can be advantageous if you expect your barcode to contain low-complexity sequences which could 
prevent getting any match at all.
This behaviour can be activated/deactivated by changing the `blast_filter_low_complexity` from `True` to `False`.

The default behaviour (`False`) uses the default 'DUST' filter of the blast tool:
`-dust 20 64 1 -soft_masking true`.

### 1.6.6

#### Fixes

- Corrects parsing of the `trim_primers_3end` parameter (#64)
- Added a primer disambiguation step that converts primers sequences in the IUPAC ambiguous nucleotide format to their explicit forms (#63)

### 1.6.5

#### Fixes

- BLAST rule now correctly uses the `threads_sample` parameter instead of `threads`. This results in better ressource management for the BLAST rule.

### 1.6.4

Benchmarking paper for 16S Metabarcoding of meat products is online:

Denay, G.; Preckel, L.; Petersen, H.; Pietsch, K.; Wöhlke, A.; Brünen-Nieweler, C. 
Benchmarking and Validation of a Bioinformatics Workflow for Meat Species Identification Using 16S rDNA Metabarcoding. 
Foods 2023, 12, 968. https://doi.org/10.3390/foods12050968 

#### New features

- It is now possible to filter specific sequences form the database using a list of accession provided in a text file using the `seq_blocklist` parameter (#60).

#### Improvements

- Taxa names are now displayed in the benchmarking report
- Fixed many typos and errors in the documentation
- Improved report aesthetics
- Improved the fetch_nt_blast.sh script to make it easier to resume interupted processes, also more robust

#### Fixes

- Fixed confusion matrix calculation when the expected taxid rank is above the target rank

### 1.6.3

#### Improvements
- Benchmarking: confusion table now reports prediction rank as well

#### Bug fix

- Fixed major problem in confusion matrix determnination for the benchmark module. Prior to this fix, False negatives were not correctly reported.
- Small fix to report aggregation rules for a rarely happening failure
- Fixed composition summary when input is empty

### 1.6.2

#### Bug fix

- Moved the benchmark module form a rule to a workflow. THis allows to ignore the benchmark arguments when running routine
  analysis, as originally intended. Benchmarking is now called with the  `-s path/to/FooDMe/woorkflow/benchmark` argument.
  This has no impact on the basic analysis (no `-s` argument) or the paramspace analysis (`-s path/to/FooDMe/woorkflow/paramspace`).

#### Documentation

- Fixed a few errors in the `config.yaml` comments
- Modified documentation of the `benchmark` module

### 1.6.1

#### Bug fix

- Fixed handling of last common node calculation in confusion matrix where last common ancesotr is the root node.
- Modified dependencies in environement definition files to solve some issues in conda solving

### 1.6.0

**This update is not backwards compatible.**
**A configuration file update is nescessary.**

#### New features

- Added a new Snakefile for parameter space exploration.
  Basically acts as a wrapper around the foodme benchmark workflow 
  for parameter grid search using snakemake's `Paramspcae` utility.

#### Non backward compatible changes

- Flattened the parameter structure in the configuration.
  This is more compatible with the `--config` CLI argument and
  was required for the implementation of the parameter space exploration workflow.
  This requires users to update their configurations.

#### Documentation

- Added documentation for the `paramspace` workflow.

### 1.5.1

#### Bug fix

- Added missing parameters to meat config file
- Fix dtype parsing in confusion matrix calculation
- Fix package version reporting
- Fix Error calculation in benchmarking module
- Fix multiple plotting in benchmarking report

### 1.5.0

**This update is not backwards compatible.**
**A configuration file update is nescessary.**

#### New features

- Benchmark module is live with possibility to compare results to an expected sample
  composition. The benchmark module will output the comparison results and several useful 
  metrics in an HTML report. It can be used directly for validation or parameter space exploration.

#### Non backward compatible changes

- Added the benchmarking module which can be called with `snakemake benchmark`
- Added required parameters in the config file
- The python wrapper is now deprecated. See the documentation on how to use configuration files.

#### Dependencies

- Added a dependency to Scikit-learn
- Updated R packages dependency in the `rmarkdown` environment

#### Documentation

- Added documentation for the benchmark module

### 1.4.9

#### Bug fix

- Added missing `pandas` dependency in `taxidTools` environment
- Moved log directive to top of python scripts to catch import errors
- Replaced all `bc` callsby `printf` statements in `vsearch.smk` (#52)
- Improved logging for OTU workflow
- Added test suite for OTU workflow

#### Documentation

- Improved Conda installation guide by quoting the Bioconda guide and adding new snakemake requirement to set sstrict channel priority

### 1.4.8

#### Bug fix

- Fixed header in `consensus-table.tsv`
- Fixed a bash synthax misuse in the calculation of VSearch statistics

#### Documentation

- Moved the documentation to the homepage at https://cvua-rrw.github.io/FooDMe/

### 1.4.7

#### Bug fix

- Fixed missing report-wise reports

### 1.4.6

#### Bug fix

- The parameter `taxid_filter` now only accepts integers, default config values have been changed (#42).

#### Improvements

- Now correctly reports composition as both percentage of total usable reads and assigned reads (#41) 
- Added a configuration file for 16S birds and mammals experiments

### 1.4.5

#### Documentation

- The usage of the python warpper is not recommended. Prefer the use a yaml configuration file. 
- Pending deprecation warning added to the python wrapper
- Expanded documentation on the use of the config file 

#### Improvements

- Improved error handling and logging for the DADA2 steps. Will now correctly output number of reads and denoisin/merging results for failing samples.

### 1.4.4

#### Improvements

- Now unpacks trimmed read files on a sample wise fashion prior to Dada2 denoising instead of unpacking all samples at once. This should reduce the memory use during the analysis.
- Preventively fixed a pandas CopyWarning (#31) and FutureWarning
- Updated dependencies to newer versions. NCBI's upcoming new identifier definitions should be supported (#33).
- Check compatibility with snakemake v7 (#34)
- Dependency taxidTools now handled through conda environment and therefore not needed in the base environment anymore (#36). 
- Reorganised logging (#27)
- Fully linted and formatted (#28)

### 1.4.3

#### Bug Fix

- Fixed a variable refernece breaking Vsearch pipeline
- Fixed time display upon pipeline completion on success or error

### 1.4.2

#### Bug Fix

- Fixed wrapper

### 1.4.1

#### Improvements:

- Moderate performance improvements due to saving taxonomy as a filtered JSON file. Expect the workflow to be about 1 min faster per sample.

#### Bug fixes:

- Fixed Github version paring for lightweight tags. 

#### Standardization

- `tests` was renamed `.tests`
- Linting and reorganize workflow to match be closer to snakemake standards
- Added JSON-Schema validation for the config and sample sheet files
- Added the possibility to export a Snakemake report containing QC summaries and results as well as the workflow runtime and DAG using the `--report` argument (snakemake CLI only)

### Version 1.4.0

#### Implementation changes:

- Migrated to TaxidTools version 2. The taxidTools package must now be installed via conda or pip before starting th epipeline (See README.md).
- Modified default parameters of the config file and python laucher with more sensible values

#### Improvements:

- Expand disambiguation info with the frequency of each species (#17)
- Add minimum consensus filter as an alternative to last common ancestor. Use it with the parameter `--min_consensus`. The value be be in the interval (0.5;1], 1 being a last common ancestro behavior and 0.51 a simple majority vote.
- Added blocklist of taxids to mask (#13). Default behaviour is to mask extinct taxids. Users can skip this steps or provide their own blocklist with the `--blocklist` parameter.

#### Bug fixes

- Cluster that do not find a matching reference in BLAST are not counted towards the compoisiton total anymore. Additionnaly the number of assigned reads is now visible in the summary report(#12)
- Fixed the calculation of the "No primer found" field under the triming statistics (#19)

### Version 1.3.5

- Upgraded Dada2 dependency to version 1.20

### Version 1.3.4

- Upgraded dependencies to last (conda) version

### Version 1.3.3

- Test now runs with just `snakemake --cores 1 --use-conda`
- Added CI in github actions
- Reworked environments definition files, environments should build correctly.

### Version 1.3.2

- Added a very basic test script. This is meant to test the installation - not provide unit testing
- Added an example of expected output
- Fixed Snakemake version
- Fixed summary report

### Version 1.3.1

- Workflow will no longer crash on blank samples
- Now reports the proportion of reads discarded during primer trimming

### Version 1.3.0

#### Incompatible changes

- Now requires user to provide a fasta file with primer sequences

#### New features

- Taxonomic reports now include a 'disambiguation' field summarizing the different blast hit for each cluster
- Primers will now be trimmed for the reads before quality trimming. It is possible to trim primers on both ends

#### Minor changes

- Performance fix for the display of large tables in the html report

### Version 1.2.2:

- Updated BLAST+ and Fastp to the latest version
- Report now includes links to BLast reports
- Blast report now includes number of mismatch, gaps and alignment length
- Added the --skip_adapter_trimming option to disable adapter trimming in fastp (only recommended for artificial dataset)

### Version 1.2.1:

- taxidTools is now a submodule
- Cloning the repository should now be done with '--recurse-submodules'
- taxidtools updated to version 2
- Adapted scripts to the new version of taxidtools
- Changed BLAST database masking to not be silent about taxids missing from the Taxdump definition files

### Version 1.2.0:

#### New features

- Added the option to filter the BLAST search by taxid

#### Fixes

- Fixed a performance issue for BLAST filtering
- Better error handling for LCA determination
- Snakemake logging has been moved to the logs folder

### Version 1.1.0:

#### New features

- Added primer trimming option (experimental)
- Added subspecies to taxonomy levels
- Added a helper script to fetch the BLAST nt database

#### Fixes

- Fixed Krona broken link in report
- Fixed BLAST filtering for floating point values of bitscores
- Fixed crash upon absence of BLAST hits
- Fixed BLAST database version reporting

### Version 1.0.0:

- initial release
