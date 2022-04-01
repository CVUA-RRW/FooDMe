### Dev

#### Improvements
* Now unpacks trimmed read files on a sample wise fashion prior to Dada2 denoising instead of unpacking all samples at once. This should reduce the memory use during the analysis.
* Preventively fixed a pandas CopyWarning (#31) and FutureWarning

### 1.4.3

#### Bug Fix
* Fixed a variable refernece breaking Vsearch pipeline
* Fixed time display upon pipeline completion on success or error

### 1.4.2

#### Bug Fix
* Fixed wrapper

### 1.4.1

#### Improvements:
* Moderate performance improvements due to saving taxonomy as a filtered JSON file. Expect the workflow to be about 1 min faster per sample.

#### Bug fixes:
* Fixed Github version paring for lightweight tags. 

#### Standardization
* `tests` was renamed `.tests`
* Linting and reorganize workflow to match be closer to snakemake standards
* Added JSON-Schema validation for the config and sample sheet files
* Added the possibility to export a Snakemake report containing QC summaries and results as well as the workflow runtime and DAG using the `--report` argument (snakemake CLI only)

### Version 1.4.0

#### Implementation changes:

* Migrated to TaxidTools version 2. The taxidTools package must now be installed via conda or pip before starting th epipeline (See README.md).
* Modified default parameters of the config file and python laucher with more sensible values

#### Improvements:
* Expand disambiguation info with the frequency of each species (#17)
* Add minimum consensus filter as an alternative to last common ancestor. Use it with the parameter `--min_consensus`. The value be be in the interval (0.5;1], 1 being a last common ancestro behavior and 0.51 a simple majority vote.
* Added blocklist of taxids to mask (#13). Default behaviour is to mask extinct taxids. Users can skip this steps or provide their own blocklist with the `--blocklist` parameter.

#### Bug fixes
* Cluster that do not find a matching reference in BLAST are not counted towards the compoisiton total anymore. Additionnaly the number of assigned reads is now visible in the summary report(#12)
* Fixed the calculation of the "No primer found" field under the triming statistics (#19)

### Version 1.3.5

* Upgraded Dada2 dependency to version 1.20

### Version 1.3.4

* Upgraded dependencies to last (conda) version

### Version 1.3.3

* Test now runs with just `snakemake --cores 1 --use-conda`
* Added CI in github actions
* Reworked environments definition files, environments should build correctly.

### Version 1.3.2

* Added a very basic test script. This is meant to test the installation - not provide unit testing
* Added an example of expected output
* Fixed Snakemake version
* Fixed summary report

### Version 1.3.1

* Workflow will no longer crash on blank samples
* Now reports the proportion of reads discarded during primer trimming

### Version 1.3.0

#### Incompatible changes

* Now requires user to provide a fasta file with primer sequences

#### New features

* Taxonomic reports now include a 'disambiguation' field summarizing the different blast hit for each cluster
* Primers will now be trimmed for the reads before quality trimming. It is possible to trim primers on both ends

#### Minor changes

* Performance fix for the display of large tables in the html report

### Version 1.2.2:

* Updated BLAST+ and Fastp to the latest version
* Report now includes links to BLast reports
* Blast report now includes number of mismatch, gaps and alignment length
* Added the --skip_adapter_trimming option to disable adapter trimming in fastp (only recommended for artificial dataset)

### Version 1.2.1:

* taxidTools is now a submodule
* Cloning the repository should now be done with '--recurse-submodules'
* taxidtools updated to version 2
* Adapted scripts to the new version of taxidtools
* Changed BLAST database masking to not be silent about taxids missing from the Taxdump definition files

### Version 1.2.0:

#### New features

* Added the option to filter the BLAST search by taxid

#### Fixes

* Fixed a performance issue for BLAST filtering
* Better error handling for LCA determination
* Snakemake logging has been moved to the logs folder

### Version 1.1.0:

#### New features

* Added primer trimming option (experimental)
* Added subspecies to taxonomy levels
* Added a helper script to fetch the BLAST nt database

#### Fixes

* Fixed Krona broken link in report
* Fixed BLAST filtering for floating point values of bitscores
* Fixed crash upon absence of BLAST hits
* Fixed BLAST database version reporting

### Version 1.0.0:

* initial release
