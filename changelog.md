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
