![CI](https://github.com/CVUA-RRW/FooDMe/workflows/CI/badge.svg?branch=master)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/CVUA-RRW/FooDMe)

# FooDMe - A pipeline for Food DNA Metabarcoding

FooDMe is a snakemake workflow for taxonomic assignement of targeted sequencing reads (DNA Metabarcoding). 
It was designed with 16S amplicon sequencing of food samples (mammals and birds metabarcoding) in mind but could be applied to
 other datasets. 
FooDMe will process demultiplexed Illumina sequencing reads to:

* Trim primers
* Cut low quality 3' ends and apply basic quality filtering
* Cluster sequences in a sample-wise fashion and apply some quality filtering at the read and cluster levels
* BLAST sequences in a user-provided database
* Determine a taxonomic consensus based on sequence similarity
* Output quality reports and results

## Getting started 

### Prerequisites

FooDMe runs in a UNIX environment with BASH (tested on Debian GNU/Linux 10 (buster), may work on MacOS too) and requires conda and an internet 
connection (at least for the first run).

### Installing

Start by getting a copy of this repository on your system, either by downloading and unpacking the archive, 
or using 'git clone':

```bash
cd path/to/repo/
git clone --recurse-submodules https://github.com/CVUA-RRW/FooDMe.git
```

Set up a conda environment containing snakemake, python and the pandas library and activate it:

```bash
conda create --name foodme -c bioconda -c anaconda snakemake pandas
conda activate foodme
```

### Testing the installation

You can run an installation test with the data provided in the `tests` folder:

```bash
conda activate foodme
cd path/to/repo
snakemake --use-conda --cores 1
```

This should produce the same output as the one present in the `foodme-example_output` folder.

Note that if you provide the `--conda-prefix /path/to/conda/envs` argument, the above command will 
set-up the conda environements that you can reuse offline later.

### Getting the databases

FooDMe requires several databases to run, all are available from the NCBI ftp servers:

* [taxdump](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz)
* [taxdb](https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz)
* Any nucleotide collection you want to use, this needs to be a searchable BLAST database with taxonomy information. For this
 you can build a local database from a subset of sequences, for exemple from the BOLD database. Check the
 [BLAST documentation](https://www.ncbi.nlm.nih.gov/books/NBK279688/) to know how to do this.

#### Installing the RefSeq mitochondria database

This distribution provides a utility script to download and set up the database for **mitochondria RefSeq sequences**. Running
this script will require the BLAST command line application:

```bash
conda create --name blast -c bioconda blast
conda activate blast
bash /path/to/FooDMe/ressources/create_RefSeq_blastdb.sh -d /path/to/database -t -c
```

This will download all nescessary files to the '/path/to/database/' folder. Use the '-t' flag to include the download of the 
Taxdump files.

#### Fetching the pre-formatted BLAST nt database

We also provide a utility script, 'fetch_nt_blast.sh' to download the pre-formatted BLAST nt database.
Note that you will need ~110 Gb of available disk space.

```bash
bash /path/to/FooDme/ressources/fetch_nt_blast.sh -d /path/to/database
```

### Creating a sample sheet

FooDMe requires a tabular file linking sample names to forward and reverse read files.
This can be produced using the 'create_sampleSheet.sh' script from Carlus Deneke (BfR).

Consult the help with:

```bash
bash /path/to/FooDMe/ressources/create_sampleSheet.sh -h
```

### Running FooDMe

You can run FooDMe either by using the python wrapper or by calling directly snakemake. The python wrapper will each time 
generate a config file containing the run's parameters. 

Calling directly snakemake will allow you to reuse previously generated config files. This is especially useful to routinely 
run the pipeline with fixed parameters. 

#### Using the python wrapper

```
usage: FooDMe [-h] [-v] -l SAMPLE_LIST -d WORKING_DIRECTORY [--forceall] [-n]
              [-T THREADS] [-t THREADS_SAMPLE] [-c CONDAPREFIX] [-s SNAKEFILE]
              [--keep_temp] [--fastp_length FASTP_LENGTH]
              [--fastp_min_phred FASTP_MIN_PHRED]
              [--fastp_window FASTP_WINDOW] [--fastp_meanq FASTP_MEANQ]
              --primers_fasta PRIMERS_FASTA [--trim_3end]
              [--primer_error_rate PRIMER_ERROR_RATE]
              [--merge_minlength MERGE_MINLENGTH]
              [--merge_maxlength MERGE_MAXLENGTH] [--merge_maxee MERGE_MAXEE]
              [--merge_maxns MERGE_MAXNS] [--denoise]
              [--cluster_id CLUSTER_ID] [--cluster_minsize CLUSTER_MINSIZE]
              [--skip_chimera] [--taxdump TAXDUMP] [--nodes_dmp NODES_DMP]
              [--rankedlineage_dmp RANKEDLINEAGE_DMP] --blastdb BLASTDB
              --taxdb TAXDB [--taxid_filter TAXID_FILTER]
              [--blast_eval BLAST_EVAL] [--blast_id BLAST_ID]
              [--blast_cov BLAST_COV] [--bitscore BITSCORE]

Another pipeline for (Food) DNA metabarcoding

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Print pipeline version and exit

I/O path arguments:
  -l SAMPLE_LIST, --sample_list SAMPLE_LIST
                        Tab-delimited list of samples and paths to read files.
                        Must contain one line of header, each further line
                        contains sample_name, read1_path, read2_path (default:
                        None)
  -d WORKING_DIRECTORY, --working_directory WORKING_DIRECTORY
                        Directory to create output files (default: None)

Snakemake arguments:
  --forceall            Force the recalculation of all files (default: False)
  -n, --dryrun          Dryrun. Create config file and calculate the DAG but
                        do not execute anything (default: False)
  -T THREADS, --threads THREADS
                        Maximum number of threads to use (default: 8)
  -t THREADS_SAMPLE, --threads_sample THREADS_SAMPLE
                        Number of threads to use per concurent job (default:
                        1)
  -c CONDAPREFIX, --condaprefix CONDAPREFIX
                        Location of stored conda environment. Allows snakemake
                        to reuse environments. (default: False)
  -s SNAKEFILE, --snakefile SNAKEFILE
                        Path to the Snkefile in the FOodMe repo (default:
                        /home/debian/NGS/spezies_indev/FooDMe/Snakefile)
  --keep_temp           Keep large fasta and fastq files, mostly for debug
                        purposes (default: False)

Trimming options:
  --fastp_length FASTP_LENGTH
                        Minimum length of input reads (after primer trimming)
                        to keep (default: 50)
  --fastp_min_phred FASTP_MIN_PHRED
                        Minimal quality value per base (default: 20)
  --fastp_window FASTP_WINDOW
                        Size of the sliding window for tail quality trimming
                        (default: 4)
  --fastp_meanq FASTP_MEANQ
                        Minimum mean Phred-score in the sliding window for
                        tail quality trimming (default: 25)
  --primers_fasta PRIMERS_FASTA
                        Fasta file with primers sequences for primer trimming
                        (default: None)
  --trim_3end           Should primers be trimmed on the 3' end of reads? Only
                        relevant if sequencing through the amplicons.
                        (default: False)
  --primer_error_rate PRIMER_ERROR_RATE
                        Maximum error-rate allowed for primer matching
                        (default: 0.1)

Merged reads filtering options:
  --merge_minlength MERGE_MINLENGTH
                        Minimum length merged reads to keep (default: 70)
  --merge_maxlength MERGE_MAXLENGTH
                        Maximum length merged reads to keep (default: 100)
  --merge_maxee MERGE_MAXEE
                        Maximum expected errors in merged reads to keep. Note
                        that for the denoising procedure this filter is
                        applied on the reads BEFORE merging. (default: 2)
  --merge_maxns MERGE_MAXNS
                        Maximum number of 'N' base in merged reads. If using
                        denoising procedure this will be automatically reset
                        to 0 (default: 0)

Clustering options:
  --denoise             Use denoising instead of identity clustering (default:
                        False)
  --cluster_id CLUSTER_ID
                        Minimum identity for clustering sequences in OTUs
                        (between 0 and 1). Will be ignored if using denoising
                        (default: 0.97)
  --cluster_minsize CLUSTER_MINSIZE
                        Minimal size cutoff for OTUs. Will be ignored if using
                        denoising (default: 2)
  --skip_chimera        Skip de novo chimera detection and filtering step
                        (default: False)

Taxonomic assignement files:
  --taxdump TAXDUMP     Path to the taxump folder containing nodes.dmp and
                        rankedlineages.dmp (default: None)
  --nodes_dmp NODES_DMP
                        Path to the nodes.dmp file, needed if --taxdump is
                        omitted (default: None)
  --rankedlineage_dmp RANKEDLINEAGE_DMP
                        Path to the names.dmp file, needed if --taxdump is
                        omitted (default: None)

Options for BLAST search:
  --blastdb BLASTDB     Path to the BLAST database, including database
                        basename but no extension (e.g. '/path/to/db/nt')
                        (default: None)
  --taxdb TAXDB         Path to the BLAST taxonomy database (folder) (default:
                        None)
  --taxid_filter TAXID_FILTER
                        Limit BLAST search to the taxids under the given node
                        (default: None)
  --blast_eval BLAST_EVAL
                        E-value threshold for blast results (default: 1e-10)
  --blast_id BLAST_ID   Minimal identity between the hit and query for blast
                        results (in percent) (default: 90)
  --blast_cov BLAST_COV
                        Minimal proportion of the query covered by a hit for
                        blast results. A mismatch is still counting as
                        covering (in percent) (default: 90)
  --bitscore BITSCORE   Maximum bit-score difference with the best hit for a
                        blast result to be included in the taxonomy consensus
                        detemination (default: 0)
```

Below is a minimal exemple for using the python wrapper:

```bash
conda activate foodme

DATABASES=/path/to/database/folder

python /path/to/FooDMe/foodme.py -l /path/to/sample_sheet.tsv \
    -d /path/to/working/dir \
    --taxdump ${DATABASES} \
    --taxdb ${DATABASES} \
    --blastdb ${DATABASES}/my_blast_db \
    --primers_fasta /path/to/my/primers.fa
```

To use the pipeline offline and/or avoid recreating the conda environements use 
`-c /path/to/conda/envs`

#### Calling snakemake with a configuration file

Below is a minimal exemple for calling snakemake directly. Consult 
[snakemake's documentation](https://snakemake.readthedocs.io/en/stable/) for more details.

```bash 
snakemake -s /path/to/FooDMe/Snakefile --config path/to/config.yaml --use-conda --cores 1
```

To use the pipeline offline and/or avoid recreating the conda environements use 
`--conda-prefix /path/to/conda/envs`

## Workflow details

### Reads pre-processing

As a first analysis step, primers will be trimmed form the reads. By default primers 
are only matched on the 5' end of the reads. In some cases (e.g. Sequencing is 
longer than the amplicon length) one may want to trim primers on the 3' end as well. 
This is possible with the `--trim_3end` flag.
The reads will then be pre-processed for quality trimming on the 3'.

### Clustering methods

By default FooDMe uses identity clustering as implemented by VSearch. To change this behaviour uses the `--denoise` option in 
the python wrapper. This will results in FooDMe using the amplicon denoising procedure implemented in DADA2.

#### Identity clustering

For identity clustering, the paired reads will first be merged into pseudo reads. Quality filtered pseudo-reads will then be 
used to form Operational Taxonomic Units (OTU) based on the identity threshold specified by `--cluster_id`. OTUs containing 
less the the minimal amount of reads specified by `--cluster_minsize` will be discarded.

#### Amplicon denoising

For amplicon denoising, reads will first be quality filtered and error-corrected. Corrected reads will then be merged and 
Amplicon Sequence Variants (ASV) will be determined. As ASV infer the real composition of the sample, the number of ASV 
should be much closer to the expected number of different sequences in the sample than that of OTUs.

### Chimera filtering

FooDMe will try to determine chimeric sequences after clustering and these will be discarded. To remove this behaviour use 
the `--skip_chimera` flag.

### BLAST filtering

You can fine tune the BLAST procedure by specifying a minimal e-value, identity, and coverage of the BLAST search. 
The BLAST results will then be filtered using the bitscore of each hit sequence. The maximum allowed bitscore difference to 
the best hit can be changed with the `--bitscore` option.

It can be advisable to limit the BLAST search to the descendant of a node of interest. You can do this by providing the 
parent node to the `--taxid_filter` option. For example providing `32524` will limit the search to Amniota and `40674` will 
limit the search to Mammals.

### Taxonomic consensus determination

Consensus determination will return the lowest common node of all retrieved BLAST hits for a sequence. You should expect most 
sequences to be determined at the species or genus level. 
Additionaly a summary of the the BLAST hits will be shown in the disambiguation column. 
THis could allow you to refine the consensus determination or identify spurious BLAST results.


## Related tools

These other pipelines may be of interest:

* [MetaSeqSim](https://github.com/CVUA-RRW/MetaSeqSim) can be used to generate synthetic metabarcoding datasets
* [BAnalyzer](https://github.com/CVUA-RRW/BAnalyzer) provides quality control of nucleotide databases, specifically for metabarcoding experiments


## Credits

FooDMe is built with [Snakemake](https://snakemake.readthedocs.io/en/stable/) and uses the following tools:

* [Fastp](https://github.com/OpenGene/fastp)
* [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)
* [VSearch](https://github.com/torognes/vsearch) 
* [DADA2](https://benjjneb.github.io/dada2/)
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) 
* [Krona](https://github.com/marbl/Krona)
* [AQUAMIS' create_sampleSheet script](https://gitlab.com/bfr_bioinformatics/AQUAMIS)


## Contributing

All contributions are welcomed.
For new features or to report bugs please submit issues directly on the online repository.


## License

This project is licensed under a BSD 3-Clauses License, see the LICENSE file for details.


## Author

For questions about the pipeline, problems, suggestions or requests, feel free to contact:

Grégoire Denay, Chemisches- und Veterinär-Untersuchungsamt Rhein-Ruhr-Wupper 

<gregoire.denay@cvua-rrw.de>
