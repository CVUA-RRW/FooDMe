# FooDMe - A pipeline for Food DNA Metabarcoding

FooDMe is a pipeline for taxonomic assignement of targeted sequencing reads (DNA Metabarcoding). 
It was designed with 16S amplicon sequencing of Food samples (Animal and birds metabarcoding) in mind but could be applied to
 other datasets. 
FooDMe will process demultiplexed Illumina sequencing reads to:

* Cut low quality 3' ends and apply basic quality filtering
* Determine OTU or ASV in a sample-wise fashion and apply some quality filtering at the read and cluster levels
* BLAST sequences in a user-provided database
* Determine a taxonomic consensus based on sequence similarity
* Output quality reports and results

## Getting started 

### Prerequisites

FooDMe runs in a UNIX environment with BASH (tested on Debian GNU/Linux 10 (buster)) and requires conda and an internet connection (at 
least for the first run).

### Installing

Start by getting a copy of this repository on your system, either by downloading and unpacking and archive, or using 'git clone'.

Set up a conda environment containing snakemake, python and the pandas library and activate it:

```bash
conda create --name foodme -c bioconda -c anaconda snakemake pandas
```

### Getting the databases

FooDMe requires several databases to run, all are available from the NCBI ftp servers:

* [taxdump](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz)
* [taxdb](https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz)
* Any nucleotide collection you want to use, this needs to be a searchable BLAST database with taxonomy information. For this
 you can get your own local copy of the BLAST nucleotide database, or build a local database from a subset of sequences. Check
 the [BLAST documentation](https://www.ncbi.nlm.nih.gov/books/NBK279688/) to know how to do this.
 
This distribution provides a utility script to download and set up the database for **mitochondria RefSeq sequences**. Running
this script will require the BLAST command line application:

```bash
conda create --name blast -c bioconda blast
conda activate blast
bash /path/to/FooDMe/ressources/create_RefSeq_blastdb.sh -d /path/to/database -t -c
```

This will download all nescessary files to the '/path/to/database/' folder.

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

Below is a minimal exemple for using the python wrapper, to get the full list of arguments use `foodme.py -h`.

```bash
DATABASES=/path/to/database/folder
python /path/to/FooDMe/foodme.py -l /path/to/sample_sheet.tsv \
	-d /path/to/working/dir \
	--taxdump ${DATABASES} \
	--taxdb ${DATABASES} \
	--blastdb ${DATABASES}/my_blast_db.fasta
```

Below is a minimal exemple for calling snakemake directly. Consult [snakemake's documentation](https://snakemake.readthedocs.io/en/stable/) for more details.

```bash 
snakemake -s /path/to/FooDMe/Snakefile --config path/to/config.yaml --use-conda
```
## Workflow details

### Reads pre-processing

As first analysis step, the reads will be pre-processed for quality trimming on the 3' end and adapter trimming.
If you want to trim the primer sequences, you can do so by indicating the forward and reverse primer length with the 
`--fastp_prune1` and `--fastp_prune2` arguments (experimental).

### Clustering methods

By default FooDMe uses identity clustering as implemented by VSearch. To change this behaviour uses the `--denoise` option in 
the python wrapper. This will results in FooDMe using the amplicon denoising procedure implemented in DADA2.

#### Identity clustering

For identity clustering, the paired reads will first be merged into pseudo reads. Quality filtered pseudo-reads will then be 
used to form Operational Taxonomic Units (OTU) based on the identity threshold specified by `--cluster_id`. OTUs containing 
less the the minimal amount of reads specified by `--cluster_minsize` will be discarded.

##### Amplicon denoising

For amplicon denoising, reads will first be quality filtered and error-corrected. Corrected reads will then be merged and 
Amplicon Sequence Variants (ASV) will be determined. As ASV infer the real composition of the sample, the number of ASV 
should be much closer to the expected number of different sequences in the sample than that of OTUs.

#### Chimera filtering

FooDMe will try to determine chimeric sequences after clustering and these will be discarded. To remove this behaviour use 
the `--skip_chimera`flag.

### BLAST filtering

You can fine tune the BLAST procedure by specifying a minimal e-value, identity, and coverage of the BLAST search. I however 
recommend to keep these parameters relatively permissive and to filter the BLAST results using the bit-score difference. 

Using the bit-score gives you a stable metric that is database-independent. 

### Taxonomic consensus determination

Consensus determination will return the last common node of all retrieved BLAST hits for a sequence. You should expect most 
sequences to be determined at the species or genus level. 

## Credits

FooDMe is built with [Snakemake](https://snakemake.readthedocs.io/en/stable/) and uses the following tools:

* [Fastp](https://github.com/OpenGene/fastp)
* [VSearch](https://github.com/torognes/vsearch) 
* [DADA2](https://benjjneb.github.io/dada2/)
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) 
* [Krona](https://github.com/marbl/Krona)
* [AQUAMIS' create_sampleSheet script](https://gitlab.com/bfr_bioinformatics/AQUAMIS)

## Versionning

Stable version will be given a tag in the form 1.0.0

## Contributing

For new features or to report bugs please submit issues directly on the online repository.

## License

This project is licensed under a BSD 3-Clauses License, see the LICENSE file for details.

## Author

For questions about the pipeline, problems, suggestions or requests, feel free to contact:

Grégoire Denay, Chemisches- und Veterinär-Untersuchungsamt Rhein-Ruhr-Wupper 

<gregoire.denay@cvua-rrw.de>
