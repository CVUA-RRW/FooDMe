# FooDMe - A pipeline for Food DNA Metabarcoding

## Description

FooDMe is a pipeline for taxonomic assignement of targeted sequencing reads (DNA Metabarcoding). 
It is meant for food authenticity analysis but would most likely work for other applications too.
It works on raw paired-end Illumina reads, it will trim them, assemble read-pairs and perform a quality filtering, merge samples to filter out chimeras and find OTUs,
map the assemble reads to the OTUs, and finally, perform taxonomic assignment using a nucleotide database.

FooDMe is built on Snakemeke and uses the following tools:

* Fastp
* VSearch 
* BLAST+ 

## Source

FooDMe is not yet distributed through a centralized hub. It is however Git-versionned.
For the latest version contact, see "Contact".

## Installation

Get the latest version and upack it in the Repo directory:

```bash
REPO=path/to/repo
cd $REPO
tar -xzvf FooDme.tar
```

### Conda environment

FooDMe requires conda to manage environments of all dependencies. You can install conda from [here](https://docs.conda.io/en/latest/miniconda.html) or chose another distribution.
Then you need to setup an environment to lauch the pipeline:

```bash
conda create -n foodme snakemake biopython
```

Additional environments should be installed during the first execution of the pipeline.

### BLAST database

Taxonomic assignment of reads requires sequence alignment against a nucleotide database. Such a database can be directly downloaded from the [NCBI server](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
or built using the provided helper script. 

```bash
bash ${REPO}/FooDMe/scripts/create_blast_db.sh 
```

By default the script will download the mitochondria RefSeq files from the NCBI server and create a database with the Vertebrata large ribosomal subunit (16S) sequences in the `db/blast` subfolder.

The database location can be modified by providing it as an argument to the script:

```bash
bash ${REPO}/FooDMe/scripts/create_blast_db.sh path/to/database
```

The content of the database can be tweaked by modifying the python call to the `GenBank_to_FASTA_parser.py` script:

```bash
python GenBank_to_FASTA_Parser.py -i [input gbff file] -o [output fasta name] -t [taxa filter] -f [feature 1] -f [feature 2]
```

The `taxa filter` and `feature` arguments must be provided as strings matching the gbff nomenclature. An arbitrary number of features can be provided. 

## Usage

### Execution

To run FooDMe, first activate the conda environment and call the snakemake file:

```bash
conda activate foodme
snakemake -s ${REPO}/FoodMe/Snakefile  --conda-prefix /path/to/conda/envs --configfile /path/to/config.yaml --use-conda
```

### Sample sheet 

FooDMe execution requires a sample sheet to be provided. This is a tab separated table referencing mate pair files for each samples.
Such a table can be automatically generated with the script `create_sampleSheet.sh`, originally developped by the BfR Study Center for Genome Sequencing and Analysis :

```
bash ${REPO}/scripts/create_sampleSheet.sh -f path/to/reads
```

This script includes a range of options for dealing with no-Illumina file name formatting. 

### Configuration

Pipeline parameters can be modified in the provided `config.yaml` file.
Here is a breakdown of the parameters:

* `workdir`: Path to the folder where the files and reports will be generated
* `samples`: path to the sample list 
* Ressource allocation:
	- `threads`: Number of threads to allocate per sample
	- `cores`: Number of cores to allocate for the jobs involving pooled samples (Pipeline bottleneck)
* fastp:
	- `length_required`: Minimal read length
	- `qualified_quality_phred`: Minimal quality score
* Read filtering:
	- `min_length`: Minimal length of assembled reads
	- `max_length`: Maximal length of assembled reads
	- `max_expected_errors`: Maximal allowed number of expected errors per assembled reads 
* Clustering:
	- `cluster_identity`: Minimal identity value for OTU clustering (0 to 1)
	- `chimera_DB`: Path to database for chimera filtering. Ideally a high-quality collection of parent sequences.
* BLAST:
	- `blast_DB`: Path to the Blast database. See **Installation**.
	- `taxdb`: Path to the *folder* containing the ncbi taxonomy database files (blastdb.btd and blastdb.bti)
	- `e_value`: Maximal e-value allowed for blast hits
	- `perc_identity`: Minimal identity allowed for balst hits (0 to 100)
	- `qcov`: Minimal query coverage allowed for blast hits (0 to 100)
	- `bit_score_diff`: Maximal bit-score difference allowed to taxonomic consensus determination
* Taxonomic assignement:
	- `names_dmp`: Path to the names.dmp taxonomy file
	- `nodes_dmp`: Path to the nodes.dmp taxonomy file 

## Reports

FooDMe generates .tsv files with quality summaries of the different steps of the analysis. These files are located in the reports subfolder of the working directory.
The final result consists of a composition table summarizing the abundance of each taxa in the different samples. These can be found in the reports/results subfolder.

## Workflow details

### Read trimming, merging, and quality filtering

Raw reads are pre-processed with fastp to remove adapters and mask low quality regions. 
Mate pairs are then merged and the assembled reads are filtered for the their length and overall quality.

### Clustering and chimera filtering

Reads from all samples are pooled and dereplicated (exact match). Reads are then clustered in centroids using VSearch's cluster\_size option. 
Chimeras are then detected and removed using *de novo* and database-based approaches.

### Taxonomic assignment

Dereplicated reads are mapped to the OTUs using VSearch global alignment otpion usearch_global. 
Each OTU is blasted against a nucleotide reference database to determine the Taxa of origin. 
Hits are filtered by bit-score difference to the bast hit to remove lower similarity hits. 
A taxonomic consensus is then determined as the lowest common node of all hits.

## Credits

This pipeline includes thrird-party source-code, available under BSD-3 license:

 * The create_sampleSheet.sh script as well as code snippets were written by Carlus Deneke <https://gitlab.com/bfr_bioinformatics/AQUAMIS>
 * Taxonomic lineage extraction uses the ncbi\_taxdump\_utils.py module from Titus Brown, <https://github.com/dib-lab/2018-ncbi-lineages>

## Contact

For questions about the pipeline or problems, feel free to contact:
Grégoire Denay, Chemisches- und Veterinär-Untersuchungsamt Rhein-Ruhr-Wupper, <gregoire.denay@cvua-rrw.de>
