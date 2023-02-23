# Nucleotide database

In order to taxonomically assign sequences, a certain amount of reference data is required.
Below are the instructions to retrieve standard databases and how to create custom ones.

## BLAST database

The sequence comparison step as implemented uses the BLAST command line tools. This requires that the 
nucleotide database is indexed and formatted in a BLAST-specific way.

### Preformatted databases

A large pre-formatted nucleotide sequence database is freely available from the NCBI
servers. The BLAST NT database contains a collection of sequences from different sources
and can be downloaded directly from the [NCBI servers](https://ftp.ncbi.nlm.nih.gov/blast/db/). 

A script to fetch the BLAST NT database and additional required taxonomic definitions 
is available with FooDMe:

```bash
cd ~
mkdir blast_nt
bash ~/FooDMe/ressources/fetch_nt_blast.sh -d blast-nt
```

Running the above commands will create the blast_nt directory and retrieve all the 
nescessary files from the NCBI servers.

!!! warning 

    Downloading the BLAST NT database will require a large chunk of available memory
    and take several hours. 

### Custom databases

A collection of non-redundant reference sequences (RefSeq) is also available from the NCBI servers.
The collection is not yet in a BLAST format but we provide a script to retrieve and format it.

Running the commands below will create a folder called refseq, download the RefSeq collection and 
format it in a blast format. This requires to create a conda environment containing the blast tools.

```bash
mamba create -n blast -c conda-forge -c bioconda blast
conda activate blast
cd ~
mkdir refseq
bash ~/FooDMe/ressources/create_RefSeq_blastdb.sh -d refseq -t
```

Additionally it is possible to format any sequence collection in a BLAST compliant format.
A User Guide therefore is available from the [BLAST documentation](https://www.ncbi.nlm.nih.gov/books/NBK569841/).

## Additonal files

In addition to the nucleotide collection, several taxonomy definition files are nescessary.
These are available from the NCBI servers and can be used with other sources.

!!! warning 

    The links below are provided as information and might not link to the most recent 
    version of the files. 

!!! note 

    If using the scripts above to retrieve either the NT or RefSeq database, the files 
    below should already be included in your local database.

### Taxonomic information

The `taxdb` files are nescessary to link taxonomic identifier (taxid) and human-readable 
informations. These can be downloaded here:

> [https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz](https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz)

### Taxonomic classification

The `taxdump` files contain the taxonomy hierarchical information and are nescessary 
to determine the degree of relationship between different taxa. They can be downloaded here:

> [https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)

