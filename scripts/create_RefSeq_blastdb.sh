#!/usr/bin/env bash
#set -e
#set -u
#set -o pipefail

# Build a BLAST databse containing RefSeq Mitochondrion genomes
# Author: G. Denay, gregoire.denay@cvua-rrw.de

VERSION=1.0

# Take arguments
USAGE="Usage: $0 -d DIRECTORY [-t] [-c] [-v] [-h]"

while getopts :d:chvt opt
do
	case $opt in
	d	)	directory=$OPTARG
			;;
	t	)	taxdmp=true
			;;
	c	)	clean=true
			;;
	h 	) 	help=true
			;;
	v	)	version=true
			;;
	t	)	test=true
			;;
	:	) 	echo "Missing option argument for -$OPTARG" >&2
			echo $USAGE >&2
			exit 1
			;;
	'?'	)	echo "$0: invalid option -$OPTARG" >&2
			echo $USAGE >&2
			exit 1
			;;
	esac
done

# if no directory specified use current directory
if [ -z "$directory" ]
then 
	directory="$PWD"
fi

# version
if [[ $version == true ]]
then 
	echo "create_RefSeq_blastdb.sh version: $VERSION"
fi

# help
if [[ $help == true ]]
then 
	echo "create_RefSeq_blastdb.sh (version: $VERSION)"
	echo "Localy create a BLAST database of RefSeq mitochondrion sequences"
	echo
	echo $USAGE
	echo
	echo "Options:"
	echo "	-d: output directory for the database"
	echo "	-t: include the download of taxdump files"
	echo " 	-c: clean-up source files before exiting"
	echo "	-v: Print version and exit"
	echo "	-h: Print this help and exit"
fi	

# Create directory
if [ ! -d "$directory" ] 
then 
	mkdir -p "$directory"
fi

cd "$directory"
mkdir source

# test mode
if [[ $test == true ]]
then
	echo "$PWD"
	exit 0
fi

# Get source files
# taxdb
wget -S -P ./source/ https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz \
	&& tar -xzvf ./source/taxdb.tar.gz 
	
# fasta files
wget -S -P ./source/ https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.1.genomic.fna.gz \
	&& gunzip -c ./source/mitochondrion.1.1.genomic.fna.gz > ./source/mitochondrion.1.1.genomic.fna
	
wget -S -P ./source/ https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.2.1.genomic.fna.gz \
	&& gunzip -c ./source/mitochondrion.2.1.genomic.fna.gz > ./source/mitochondrion.2.1.genomic.fna
	
cat ./source/mitochondrion.1.1.genomic.fna ./source/mitochondrion.2.1.genomic.fna > mitochondrion.genomic.fna 

wget -S https://ftp.ncbi.nih.gov/refseq/release/RELEASE_NUMBER \
	&& echo "Downloaded RefSeq Mitochondrion genomic sequences release number: " $(cat RELEASE_NUMBER)

# taxid file
wget -S -P ./source/ https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz \
	&& gunzip -c ./source/nucl_gb.accession2taxid.gz > ./source/nucl_gb.accession2taxid \
	&& cut -f 1,3 ./source/nucl_gb.accession2taxid | tail -n +2 > genbank2taxid 


# Make db
makeblastdb -in mitochondrion.genomic.fna -parse_seqids -blastdb_version 5 -taxid_map genbank2taxid -dbtype nucl

# Taxdump
if [[ $taxdmp == true ]]
then
	wget -S -P ./source/ https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz \
	&& tar -xzvf ./source/new_taxdump.tar.gz
fi

# Clean up
if [[ $clean == true ]]
then
	echo "Cleaning up sources"
	rm -r ./source
fi
