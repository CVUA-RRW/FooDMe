#!/usr/bin/env bash
#set -e
#set -u
#set -o pipefail

# Build a BLAST databse containing RefSeq Mitochondrion genomes
# Author G. Denay, gregoire.denay@cvua-rrw.de

VERSION=1.0

# Take arguments
USAGE="Usage: $0 -d DIRECTORY [-v] [-h]"

while getopts :d:hvt opt
do
	case $opt in
	d	)	directory=$OPTARG
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
if [ -z $directory ]
then 
	directory=$PWD
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
	echo "	-v: Print version and exit"
	echo "	-h: Print this help and exit"
fi	

# test mode
if [[ $test == true ]]
then
	exit 0
fi

# Get source files
# taxdb
wget -P ./source/ https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz \
	&& tar -xzvf ./source/taxdb.tar.gz \
	&& rm ./source/taxdb.tar.gz
	
# fasta files
wget -P ./source/ https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.1.genomic.fna.gz \
	&& gunzip -c ./source/mitochondrion.1.1.genomic.fna.gz > ./source/mitochondrion.1.1.genomic.fna
wget -P ./source/ https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.2.1.genomic.fna.gz \
	&& gunzip -c ./source/mitochondrion.2.1.genomic.fna.gz > ./mitochondrion.2.1.genomic.fna
cat ./source/mitochondrion.1.1.genomic.fna ./source/mitochondrion.2.1.genomic.fna > mitochondrion.genomic.fna \
	&& rm ./source/mitochondrion.1.* ./source/mitochondrion.2.*
	
# taxid file
wget -P ./source/ https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz \
	&& gunzip -c ./source/nucl_gb.accession2taxid.gz > ./source/nucl_gb.accession2taxid \
	&& cut -f 1,3 ./source/nucl_gb.accession2taxid | tail -n +2 > genbank2taxid \
	&& rm ./source/nucl_gb*
	
# Make db
makeblastdb -in mitochondrion.genomic.fna -parse_seqids -blastdb_version 5 -taxid_map genbank2taxid -dbtype nucl