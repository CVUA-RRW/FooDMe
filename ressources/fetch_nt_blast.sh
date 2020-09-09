#!/usr/bin/env bash
#set -e
#set -u
#set -o pipefail

# Fetch the pre-built BLAST nt database and taxdump files
# Author: G. Denay, gregoire.denay@cvua-rrw.de

VERSION=1.0

# URLs -------------------------------------------------------------

BLAST="https://ftp.ncbi.nlm.nih.gov/blast/db/"
TAXDUMP="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/"

# Arguments parsing ------------------------------------------------

USAGE="Usage: $0 -d DIRECTORY [-v] [-h]"

while getopts :d:chvt opt
do
	case $opt in
	d	)	directory=$OPTARG
			;;
	h 	) 	help=true
			;;
	v	)	version=true
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

## if no directory specified use current directory
if [ -z "$directory" ]
then 
	directory="$PWD"
fi

## version
if [[ $version == true ]]
then 
	echo "fetch_nt_blast.sh version: $VERSION"
fi

# help
if [[ $help == true ]]
then 
	echo "fetch_nt_blast.sh (version: $VERSION)"
	echo "Fetch the pre-built BLAST nt database and taxdump files"
	echo
	echo $USAGE
	echo
	echo "Options:"
	echo "	-d: output directory for the database"
	echo "	-v: Print version and exit"
	echo "	-h: Print this help and exit"
fi	

# Main script ------------------------------------------------------------

echo "Fetching database to ${PWD}"
cd ${PWD}

# Get directory listing in html format

echo "Retrieving remote directory ${BLAST}"
curl ${BLAST} > ftp_dir.html 

echo "Fetching remote database"

## BLAST DB
grep -E "\"nt\.[0-9]+\.tar\.gz\"" ftp_dir.html \
	| cut -d'"' -f2 \
	| sed "s;^;wget --tries 5 -nc ${BLAST};" \
	| sh -x

grep -E "\"nt\.[0-9]+\.tar\.gz\.md5\"" ftp_dir.html \
	| cut -d'"' -f2 \
	| sed "s;^;wget --tries 5 -O ->> checksum.md5 ${BLAST};" \
	| sh -x

wget --tries 5 -nc -O README.html ${BLAST}README

## Taxdump
wget --tries 5 -nc ${TAXDUMP}new_taxdump.tar.gz
wget --tries 5 -O ->> checksum.md5 ${TAXDUMP}new_taxdump.tar.gz.md5

# Check md5 
echo "Checking md5 sums"

for f in $(cut -d" " -f3 checksum.md5); do
	if [[ $(md5sum $f) == $(grep $f checksum.md5) ]]; then
		echo "$f  ok"
	else
		echo "ERROR: md5 checksum of $f has unexpected value"
		exit 1
	fi
done

# Unpacking
echo "Unpacking archives"
for f in *.tar.gz; do
	tar -xzvf $f \
	&& rm $f
done

# clean up 
rm checksum.md5
rm ftp_dir.html

echo "done"
