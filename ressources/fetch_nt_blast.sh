#!/usr/bin/env bash
#set -e
#set -u
#set -o pipefail

# Fetch the pre-built BLAST nt database and taxdump files
# Author: G. Denay, gregoire.denay@cvua-rrw.de

VERSION=2.0

# URLs -------------------------------------------------------------

BLAST="https://ftp.ncbi.nlm.nih.gov/blast/db/"
TAXDUMP="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/"

# Arguments parsing ------------------------------------------------

USAGE="Usage: $0 -d DIRECTORY [-v] [-h]"

while getopts :d:chvt opt
do
  case $opt in
  d ) directory=$OPTARG
      ;;
  h ) help=true
      ;;
  v ) version=true
      ;;
  : ) echo "Missing option argument for -$OPTARG" >&2
      echo $USAGE >&2
      exit 1
      ;;
  '?' ) echo "$0: invalid option -$OPTARG" >&2
      echo $USAGE >&2
      exit 1
      ;;
  esac
done

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
  echo "  -d: output directory for the database"
  echo "  -v: Print version and exit"
  echo "  -h: Print this help and exit"
fi  

## Check if directory is set
if [ -v directory ]; then
  # make sure the directory exists"
  mkdir -p "$directory"
else
  # set current dir als target
  directory="$PWD"
fi

echo "Database will be created to $directory"
cd "$directory"

# Main script ------------------------------------------------------------

# Get Readme
wget --quiet --tries 3 -O README.html ${BLAST}README

# Get directory listing in html format
echo "Retrieving remote directory ${BLAST}"
if curl ${BLAST} --fail --silent > ftp_dir.html; then :
else
  echo "ERROR: URL does not exist: ${BLAST}"
  exit 1
fi

# Extract checksums and fasta links
if [ -f checksums.links ]; then
  rm checksums.links
fi
if [ -f parts.links ]; then
  rm parts.links
fi

echo "Extracting links"
grep -E "\"nt\.[0-9]+\.tar\.gz\"" ftp_dir.html \
  | cut -d'"' -f2 \
  | awk  -v blast=${BLAST} '{print blast$0}' \
  | sort -d \
  >> parts.links
grep -E "\"nt\.[0-9]+\.tar\.gz\.md5\"" ftp_dir.html \
  | cut -d'"' -f2 \
  | awk  -v blast=${BLAST} '{print blast$0}' \
  | sort -d \
  >> checksums.links

echo "Fetching cheksums"
wget --quiet --tries 3 -i checksums.links -O checksums.md5

# Fetching parts one by one
echo "Fetching remote database"
while IFS= read -r part; do
  # check if file exist
  echo $(basename ${part})
  if [ -f $(basename ${part}) ]; then
    echo "$part already exists"
    # check md5
    md5sum -c --ignore-missing checksums.md5
    # md5 exit status should be 0 if everything is ok
    if [ $? -ne 0 ]; then
      # Re download and check md5
      echo "Checksum invalid, redownlaoding $part"
      wget --quiet --tries 3 $part
      md5sum -c --ignore-missing checksums.md5
    fi
  else
    # download and check md5
    echo "Downloading $part"
    wget --quiet --tries 3 $part
    md5sum -c --ignore-missing checksums.md5
  fi
  
  # Check md5 status and exit on error
  if [ $? -ne 0 ]; then
    echo "ERROR: md5 checksum invalid for file $part"
    exit 1
  fi
  
  # Unpack and clean up
  tar -xzf $part \
    && rm $part
done < parts.links

## Taxdump
echo "Fetching Taxdump"
wget --quiet --tries 3 -nc ${TAXDUMP}new_taxdump.tar.gz
wget --quiet --tries 3 -O txd_checksum.md5 ${TAXDUMP}new_taxdump.tar.gz.md5

md5sum -c txd_checksum.md5
if [ $? -ne 0 ]; then
   "ERROR: md5 checksum invalid for file taxdump"
fi

echo "done"
