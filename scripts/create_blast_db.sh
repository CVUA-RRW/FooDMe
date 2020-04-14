#!/usr/bin/env bash

if [ $# -eq 0 ]
then
	SCRIPT=$(dirname $(readlink -f "$0"))
	BLASTDB=$SCRIPT/../db/blast
else
	BLASTDB=$1
fi

mkdir -p ${BLASTDB} && cd ${BLASTDB}

# Get datasets
wget https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz && tar -xzvf taxdb.tar.gz && rm taxdb.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz && gunzip nucl_gb.accession2taxid.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.genomic.gbff.gz && gunzip mitochondrion.1.genomic.gbff.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.2.genomic.gbff.gz && gunzip mitochondrion.2.genomic.gbff.gz

cat mitochondrion.1.genomic.gbff mitochondrion.2.genomic.gbff > mitochondrion.merged.genomic.gbff && rm mitochondrion.1.genomic.gbff mitochondrion.2.genomic.gbff

# Parse gbff to extract 16S sequences from Metazoans
python GenBank_to_FASTA_Parser.py -i mitochondrion.merged.genomic.gbff -o mitochondrion.LSU.raw.faa -t Vertebrata -f 16S -f l-rRNA

# Remove duplicate entries (Same Accession and sequence)
sed '/^>/s/$/@/' < mitochondrion.LSU.raw.faa |\
sed 's/^>/#/' |\
tr -d '\n' | tr "#" "\n" |
sort -u -f -k 1,1 |\
sed -e 's/^/>/' |\
tr '@' '\n/' | tail -n +2 > mitochondrion.LSU.faa && rm mitochondrion.LSU.raw.faa

# Clean Accession to taxid table
cut -f 1,3 nucl_gb.accession2taxid | tail -n +2 > acc2taxid.tsv && rm nucl_gb.accession2taxid

# Make BLAST database
makeblastdb -in mitochondrion.LSU.faa -parse_seqids -blastdb_version 5 -taxid_map acc2taxid.tsv -title "Vertebrata_16S" -dbtype nucl

echo "Database created at " ${BLASTDB}