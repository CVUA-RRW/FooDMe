#! bin/bash

cd ${BLASTDB}

mkdir ncbi_dump
cd ncbi_dump

# Get datasets
wget https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.genomic.gbff.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.2.genomic.gbff.gz

cd ${BLASTDB}

# Extract taxid database
tar -xzvf ncbi_dump/taxdb.tar.gz

# Parse gbff to extract 16S sequences from Metazoans
python GenBank_to_FASTA_Parser.py

# Remove duplicates (Multtiple entries of the 16s rRNA for the same organism)
sed '/^>/s/$/@/' < ncbi_dump/mitochondrion.16S.metazoan.faa |\
sed 's/^>/$/' |\
tr -d '\n' | tr "$" "\n" | tr "@" "\t" |\
sort -u -f -k 1,1  |\
sed -e 's/^/>/' |\
tr '\t' '\n/' |\
tail -n +2 > mitochondrion.16S.metazoan.faa

# Clean Accession to taxid table
cut -f 1,3 ncbi_dump/nucl_gb.accession2taxid.txt | tail -n +2 > ncbi_dump/acc2taxid.tsv

# Make BLAST database
makeblastdb -in mitochondrion.16S.metazoan.faa -parse_seqids -blastdb_version 5 -taxid_map ncbi_dump/acc2taxid.tsv -title "Metazoa_16S" -dbtype nucl


