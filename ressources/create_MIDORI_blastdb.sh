# Get MIDORI Fasta and reformat header
wget -S -P ./source/ https://www.reference-midori.info/forceDownload.php?fName=download/Databases/GenBank259_2023-12-17/BLAST/uniq/fasta/MIDORI2_UNIQ_NUC_GB259_lrRNA_BLAST.fasta.zip \
    && gunzip -c ./source/MIDORI2_UNIQ_NUC_GB259_lrRNA_BLAST.fasta.zip > ./source/MIDORI2_UNIQ_NUC_GB259_lrRNA_BLAST.fasta \
    && cut -d '#' -f1 ./source/MIDORI2_UNIQ_NUC_GB259_lrRNA_BLAST.fasta \
        | tr -d '<' \
        | sed 's/^>/@/' \
        | tr -d '>' \
        | tr '@' '>' \
        | cut -d ',' -f1,2 \
        | tr ',' '_' \
        > MIDORI2_UNIQ_NUC_GB259_lrRNA_BLAST.fasta

# taxdb
wget -S -P ./source/ https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz \
	&& tar -xzvf ./source/taxdb.tar.gz 

# taxdump
wget -S -P ./source/ https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz \
    && tar -xzvf ./source/new_taxdump.tar.gz

# taxid file
# Get taxid from fasta headers
paste \
    <(cat ./source/MIDORI2_UNIQ_NUC_GB259_lrRNA_BLAST.fasta | grep '^>' | cut -d '#' -f1 | tr -d '<' | tr -d '>' | cut -d ',' -f1,2 | tr ',' '_') \
    <(cat ./source/MIDORI2_UNIQ_NUC_GB259_lrRNA_BLAST.fasta | grep '^>' | rev | cut -d '_' -f1 | rev) \
    > ids2taxid

# Make db
makeblastdb -in MIDORI2_UNIQ_NUC_GB259_lrRNA_BLAST.fasta -parse_seqids -blastdb_version 5 -taxid_map ids2taxid -dbtype nucl
