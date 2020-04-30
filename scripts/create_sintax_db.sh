#!/usr/bin/env bash

fasta=../db/blast/mitochondrion.LSU.faa
acc2taxid=../db/blast/acc2taxid.tsv

acc=../db/sintax/acc.txt

# filter and sort acc2taxid
grep -E "^>" $fasta | cut -d " " -f 1 | tr -d ">" > $acc
#grep -f $acc $acc2taxid | sort -k1,1 > ../db/sintax/acc2taxid.filtered.txt
grep -f $acc $acc2taxid > ../db/sintax/acc2taxid.filtered.txt

# sort fasta
# sed '/^>/s/$/@/' < $fasta |\
# sed 's/^>/#/' |\
# tr -d '\n' | tr "#" "\n" |
# sort -k 1,1 |\
# sed -e 's/^/>/' > ../db/sintax/mitochondrion.LSU.sorted.faa 

python create_sintax_db.py