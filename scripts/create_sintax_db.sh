#!/usr/bin/env bash

# First create BLast DB!

fasta=../db/blast/mitochondrion.LSU.faa
acc2taxid=../db/blast/acc2taxid.tsv

acc=../db/sintax/acc.txt

# filter acc2taxid to spare ~4GB memory
grep -E "^>" $fasta | cut -d " " -f 1 | tr -d ">" > $acc
grep -f $acc $acc2taxid > ../db/sintax/acc2taxid.filtered.txt

python create_sintax_db.py