#!/usr/bin/env python3

import ncbi_taxdump_utils
import subprocess
from collections import defaultdict
from Bio import SeqIO

want_taxonomy = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']	

def init_tax(names_dmp, nodes_dmp):
	taxfoo = ncbi_taxdump_utils.NCBI_TaxonomyFoo()  
	taxfoo.load_nodes_dmp(nodes_dmp)
	taxfoo.load_names_dmp(names_dmp)
	return taxfoo

def get_description(taxid, taxfoo):
	lin = taxfoo.get_lineage_as_dict(taxid, want_taxonomy)
	desc = ";tax="
	for rank, letter in zip(want_taxonomy, ["k:", ",p:", ",c:", ",o:", ",f:", ",g:", ",s:"]):
		try:
			desc += letter + str(lin[rank])
		except KeyError:
			continue
	return(desc)

if __name__ == '__main__':

	names_dmp = "../db/taxdump/names.dmp"
	nodes_dmp = "../db/taxdump/nodes.dmp"
	fasta = "../db/blast/mitochondrion.LSU.faa"
	sintax = "../db/sintax/mitochondrion.LSU.sintax.faa"
	acc2taxid = "../db/sintax/acc2taxid.filtered.txt"

	taxfoo= init_tax(names_dmp, nodes_dmp)

	# Load taxid records
	with open(acc2taxid, 'r') as fi:
		rows = (line.split('\t') for line in fi)
		taxid = {row[0]:row[1] for row in rows}

	# Load fasta records
	with open(fasta, 'r') as handle:
		records = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

	# Create lineage
	lineages = defaultdict()
	for acc, tax in taxid.items():
		lineages[acc] = get_description(tax, taxfoo)
		
	# Write sintax-formatted file
	with open(sintax, 'w') as ouf:
		for acc, lin in lineages.items():
			ouf.write(">" + acc +lin +"\n")
			ouf.write(str(records[acc].seq) + "\n")