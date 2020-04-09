#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, os
from collections import defaultdict
import ncbi_taxdump_utils

want_taxonomy = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']	
	
def parse_blast(blast):
	"""
	Parse a BLAST report and returns a dictionnary where Keys are query sequence names and values list of taxids for each hit.
	BLAST report must have the following formatting:
		'6 qseqid sseqid evalue pident bitscore sacc staxids sscinames scomnames stitle'
	"""
	dict = defaultdict()
	with open(blast, 'r') as fi:
		for line in fi:
			l = line.split()
			if l[0] in dict.keys():
				dict[l[0]].append(l[6])
			else:
				dict[l[0]] = [l[6]]
	return dict

def get_lineage(taxid, nodes_dmp, names_dmp):
	"""
	Returns a lineage dictionnary.
	Keys are taxonomic levels and values are taxonomic names.
	"""
	taxfoo = ncbi_taxdump_utils.NCBI_TaxonomyFoo()    
	taxfoo.load_nodes_dmp(nodes_dmp)
	taxfoo.load_names_dmp(names_dmp)
	
	return taxfoo.get_lineage_as_dict(taxid, want_taxonomy)

def get_consensus(entry):
	"""
	Takes a list of lineage dictionnary and returns the taxnomic consensus level and name.
	"""
	for level in want_taxonomy[::-1]:
		cons = False
		agree = True
		for e in entry:
			if not cons:
				cons = e[level]
			elif e[level] != cons:
				agree = False
		if agree == True:	
			return level, cons

def main(blast_report, output, names_dmp, nodes_dmp):
	otu_dict = parse_blast(blast_report)
	with open(output, 'w') as out:
		out.write("queryID\tConsensus\tRank\n")
	for queryID, taxid_list in otu_dict.items():
		lineages = []
		for taxid in taxid_list:
			lineages.append( get_lineage(taxid, nodes_dmp, names_dmp) )
		level, name = get_consensus(lineages)
		with open(output, 'a') as out:
			out.write("{0}\t{1}\t{2}\n".format(queryID, name, level))
	

main(snakemake.input[0], snakemake.output[0], snakemake.params["names"], snakemake.params["nodes"])