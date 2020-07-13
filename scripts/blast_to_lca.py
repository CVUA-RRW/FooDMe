#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from collections import defaultdict
from taxidTools import Taxdump

want_taxonomy = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']	
	
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

def main(blast_report, output, rankedlineage_dmp, nodes_dmp):
	txd = Taxdump(rankedlineage_dmp, nodes_dmp, want_taxonomy)
	otu_dict = parse_blast(blast_report)
	with open(output, 'w') as out:
		out.write("queryID\tConsensus\tRank\tTaxid\n")
		for queryID, taxid_list in otu_dict.items():
			lca = txd.lowestCommonNode(taxid_list)
			rank = txd.getRank(lca)
			name = txd.getName(lca)
			out.write("{0}\t{1}\t{2}\t{3}\n".format(queryID, name, rank, lca))
			
if __name__ == '__main__':
	main(snakemake.input[0], snakemake.output[0], snakemake.params["lineage"], snakemake.params["nodes"])