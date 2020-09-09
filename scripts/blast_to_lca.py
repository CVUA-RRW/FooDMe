#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from collections import defaultdict
from taxidTools import Taxdump

want_taxonomy = ['subspecies', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']	
# taxon_filter = '7711' 
	
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
		# Set taxon_filter to one of want_rank
		for queryID, taxid_list in otu_dict.items():
			# Taxon filtering
			# filtered = []
			# for taxid in taxid_list:
				# if taxon_filter in txd.getLineageAsList(taxid):
					# filtered.append(taxid)
			# taxid_list = filtered
			# LCA determination
			try:
				lca = txd.lowestCommonNode(taxid_list)
				rank = txd.getRank(lca)
				name = txd.getName(lca)
				out.write("{0}\t{1}\t{2}\t{3}\n".format(queryID, name, rank, lca))
			except KeyError:
				out.write("{0}\t{1}\t{2}\t{3}\n".format(queryID, "NA", "NA", "NA"))
			
if __name__ == '__main__':
	main(snakemake.input[0], snakemake.output[0], snakemake.params["lineage"], snakemake.params["nodes"])