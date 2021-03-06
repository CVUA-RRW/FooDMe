#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from collections import defaultdict
from taxidTools.taxidTools import Taxdump
	
def parse_blast(blast):
	"""
	Parse a BLAST report and returns a dictionnary where Keys are query sequence names and values list of taxids for each hit.
	BLAST report must have the following formatting:
		'6 qseqid sseqid evalue pident bitscore sacc staxids sscinames scomnames stitle'
	"""
	dict = defaultdict()
	with open(blast, 'r') as fi:
		header = fi.readline()
		for line in fi:
			l = line.split()
			taxids = l[6].split(";") # split taxids if nescessary
			# extend taxids list for this OTU
			if l[0] in dict.keys():
				dict[l[0]].extend(taxids)
			# or inititate the list
			else:
				dict[l[0]] = taxids
				
	# Make sure everything is str formated and remove duplicates
	dict = {k: [str(e) for e in set(v)] for k,v in dict.items()}

	return dict

def main(blast_report, output, rankedlineage_dmp, nodes_dmp):
	txd = Taxdump(rankedlineage_dmp, nodes_dmp)
	otu_dict = parse_blast(blast_report)
	with open(output, 'w') as out:
		out.write("queryID\tConsensus\tRank\tTaxid\tDisambiguation\n")
		for queryID, taxid_list in otu_dict.items():	
			# LCA determination
			try:
				lca = txd.lowestCommonNode(taxid_list)
				rank = txd.getRank(lca)
				name = txd.getName(lca)
			
			except KeyError:
				# Taxid not present in the Taxdump version used raises a KeyError
				# Filter out missing sequences (verbose)
				taxid_list_new =[]
				for taxid in taxid_list:
					if taxid not in txd.keys():
						print("WARNING: taxid %s missing from Taxonomy reference, it will be ignored" % taxid)
					else:
						taxid_list_new.append(taxid)
				
				# Update list
				taxid_list = taxid_list_new
							
				# Empty list case:
				if not taxid_list:
					lca = "Undetermined"
					rank = "Undetermined"
					name = "Undetermined"
				else:
					# Get the LCA with the filtered taxids
					lca = txd.lowestCommonNode(taxid_list)
					rank = txd.getRank(lca)
					name = txd.getName(lca)
			
			finally:
				# Format disambiguation list
				name_list = txd.getName(list(set(taxid_list)))
				names = "; ".join(name_list)
				out.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(queryID, name, rank, lca, names))
			
if __name__ == '__main__':
	main(snakemake.input[0], snakemake.output[0], snakemake.params["lineage"], snakemake.params["nodes"])