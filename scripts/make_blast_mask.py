#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from taxidTools import Taxdump

want_taxonomy = ["forma", "varietas", "subspecies", "species"]

def get_children_list(txd, taxid):
	"""
	Returns a list of unque taxid lower than the node given as an input
	"""
	d = txd.getChildren(str(taxid), want_ranks=want_taxonomy)
	d2 = [e for dict in d for e in dict.values()]
	d3 = set([e for e in d2 if e])
	return d3
	
def main(taxid, output, rankedlineage_dmp, nodes_dmp):
	txd = Taxdump(rankedlineage_dmp, nodes_dmp, want_taxonomy)
	with open(output, "w") as fout:
		for e in get_children_list(txd, taxid):
			fout.write(e + "\n")

if __name__ == '__main__':
	main(snakemake.params["taxid"], snakemake.output[0], snakemake.params["lineage"], snakemake.params["nodes"])