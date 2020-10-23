#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from taxidTools import Taxdump

want_taxonomy = ["forma", "varietas", "subspecies", "species"]
	
def main(taxid_file, parent, output, rankedlineage_dmp, nodes_dmp):
	txd = Taxdump(rankedlineage_dmp, nodes_dmp, want_taxonomy)
	
	with open(taxid_file, "r") as fin:
		db_entries = set(fin.read().splitlines()[1:])
	
	with open(output, "w") as fout:
		for taxid in db_entries:
			try:
				if str(parent) in txd.getFullLineage(str(taxid).strip()).values():
					fout.write(taxid + "\n")
				else:
					pass
			except KeyError:
				#print("WARNING: taxid %s missing from Taxonomy reference, it will be ignored" % taxid)
				continue
				
if __name__ == '__main__':
	main(snakemake.input[0], snakemake.params["taxid"], snakemake.output[0], snakemake.params["lineage"], snakemake.params["nodes"])