#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# from collections import defaultdict
from taxidTools import Taxdump
import re

want_taxonomy = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']
	
def main(sintax_report, output, rankedlineage_dmp, nodes_dmp):
	txd = Taxdump(rankedlineage_dmp, nodes_dmp, want_taxonomy)
	with open(sintax_report, 'r') as sin, open(output, 'w') as out:
		out.write("Query\tCount\tBootstrap support\tConsensus\tRank\tTaxid\n")
		for line in sin:
			record = line.strip().split("\t")
			seq, size = record[0].split(";size=")
			if len(record) > 1:
				taxid = record[-1].split(",")[-1].split(":")[-1]
				cons = txd.getName(taxid)
				rank = txd.getRank(taxid)
				bootstrap = re.search(rf"{taxid}\((\d\.\d+)\)", record[1]).group(1)
			else:
				taxid, cons, rank, bootstrap = "-" , "-" , "-" , "0"
				
			out.write("\t".join([seq, size, bootstrap, cons, rank, str(taxid)]) +"\n")
				
if __name__ == '__main__':
	main(snakemake.input[0], snakemake.output[0], snakemake.params["lineage"], snakemake.params["nodes"])