#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# from collections import defaultdict
import ncbi_taxdump_utils
import re

# want_taxonomy = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']	

def init_tax(names_dmp, nodes_dmp):
	taxfoo = ncbi_taxdump_utils.NCBI_TaxonomyFoo()  
	taxfoo.load_nodes_dmp(nodes_dmp)
	taxfoo.load_names_dmp(names_dmp)
	return taxfoo
	
def main(sintax_report, output, names_dmp, nodes_dmp):
	taxfoo = init_tax(names_dmp, nodes_dmp)
	with open(sintax_report, 'r') as sin, open(output, 'w') as out:
		out.write("Query\tCount\tBootstrap support\tConsensus\tRank\tTaxid\n")
		for line in sin:
			record = line.strip().split("\t")
			seq, size = record[0].split(";size=")
			if len(record) > 1:
				taxid = record[-1].split(",")[-1].split(":")[-1]
				cons = taxfoo.get_taxid_name(int(taxid))
				rank = taxfoo.get_taxid_rank(int(taxid))
				bootstrap = re.search(rf"{taxid}\((\d\.\d+)\)", record[1]).group(1)
			else:
				taxid, cons, rank, bootstrap = "-" , "-" , "-" , "0"
				
			out.write("\t".join([seq, size, bootstrap, cons, rank, str(taxid)]) +"\n")
				
if __name__ == '__main__':
	main(snakemake.input[0], snakemake.output[0], snakemake.params["names"], snakemake.params["nodes"])