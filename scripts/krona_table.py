#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import ncbi_taxdump_utils
import pandas as pd

want_taxonomy = ['phylum', 'class', 'order', 'family', 'genus', 'species']	

def init_tax(names_dmp, nodes_dmp):
	taxfoo = ncbi_taxdump_utils.NCBI_TaxonomyFoo()  
	taxfoo.load_nodes_dmp(nodes_dmp)
	taxfoo.load_names_dmp(names_dmp)
	return taxfoo
	
def get_lineage(taxid, taxfoo):
	if taxid == "-":
		return ["Unassigned"]
	else:
		return taxfoo.get_lineage(taxid, want_taxonomy)

def main(input, output, names_dmp, nodes_dmp):
	taxfoo = init_tax(names_dmp, nodes_dmp)
	df = pd.read_csv(input, sep='\t', header=0)
	with open(output, "w") as out:
#		out.write("\t".join(["count"] + want_taxonomy) + "\n")
		for index, row in df.iterrows():
			out.write("\t".join([str(row["Count"])] + get_lineage(row['Taxid'], taxfoo)) + "\n")
	
if __name__ == '__main__':	
	main(snakemake.input[0], snakemake.output[0], snakemake.params["names"], snakemake.params["nodes"])