#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from taxidTools.taxidTools import Taxdump
import pandas as pd

def get_lineage(taxid, txd):
	if taxid == "-":
		return ["Unassigned"]
	elif taxid == "Undetermined":
		return ["Undetermined"]
	else:
		return txd.getName(txd.getAncestry(str(taxid)))[::-1] # inverting list to have the lineage descending for Krona

def main(input, output, rankedlineage_dmp, nodes_dmp):
	txd = Taxdump(rankedlineage_dmp, nodes_dmp)
	df = pd.read_csv(input, sep='\t', header=0)
	with open(output, "w") as out:
		for index, row in df.iterrows():
			out.write("\t".join([str(row["Count"])] + get_lineage(row['Taxid'], txd)) + "\n")
	
if __name__ == '__main__':	
	main(snakemake.input[0], snakemake.output[0], snakemake.params["lineage"], snakemake.params["nodes"])