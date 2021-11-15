#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import taxidTools as txd
import pandas as pd


def get_lineage(taxid, tax):
    if taxid == "-":
        return ["Unassigned"]
    elif taxid == "Undetermined":
        return ["Undetermined"]
    else:
        return [node.name for node in tax.getAncestry(taxid)] # inverting list to have the lineage descending for Krona

def main(input, output, rankedlineage_dmp, nodes_dmp):
    tax = txd.Taxonomy.from_taxdump(nodes_dmp, rankedlineage_dmp)
    df = pd.read_csv(input, sep='\t', header=0)
    with open(output, "w") as out:
        for index, row in df.iterrows():
            out.write("\t".join([str(row["Count"])] + get_lineage(row['Taxid'], tax)) + "\n")
    
if __name__ == '__main__':  
    main(snakemake.input[0], snakemake.output[0], snakemake.params["lineage"], snakemake.params["nodes"])