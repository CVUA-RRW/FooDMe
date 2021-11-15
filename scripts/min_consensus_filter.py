#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import taxidTools as txd
from collections import Counter, defaultdict


def parse_blast(blast_file):
    """
    Parse a BLAST report and returns a dictionnary where Keys are query sequence names and values list of taxids for each hit.
    BLAST report must have the following formatting:
        '6 qseqid sseqid evalue pident bitscore sacc staxids sscinames scomnames stitle'
    """
    dictout = defaultdict()
    with open(blast_file, 'r') as fi:
        header = fi.readline()
        for line in fi:
            l = line.split()
            taxids = l[6].split(";") # split taxids if nescessary
            # extend taxids list for this OTU
            if l[0] in dictout.keys():
                dictout[l[0]].extend(taxids)
            # or inititate the list
            else:
                dictout[l[0]] = taxids
                
    # Make sure everything is str formated
    dictout = {k: [str(e) for e in v] for k,v in dictout.items()}
    
    return dictout

def main(blast_report, output, min_consensus, rankedlineage_dmp, nodes_dmp):
    if min_consensus <= 0.5 or min_consensus >1:
        raise ValueError("'min_consensus' must be in the interval (0.5 , 1]")
    
    tax = txd.Taxonomy.from_taxdump(nodes_dmp, rankedlineage_dmp)
    otu_dict = parse_blast(blast_report)
    with open(output, 'w') as out:
        out.write("queryID\tConsensus\tRank\tTaxid\tDisambiguation\n")
        
        for queryID, taxid_list in otu_dict.items():
            try:
                consensus = tax.consensus(taxid_list, min_consensus)
            
            except KeyError:
                # Taxid not present in the Taxdump version used raises a KeyError
                # Filter out missing sequences (verbose)
                taxid_list_new =[]
                for taxid in taxid_list:
                    if taxid not in tax.keys():
                        print("WARNING: taxid %s missing from Taxonomy reference, it will be ignored" % taxid)
                    else:
                        taxid_list_new.append(taxid)
                
                # Update list
                taxid_list = taxid_list_new
                            
                # Empty list case:
                if not taxid_list:
                    consensus = "Undetermined"
                else:
                    # Get the consensus with the filtered taxids
                    consensus = tax.consensus(taxid_list, min_consensus)
            
            finally:
                # Get consensus Info
                try:
                    rank = consensus.rank
                    name = consensus.name
                    taxid = consensus.taxid
                except KeyError:
                    taxid = "Undetermined"
                    rank = "Undetermined"
                    name = "Undetermined"
                
                # Format disambiguation list
                name_list = [tax.getName(node) for node in list(set(taxid_list))]
                freqs = [((v/len(taxid_list)),tax.getName(k)) for k,v in Counter(taxid_list).items()] # (freq, name) tuple to sort
                sorted_freqs = sorted(freqs, reverse = True)

                names = "; ".join([f"{f} ({round(n,2)})" for (n,f) in sorted_freqs])
                out.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(queryID, name, rank, taxid, names))

if __name__ == '__main__':
    main(snakemake.input[0], snakemake.output[0], snakemake.params["min_consensus"], snakemake.params["lineage"], snakemake.params["nodes"])
