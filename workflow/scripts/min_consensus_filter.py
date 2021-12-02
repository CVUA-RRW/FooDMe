#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import taxidTools as txd
from collections import Counter, defaultdict


def parse_blast(blast_file):
    """
    Parse a BLAST report and returns a dictionnary where Keys are query
    sequence names and values list of taxids for each hit.
    BLAST report must have the following formatting:
        '6 qseqid sseqid evalue pident bitscore sacc
        staxids sscinames scomnames stitle'
    """
    dictout = defaultdict()
    with open(blast_file, 'r') as fi:
        for line in fi:
            ls = line.split()
            taxids = ls[6].split(";")  # split taxids if nescessary
            # extend taxids list for this OTU
            if ls[0] in dictout.keys():
                dictout[ls[0]].extend(taxids)
            # or inititate the list
            else:
                dictout[ls[0]] = taxids

    # Make sure everything is str formated
    dictout = {k: [str(e) for e in v] for k, v in dictout.items()}

    return dictout


def main(blast_report, output, min_consensus, taxonomy):
    if min_consensus <= 0.5 or min_consensus > 1:
        raise ValueError("'min_consensus' must be in the interval (0.5 , 1]")

    tax = txd.load(taxonomy)
    otu_dict = parse_blast(blast_report)
    with open(output, 'w') as out:
        out.write("queryID\tConsensus\tRank\tTaxid\tDisambiguation\n")

        for queryID, taxid_list in otu_dict.items():
            try:
                consensus = tax.consensus(taxid_list, min_consensus)

            except KeyError:
                # Taxid not present in the Taxdump version
                # used raises a KeyError
                # Filter out missing sequences (verbose)
                taxid_list_new = []
                for taxid in taxid_list:
                    if taxid not in tax.keys():
                        print(f"WARNING: taxid {taxid} missing from Taxonomy "
                              f"reference, it will be ignored")
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

                # (freq, name) tuple to sort
                freqs = [((v/len(taxid_list)), tax.getName(k))
                         for k, v in Counter(taxid_list).items()]
                sorted_freqs = sorted(freqs, reverse=True)

                names = "; ".join([f"{f} ({round(n,2)})"
                                   for (n, f) in sorted_freqs])
                out.write(f"{queryID}\t{name}\t{rank}\t{taxid}\t{names}\n")


if __name__ == '__main__':
    main(snakemake.input['blast'],
         snakemake.output['consensus'],
         snakemake.params["min_consensus"],
         snakemake.input['tax'])
