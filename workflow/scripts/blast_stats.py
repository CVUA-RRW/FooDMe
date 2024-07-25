#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


sys.stderr = open(snakemake.log[0], "w")


import os
import pandas as pd


HEADER = "\t".join([
    "Sample",
    "Query",
    "Count",
    "Blast hits",
    "Best bit-score",
    "Lowest bit-score",
    "Bit-score threshold",
    "Saved Blast hits",
    "Consensus",
    "Rank",
    "Taxid",
    "Disambiguation",
    "link_report",
    "link_filtered"
    ])


def parse_fasta(fasta):
    ids, sizes = [], []
    with open(fasta, "r") as fi:
        for l in fi:
            if l[0] ==">":
                ids.append(l.split(";")[0][1:].strip())
                sizes.append(l.split("=")[1].strip())
    return ids, sizes


def count_blasthits(id, report):
    i = 0
    with open(report, "r") as fi:
        for l in fi:
            if l.split(";")[0] == id:
                i += 1
    return i


def main(otus_in, blast_in, filtered_in, lca_in, report_out, bit_diff, sample):
        with open(report_out, "w") as fo:
                fo.write(f"{HEADER}\n")
        if not os.path.isfile(blast_in) or os.stat(blast_in).st_size == 0:
            with open(report_out, "a") as fo:
                fo.write("\t".join([
                    str(sample),
                    "-",
                    "-",
                    "0",
                    "0",
                    "0",
                    "0",
                    "0",
                    "-",
                    "-",
                    "-",
                    "-",
                    f"../{blast_in}",
                    f"../{filtered_in}"
                ])+"\n")
        
        else:
            otus, sizes = parse_fasta(otus_in)
            blast_df = pd.read_csv(blast_in, sep="\t")
            filtered_df = pd.read_csv(filtered_in, sep="\t")
            lca_df = pd.read_csv(lca_in, sep="\t")

            for otu, size in zip(otus, sizes):
                bhits = count_blasthits(otu, blast_in)
                if bhits == 0:
                    with open(report_out, "a") as fo:
                        fo.write("\t".join([
                            str(sample),
                            str(otu),
                            str(size),
                            "0",
                            "0",
                            "0",
                            "0",
                            "0",
                            "-",
                            "-",
                            "-",
                            "- (1.0)",
                            f"../{blast_in}",
                            f"../{filtered_in}"
                        ])+"\n")
                else:
                    bit_best = blast_df["bitscore"].max()
                    bit_low = blast_df["bitscore"].min()
                    bit_thr = bit_best - bit_diff
                    shits = filtered_df["query"].str.count(otu).sum()
                    cons = lca_df[lca_df["queryID"] == otu]

                    with open(report_out, "a") as fo:
                        fo.write("\t".join([
                            str(sample),
                            str(otu),
                            str(size),
                            str(bhits),
                            str(bit_best),
                            str(bit_low),
                            str(bit_thr),
                            str(shits),
                            str(cons["Consensus"].values[0]),
                            str(cons["Rank"].values[0]),
                            str(cons["Taxid"].values[0]),
                            str(cons["Disambiguation"].values[0]),
                            f"../{blast_in}",
                            f"../{filtered_in}"
                        ])+"\n")


if __name__ == '__main__':
        main(
                snakemake.input['otus'],
                snakemake.input['blast'],
                snakemake.input['filtered'],
                snakemake.input['lca'],
                snakemake.output[0],
                snakemake.params['bit_diff'],
                snakemake.params['sample'],
        )
