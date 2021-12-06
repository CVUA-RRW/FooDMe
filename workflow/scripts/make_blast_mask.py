#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import taxidTools as txd


def main(taxid_file, parent, output, taxonomy):

    tax = txd.load(taxonomy)

    with open(taxid_file, "r") as fin:
        db_entries = set(fin.read().splitlines()[1:])

    with open(output, "w") as fout:
        for taxid in db_entries:
            try:
                if tax.isDescendantOf(str(taxid).strip(), str(parent).strip()):
                    fout.write(taxid + "\n")
                else:
                    pass
            except KeyError:
                print(f"WARNING: taxid {taxid} missing from Taxonomy reference"
                      f", it will be ignored")


if __name__ == '__main__':
    main(snakemake.input['taxlist'],
         snakemake.params["taxid"],
         snakemake.output['mask'],
         snakemake.input['tax'])
