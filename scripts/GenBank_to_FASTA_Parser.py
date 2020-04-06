#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import os, sys, argparse

def gbff_2_faa(input_file, output_file_name, tax_name, targets):
	for rec in SeqIO.parse(input_file, "gb"):
		if any(tax_name in s for s in rec.annotations["taxonomy"]):
			for feature in rec.features: #looks for features in the genbank
				for key, val in feature.qualifiers.items(): #looks for val in the feature qualifiers
					for t in targets:
						if any(t in s for s in val):
							seq = rec.seq[feature.location.start:feature.location.end] # retrieve Sequence
							with open(output_file_name, "a") as ofile: #opens the output file and "a" designates it for appending
								ofile.write(">{0} {1}, {2}\n{3}\n".format(rec.name, rec.annotations["organism"], t, seq)) #Writes my FASTA format sequences to the output file
								
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Parse a gbff file to extract sequence entries corresponding to the given taxa and features and outputs a fasta file")
	parser.add_argument("-i", dest = "input_file", type = str, required = True, help ="<required> Input GBFF file path")
	parser.add_argument("-o", dest = "output_file_name", type = str, required = True, help="<required> Output FASTA file path")
	parser.add_argument("-t", dest= "tax_name", type = str, required = True, help = "<required> Name of the taxa to select, use a high rank name to filter several taxa (e.g 'Vertebrata', 'Metazoa')")
	parser.add_argument("-f", dest = "targets", type = str, action = 'append', required = True, help = "<required> Features to extract. e.g. -f 16S -f l-rRNA")
	
	args = parser.parse_args()

	gbff_2_faa(args.input_file, args.output_file_name, args.tax_name, args.targets)