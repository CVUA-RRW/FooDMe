#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import os, sys

# Extracts all sequences matching 16S and Metazoa from a gbff file and exports a fasta file
# It will create one fasta record for each occurence of 16S feature. Because of messy annotation of gbff files this can lead to duplicate records


input_file = sys.argv[1] #Your GenBank file location. e.g C:\\Sequences\\my_genbank.gb
output_file_name = sys.argv[2] #The name out your fasta output
#accession_numbers = [line.strip() for line in open(input_file)] #the same as your input file, defines the headers for each sequence

if not os.path.exists(output_file_name): #checks for a pre-existing file with the same name as the output
	for rec in SeqIO.parse(input_file, "gb"): #calls the record for the genbank file and SeqIO (BioPython module) to parse it
#		acc = rec.annotations['accessions'][0] #Defines your accession numbers
#		organism = rec.annotations['organism'] #defines your organism ID
#		tax_line = ("| ").join(rec.annotations['taxonomy']) #defines your taxonomy and seperates entries with a |, remove this line, the 'tax_line', and the {2} in your save for a simpler output
		try:
			if rec.annotations["taxonomy"][1] == "Metazoa": #Taxonomy filter
				for feature in rec.features: #looks for features in the genbank
					for key, val in feature.qualifiers.items(): #looks for val in the feature qualifiers
						if any("16S" in s for s in val): #Finds all the features with 16S in them
							seq = rec.seq[feature.location.start:feature.location.end] # retrieve Sequence
							description = rec.description.split(",")[0]+", 16S rRNA"  #Sequence descriptor
							with open(output_file_name, "a") as ofile: #opens the output file and "a" designates it for appending
								ofile.write(">{0} {1}\n{2}\n".format(rec.id, description, seq)) #Writes my FASTA format sequences to the output file
		except:
			continue

else:
	print ("The output file already seem to exist in the current working directory {0}. Please change the name of the output file".format(os.getcwd())) #error code, so you don't overwirite your files