#!/usr/bin/env python3
"""
Python wrapper for the FoodMe Snakemake pipeline.
Creates config file from command-line arguments and starts the pipeline.
"""

import argparse
import os, sys
import datetime
import subprocess


def git_version():
	try:
		version = subprocess.check_output(["git", "describe", "--always"], cwd= os.path.join(sys.path[0],"")).strip().decode("utf-8")
	except subprocess.CalledProcessError:
		version = "Unknown"
	finally:
		return(version)

def create_config(config_file, args):
	if os.path.exists(config_file):
		print("\nWARNING: The file "+config_file+" already exists.")
		print("Overwritting config file\n")
		os.remove(config_file)
	else:
		print("\nCreating config file.\n")
	
	indent1 = " "*4

	with open(config_file, 'w') as conf:
		# File metadata
		conf.write("# This config file was automatically generated by foodme.py\n")
		conf.write("# Version : {}\n".format(git_version()))
		conf.write("# Date: {}\n".format(datetime.datetime.now()))
				
		# Workflow parameters
		conf.write("workdir: {}\n".format(args.working_directory))
		conf.write("samples: {}\n".format(args.sample_list))
		conf.write("threads_sample: {}\n".format(args.threads_sample))
		conf.write("threads: {}\n".format(args.threads))
		
		# Fastp
		conf.write("fastp:\n")
		conf.write("{}length_required: {}\n".format(indent1, args.fastp_length))
		conf.write("{}qualified_quality_phred: {}\n".format(indent1, args.fastp_min_phred))
		conf.write("{}window_size: {}\n".format(indent1, args.fastp_window))
		conf.write("{}mean_quality: {}\n".format(indent1, args.fastp_meanq))
		conf.write("{}primer_fwd: {}\n".format(indent1, args.fastp_prune1))
		conf.write("{}primer_rev: {}\n".format(indent1, args.fastp_prune2))
		conf.write("{}skip_adapter_trimming: {}\n".format(indent1, args.skip_adapter_trimming))
		
		# Read filter
		conf.write("read_filter:\n")
		conf.write("{}min_length: {}\n".format(indent1, args.merge_minlength))
		conf.write("{}max_length: {}\n".format(indent1, args.merge_maxlength))
		conf.write("{}max_expected_errors: {}\n".format(indent1, args.merge_maxee))
		conf.write("{}max_ns: {}\n".format(indent1, args.merge_maxns))
		
		# Cluster
		conf.write("cluster:\n")
		if args.denoise:
			conf.write("{}method: {}\n".format(indent1, "asv"))
		else:
			conf.write("{}method: {}\n".format(indent1, "otu"))
			conf.write("{}cluster_identity: {}\n".format(indent1, args.cluster_id))
			conf.write("{}cluster_minsize: {}\n".format(indent1, args.cluster_minsize))
		
		# Chimera
		conf.write("chimera: {}\n".format(not args.skip_chimera))
		
		# Taxonomy
		conf.write("taxonomy:\n")
		conf.write("{}rankedlineage_dmp: {}\n".format(indent1, args.rankedlineage_dmp))
		conf.write("{}nodes_dmp: {}\n".format(indent1, args.nodes_dmp))

		# Blast
		conf.write("blast:\n")
		conf.write("{}blast_DB: {}\n".format(indent1, args.blastdb))
		conf.write("{}taxdb: {}\n".format(indent1, args.taxdb))
		conf.write("{}taxid_filter: {}\n".format(indent1, args.taxid_filter))
		conf.write("{}e_value: {}\n".format(indent1, args.blast_eval))
		conf.write("{}perc_identity: {}\n".format(indent1, args.blast_id))
		conf.write("{}qcov: {}\n".format(indent1, args.blast_cov))
		conf.write("{}bit_score_diff: {}\n".format(indent1, args.bitscore))
	
def run_snakemake(config_file, args):
	# go to working directory (for proper location of log file)
	os.chdir(args.working_directory)
	
	forceall = ("--forceall" if args.forceall else "")
	dryrun = ("-n" if args.dryrun else "")
	conda_prefix= ("--conda-prefix {}".format(args.condaprefix) if args.condaprefix else "")
	notemp = ("--notemp" if args.keep_temp else "")
	call = "snakemake -s {snakefile} --configfile {config_file} --use-conda --keep-going --cores {cores} {conda_prefix} {notemp} {forceall} {dryrun}".format(snakefile= args.snakefile,
																																				config_file= config_file,
																																				conda_prefix= conda_prefix,
																																				forceall= forceall,
																																				dryrun = dryrun,
																																				cores=args.threads,
																																				notemp=notemp)
	print(call)
	subprocess.call(call, shell=True)
	
def main(): 
	# Action classes to check input parameters
	class FractionType(argparse.Action):
		def __call__(self, parser, namespace, values, option_string=None):
			if values < 0 or values > 1:
				parser.error("Invalid value: '" + self.dest + "' must be between 0 and 1.")
			setattr(namespace, self.dest, values)
	
	class PercentType(argparse.Action):
		def __call__(self, parser, namespace, values, option_string=None):
			if values < 0 or values > 100:
				parser.error("Invalid value: '" + self.dest + "' must be between 0 and 100.")
			setattr(namespace, self.dest, values)
	
	class DatabaseType(argparse.Action):
		def __call__(self, parser, namespace, values, option_string=None):
			if all( [values, not os.path.exists(values), not os.path.exists(values+".nto")] ):
				parser.error("'" + self.dest + "' not found: '" + values + "' does not exist.")
			setattr(namespace, self.dest, values)
	
	# Parser
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, prog = "FooDMe", description= "Another pipeline for (Food) DNA metabarcoding")
	parser.add_argument('-v', '--version', action='version', version="FooDMe version: "+ git_version(), help="Print pipeline version and exit")
	
	# Path arguments
	ioargs = parser.add_argument_group('I/O path arguments')
	ioargs.add_argument('-l', '--sample_list', required=True, type=os.path.abspath, 
						help="Tab-delimited list of samples and paths to read files. Must contain one line of header, each further line contains sample_name, read1_path, read2_path")
	ioargs.add_argument('-d', '--working_directory', required=True, type=os.path.abspath,
						help="Directory to create output files")
						
	# Snakemake arguments
	smkargs = parser.add_argument_group('Snakemake arguments')
	smkargs.add_argument('--forceall', required=False, default=False, action='store_true',
						help="Force the recalculation of all files")
	smkargs.add_argument('-n', '--dryrun', required=False, default=False, action='store_true',
						help="Dryrun. Create config file and calculate the DAG but do not execute anything")
	smkargs.add_argument('-T', '--threads', required=False, default=8, type=int,
						help="Maximum number of threads to use")
	smkargs.add_argument('-t', '--threads_sample', required=False, default=1, type=int,
						help="Number of threads to use per concurent job")
	smkargs.add_argument('-c', '--condaprefix', required=False, type=os.path.abspath, default=False,
						help="Location of stored conda environment. Allows snakemake to reuse environments.")
	smkargs.add_argument('-s', '--snakefile', required=False, type=os.path.abspath, default=os.path.join(sys.path[0], "Snakefile"),
						help="Path to the Snkefile in the FOodMe repo")
	smkargs.add_argument('--keep_temp', required=False, default= False, action='store_true',
						help="Keep large fasta and fastq files, mostly for debug purposes")

	# Fastp
	fastpargs = parser.add_argument_group('Fastp options')
	fastpargs.add_argument('--fastp_length', required=False, default=50, type=int,
						help="Minimum length of input reads to keep")
	fastpargs.add_argument('--fastp_min_phred', required=False, default=15, type=int,
						help="Minimal quality value per base")
	fastpargs.add_argument('--fastp_window', required=False, default=4, type=int,
						help="Size of the sliding window for tail quality trimming")
	fastpargs.add_argument('--fastp_meanq', required=False, default=20, type=int,
						help="Minimum mean Phred-score in the sliding window for tail quality trimming")
	fastpargs.add_argument('--fastp_prune1', required=False, default=0, type=int,
						help="Length of forward primer to prune from 5' end of forward reads (R1)")
	fastpargs.add_argument('--fastp_prune2', required=False, default=0, type=int,
						help="Length of reverse primer to prune from 5' end of reverse reads (R2)")
	fastpargs.add_argument('--skip_adapter_trimming', required=False, default=False, action='store_true',
						help="Skip adapter trimming. Primers provided as an extra fasta files will still be trimmed.")
						
	# Read filter
	readargs = parser.add_argument_group('Merged reads filtering options')
	readargs.add_argument('--merge_minlength', required=False, default=100, type=int,
						help="Minimum length merged reads to keep")
	readargs.add_argument('--merge_maxlength', required=False, default=125, type=int,
						help="Maximum length merged reads to keep")
	readargs.add_argument('--merge_maxee', required=False, default=2, type=int,
						help="Maximum expected errors in merged reads to keep")
	readargs.add_argument('--merge_maxns', required=False, default = 0, type=int,
						help="Maximum number of 'N' base in merged reads. If using denoising procedure this will be automatically reset to 0")
	
	# Cluster
	clsargs = parser.add_argument_group('Clustering options')
	clsargs.add_argument('--denoise', required = False, default=False, action='store_true',
						help="Use denoising instead of identity clustering")
	clsargs.add_argument('--cluster_id', required=False, default=0.97, type=float, action= FractionType,
						help="Minimum identity for clustering sequences in OTUs (between 0 and 1). Will be ignored if using denoising")
	clsargs.add_argument('--cluster_minsize', required=False, default=2, type=int,
						help="Minimal size cutoff for OTUs. Will be ignored if using denoising")
	clsargs.add_argument('--skip_chimera', required=False, default=False, action='store_true',
						help="Skip de novo chimera detection and filtering step")

	#Taxonomy 
	taxo = parser.add_argument_group('Taxonomic assignement files')
	taxo.add_argument('--taxdump', required=False, default=None, type=os.path.abspath, action=DatabaseType,
						help="Path to the taxump folder containing nodes.dmp and rankedlineages.dmp")
	taxo.add_argument('--nodes_dmp', required=False, default=None, type=os.path.abspath, action = DatabaseType,
						help="Path to the nodes.dmp file, needed if --taxdump is omitted")
	taxo.add_argument('--rankedlineage_dmp', required=False, default=None, type=os.path.abspath, action = DatabaseType,
						help="Path to the names.dmp file, needed if --taxdump is omitted")

	# Blast
	blastargs = parser.add_argument_group('Options for BLAST search')
	blastargs.add_argument('--blastdb', required=True, type=os.path.abspath, action = DatabaseType,
						help="Path to the BLAST database, including database basename but no extension (e.g. '/path/to/db/nt')")
	blastargs.add_argument('--taxdb', required=True, type=os.path.abspath, action = DatabaseType,
						help="Path to the BLAST taxonomy database (folder)")
	blastargs.add_argument('--taxid_filter', required=False, type=str, default=None,
						help="Limit BLAST search to the taxids under the given node")
	blastargs.add_argument('--blast_eval', required=False, default=1e-10, type=float,
						help="E-value threshold for blast results")
	blastargs.add_argument('--blast_id', required=False, default=90, type=float, action= PercentType,
						help="Minimal identity between the hit and query for blast results (in percent)")
	blastargs.add_argument('--blast_cov', required=False, default=90, type=float, action= PercentType,
						help="Minimal proportion of the query covered by a hit for blast results. A mismatch is still counting as covering (in percent)")
	blastargs.add_argument('--bitscore', required=False, default=4, type=int,
						help="Maximum bit-score difference with the best hit for a blast result to be included in the taxonomy consensus detemination")
	
	args = parser.parse_args()
	
	# Check Taxdump files
	if not args.taxdump:
		if (not args.nodes_dmp) or (not args.rankedlineage_dmp):
			raise parser.error(message="Argument error: Specify either --taxdump or --nodes_dmp and --rankedlineage_dmp")
	else:
		if not (os.path.exists(os.path.join(args.taxdump, "nodes.dmp")) and os.path.exists(os.path.join(args.taxdump, "rankedlineage.dmp"))):
			raise parser.error(message="Could not find nodes.dmp or rankedlineage.dmp in "+args.taxdump+"\nSpecify the paths using --nodes_dmp and --rankedlineage_dmp")
		args.__dict__["nodes_dmp"] = os.path.join(args.taxdump, "nodes.dmp")
		args.__dict__["rankedlineage_dmp"] = os.path.join(args.taxdump, "rankedlineage.dmp")

	# Create workdir
	if not os.path.exists(args.working_directory):
		os.makedirs(args.working_directory)
	
	# Create config.yaml
	config_file = os.path.join(args.working_directory, "config.yaml")
	create_config(config_file, args)

	# Execute snakemake
	run_snakemake(config_file, args)
	
	# On quit
	print("\nThank you for using FooDMe!\n")

if __name__=='__main__':
		main()