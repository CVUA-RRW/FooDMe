#!/usr/bin/env bash
#set -e
#set -u
#set -o pipefail

#==================================================================================
# BSD 3-Clause License
#
# Copyright (c) 2019, Carlus Denecke and Simon H. Tausch
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Source: gitlab.com/bfr_bioinformatics/AQUAMIS
#=================================================================================

#Goal: Create sample sheet for all fastq files in a folder
#Author: Carlus Deneke, Carlus.Deneke@bfr.bund.de
version=0.4.1

# Read in parameters --------------------------------

getopt --test > /dev/null
if [[ $? -ne 4 ]]; then
    echo "I’m sorry, `getopt --test` failed in this environment."
    exit 1
fi


OPTIONS=hdm:f:o:iF
LONGOPTIONS=help,dryrun,nocheck,mode:,fastxDir:,outDir:,interactive,force
#,keep-undetermined

# DEFAULT
check=true
#keep-undetermined=false


## Helpfile and escape
if [ $# -eq 0 ]; then
    echo "Please provide at least one argument. For the help file, type create_sampleSheet.sh --help"
fi

# -----------------------------

PARSED=$(getopt --options=${OPTIONS} --longoptions=${LONGOPTIONS} --name "$0" -- "$@")
if [[ $? -ne 0 ]]; then
    # e.g. $? == 1
    #  then getopt has complained about wrong arguments to stdout
    exit 2
fi
# read getopt’s output this way to handle the quoting right:
eval set -- "${PARSED}"


while true; do
    case "$1" in
        -h|--help)
            help=true
            shift
            ;;
        -d|--dryrun)
            dryrun=true
            shift
            ;;
        --nocheck)
            check=false
            shift
            ;;
        #--keep-undetermined)
            #keep-undetermined=true
            #shift
            #;;
        -m|--mode)
            mode=${2,,} #convert to lower case
            shift 2
            ;;
        -f|--fastxDir)
            fastxDir="$2"
            shift 2
            ;;
        -o|--outDir)
            outDir="$2"
            shift 2
            ;;
        -i|--interactive)
            interactive=true
            shift
            ;;
        -F|--force)
            force=true
            shift
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Programming error"
            exit 3
            ;;
    esac
done


# Checks ---------------------------

# call help
if [[ ${help} == true ]]
then
    echo "You called the script create_sampleSheet.sh (version $version). Purpose: Create sample sheet for all fastq files in a specified folder"
    echo
    echo "============================="
    echo "Call: create_sampleSheet.sh --mode {illumina, trimmed, ncbi, flex, assembly} --fastxDir path/fastq/dir --outDir path/out/dir [Options]"
    echo "--mode: Choose mode from illumina, ncbi, flex, assembly  (default: illumina)"
    echo "--fastxDir: Path to existing directory containing the fastq or fasta files (default: `pwd`)"
    echo "--outDir: Path to existing outDir (default: fastxDir)"
    echo 
    echo "Options:"
    echo "--nocheck: Do not check consistency of sample sheet"
    #echo "--keep-undetermined: Also keep undetermined fastq files"
    echo "--interactive: Ask before starting the program"
    echo "--force: Overwrite existing samples.tsv files in OutDir"
    echo "--help: Display this help message"
    echo "--dryrun: Perform a dry-run"
    echo 
    echo "Details about --mode parameters"
    echo "illumina: samples are in illumina format: {samplename}_S*_R{1,2}_001.fastq*"
    echo "ncbi: samples are in ncbi format: {samplename}_{1,2}.fastq.gz"
    echo "flex: samples are in the following format: {samplename}*_R{1,2}*.fastq*. The sample name is cut after the first \"_\". If your sample name contains \"_\" the sample name will be cropped!"
    echo "assembly: samples are in format {samplename}.fasta Note that currently only the extension \".fasta\" is supported"
    echo


    echo "============================="
    exit 1

fi



# check if empty
if [[ ${mode} == "" ]]
then
  mode="illumina"
fi

# check if empty
if [[ ${fastxDir} == "" ]]
then
  fastxDir=`pwd`
fi

# check if emtpy
if [[ ${outDir} == "" ]]
then
  outDir=${fastxDir}
fi

# check if outdir exists
if [[ ! -d ${outDir} ]]
then
  echo "The output Dir $outDir does NOT exist"; 
  [[ "$0" = "$BASH_SOURCE" ]] && exit 1 || return 1
fi

# check if fastxDir exists
if [[ ! -d ${fastxDir} ]]
then
  echo "The fastq Dir $fastxDir does NOT exist"
  [[ "$0" = "$BASH_SOURCE" ]] && exit 1 || return 1
fi

outfile="$outDir/samples.tsv"

# check if outfile already present
if [[ -f ${outfile} && ${force} != true ]]
then 
   echo "The file $outfile already exists. Use --force option to force sample sheet creation"
   [[ "$0" = "$BASH_SOURCE" ]] && exit 1 || return 1
fi

# check dry-run
if [[ ${dryrun} == true ]]
then
   echo "Performing a dry-run"
   outfile="/dev/null"
fi

# check mode
if [[ ${mode} =~ "illu" ]]
then
  mode="illumina"
elif [[ ${mode} =~ "fle" ]]
then
  mode="flex"
elif [[ ${mode} =~ "name" ]]
then
  mode="names"
elif [[ ${mode} =~ "nc" ]]
then
  mode="ncbi"
elif [[ ${mode} =~ "ass" ]]
then
  mode="assembly"
else
  echo "No correct mode selected. Choose between illumina, trimmed, ncbi, flex and name";
  [[ "$0" = "$BASH_SOURCE" ]] && exit 1 || return 1
fi

echo "You chose mode $mode"




echo "interactive: $interactive, check: $check, force: $force, mode: $mode, out: $outDir, fastxDir: $fastxDir"


# ask if in interactive mode
if [[ ${interactive} == true ]]
then
	read -p "Do you want to continue? " -n 1 -r
	echo
	if [[ ! $REPLY =~ ^[Yy]$ ]]
	then
		[[ "$0" = "$BASH_SOURCE" ]] && exit 1 || return 1
	fi
fi



# run -----------------------------


if [ ${mode} = "names" ]
then
	# create simple sample list
	ls --color=never *.fastq.gz | cut -d '_' -f 1 | sort | uniq > ${outfile}
elif [ ${mode} = "illumina" ]
then
	# using find for illumina names data
    echo -e "sample\tfq1\tfq2" > ${outfile}
    for file in ${fastxDir}/*_S*_R1_001.fastq*; do
        sample=`basename ${file} | awk -F '_S[0-9]' '{print $1}'` # more flexible with wildcard in separator
        R2=`echo ${file} | sed 's/_R1_001.fastq/_R2_001.fastq/'`
        if [[ ! -f ${R2} ]]; then
            echo "Reverse read $R2 does NOT exist"
            [[ "$0" = "$BASH_SOURCE" ]] && exit 1 || return 1
        fi
        echo -e "$sample\t$file\t$R2" >> ${outfile}
    done
    
elif [ ${mode} = "flex" ]; then
    echo -e "sample\tfq1\tfq2" > ${outfile}
    for file in ${fastxDir}/*_R1*.fastq*; do
        sample=`basename ${file} | cut -f 1 -d '_'`
        R2=`echo ${file} | sed 's/_R1/_R2/'`
        if [[ ! -f ${R2} ]]; then
            echo "Reverse read $R2 does NOT exist"
            [[ "$0" = "$BASH_SOURCE" ]] && exit 1 || return 1
        fi
        echo -e "$sample\t$file\t$R2" >> ${outfile}
    done
elif [ ${mode} = "ncbi" ]
then
	# using find for ncbi names "_1", "_2"
    echo -e "sample\tfq1\tfq2" > ${outfile}
    for file in ${fastxDir}/*_1.fastq*; do
        sample=`basename ${file} | cut -f 1 -d '_'`
        R2=`echo ${file} | sed 's/_1.fastq/_2.fastq/'`
        if [[ ! -f ${R2} ]]; then
            echo "Reverse read $R2 does NOT exist"
            [[ "$0" = "$BASH_SOURCE" ]] && exit 1 || return 1
        fi
        echo -e "$sample\t$file\t$R2" >> ${outfile}
    done
elif [ ${mode} = "assembly" ]
then
	echo -e "sample\tassembly" > ${outfile}
	for file in ${fastxDir}/*.fasta; do sample=`basename -s .fasta ${file}`; echo -e "$sample\t$file" >> ${outfile}; done
	#for file in $fastxDir/*.fna; do sample=`basename -s .fna $file`; echo -e "$sample\t$file" >> $outfile; done
    #for file in $fastxDir/*.{fasta,fna}; do sample=$(basename -s .fasta $(basename -s .fna $file) ); echo -e "$sample\t$file" >> $outfile; done
else
	echo "mode not specified or recognized.Please specify mode in command line."
	echo "Choose from illumina, ncbi, flex, assembly"
	exit 1

fi



# --------------
# remove undetermined 
# if [[ ! $keep-undetermined = true && $dryrun != true ]]
if [[ ${dryrun} != true ]]
then
	sed -i '/Undetermined_/d' ${outfile}
fi


# check ---------------
if [[ ${check} == true ]]; then
    status="PASS"
    while read sample read1 read2; do
        #if ! [[ ${read1} =~ "${sample}_" ]]; then # only for fastq files
        if ! [[ ${read1} =~ "${sample}_" || ${read1} =~ "${sample}.fasta"  || ${read1} =~ "${sample}.fna" ]]; then
            echo "Error: $sample not contained in $read1"
            status="FAIL"
            #[[ "$0" = "$BASH_SOURCE" ]] && exit 1 || return 1
        fi
    done < <(tail -n +2 ${outfile})
    
    #check if empty or same sample name more than once
    dups=`tail -n +2 ${outfile} | cut -f 1 | sort | uniq -c | sed -E 's/^[ ]+//' | cut -f 1 -d ' ' | awk '$1 != 1' | wc -l`
    if [[ $dups > 0 ]]; then 
        echo "Error: Duplicated sample names. Check if your samples are formated according to the selected mode"; 
        status="FAIL"
    fi

    if [[ $status == "PASS" ]]; then
        echo "Check successful"
    else
        echo "Check NOT successful"
        rm ${outfile}
    fi


fi


# ---------------


if [[ ${dryrun} != true && $status != "FAIL" ]]
then
	samplesfound=`tail -n +2 ${outfile} | wc -l `
	echo "Found $samplesfound samples"
fi



