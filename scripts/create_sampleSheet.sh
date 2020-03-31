#!/usr/bin/env bash
#set -e
#set -u
#set -o pipefail

#Goal: Create sample sheet for all fastq files in a folder
#Author: Carlus Deneke, Carlus.Deneke@bfr.bund.de
version=0.4

# Read in parameters --------------------------------

getopt --test > /dev/null
if [[ $? -ne 4 ]]; then
    echo "I’m sorry, `getopt --test` failed in this environment."
    exit 1
fi


OPTIONS=hdcm:f:o:iF
LONGOPTIONS=help,dryrun,check,mode:,fastqDir:,outDir:,interactive,force

#TODO remove header --noheader
#TODO make this optional: keep-undetermined

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
        -c|--check)
            check=true
            shift
            ;;
        -m|--mode)
            mode=${2,,} #convert to lower case
            shift 2
            ;;
        -f|--fastqDir)
            fastqDir="$2"
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
    echo "Call: create_sampleSheet.sh --mode {illumina, trimmed, ncbi, flex, assembly} --fastqDir path/fastq/dir --outDir path/out/dir [Options]"
    echo "--mode: Choose mode from illumina, trimmed, ncbi, flex, assembly  (default: illumina)"
    echo "--fastqDir: Path to existing directory containing the fastq files (default: `pwd`)"
    echo "--outDir: Path to existing outDir (default: fastqDir)"
    echo 
    echo "Options:"
    echo "--check: Check consistency of sample sheet"
    echo "--interactive: Ask before starting the program"
    echo "--force: Overwrite existing samples.tsv files in OutDir"
    echo "--help: Display this help message"
    echo "--dryrun: Perform a dry-run"

    echo "============================="
    exit 1

fi



# check if empty
if [[ ${mode} == "" ]]
then
  mode="illumina"
fi

# check if empty
if [[ ${fastqDir} == "" ]]
then
  fastqDir=`pwd`
fi

# check if emtpy
if [[ ${outDir} == "" ]]
then
  outDir=${fastqDir}
fi

# check if outdir exists
if [[ ! -d ${outDir} ]]
then
  echo "The output Dir $outDir does NOT exist"; 
  [[ "$0" = "$BASH_SOURCE" ]] && exit 1 || return 1
fi

# check if fastqDir exists
if [[ ! -d ${fastqDir} ]]
then
  echo "The fastq Dir $fastqDir does NOT exist"
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
elif [[ ${mode} =~ "trim" ]]
then
  mode="trimmed"
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


# always check results
check=true


echo "interactive: $interactive, check: $check, force: $force, mode: $mode, out: $outDir, fastqDir: $fastqDir"


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
	#echo -e "sample\tfq1\tfq2" > $outfile
	#paste <(find $fastqDir -name "*_R1*.fastq.gz" -exec basename {} \; | cut -d '_' -f 1 | sort | uniq) <(find $fastqDir -name "*_R1*.fastq.gz" | sort) <(find $fastqDir -name "*_R2*.fastq.gz" | sort) --delimiters '\t'  >> $outfile
    echo -e "sample\tfq1\tfq2" > ${outfile}
    for file in ${fastqDir}/*_S*_R1_001.fastq*; do
        sample=`basename ${file} | awk -F '_S' '{print $1}'`
        R2=`echo ${file} | sed 's/_R1_001.fastq/_R2_001.fastq/'`
        if [[ ! -f ${R2} ]]; then
            echo "Reverse read $R2 does NOT exist"
            [[ "$0" = "$BASH_SOURCE" ]] && exit 1 || return 1
        fi
        echo -e "$sample\t$file\t$R2" >> ${outfile}
    done
    
elif [ ${mode} = "flex" ]; then
    echo -e "sample\tfq1\tfq2" > ${outfile}
    for file in ${fastqDir}/*_R1*.fastq*; do
        sample=`basename ${file} | cut -f 1 -d '_'`
        R2=`echo ${file} | sed 's/_R1/_R2/'`
        if [[ ! -f ${R2} ]]; then
            echo "Reverse read $R2 does NOT exist"
            [[ "$0" = "$BASH_SOURCE" ]] && exit 1 || return 1
        fi
        echo -e "$sample\t$file\t$R2" >> ${outfile}
    done
elif [ ${mode} = "trimmed" ]
then
	# for trimmed data
	echo -e "sample\tfq1\tfq2\tfq1U\tfq2U" > ${outfile}
    for file in ${fastqDir}/*1P.fastq*; do
        sample=`basename ${file} | cut -f 1 -d '_'`
        R2=`echo ${file} | sed 's/1P.fastq/2P.fastq/'`
        if [[ ! -f ${R2} ]]; then
            echo "Reverse read $R2 does NOT exist"
            [[ "$0" = "$BASH_SOURCE" ]] && exit 1 || return 1
        fi
        U1=`echo ${file} | sed 's/1P.fastq/1U.fastq/'`
        if [[ ! -f ${U1} ]]; then
            echo "Reverse read $U1 does NOT exist"
            [[ "$0" = "$BASH_SOURCE" ]] && exit 1 || return 1
        fi
        U2=`echo ${file} | sed 's/1P.fastq/2U.fastq/'`
        if [[ ! -f ${U2} ]]; then
            echo "Reverse read $U2 does NOT exist"
            [[ "$0" = "$BASH_SOURCE" ]] && exit 1 || return 1
        fi

        echo -e "$sample\t$file\t$R2\t$U1\t$U2" >> ${outfile}
    done
  #paste <(find $fastqDir -name "*1P.fastq.gz" -exec basename {} \; | cut -d '.' -f 1 | cut -d '_' -f 1 | sort | uniq) <(find $fastqDir -name "*1P.fastq.gz" | sort) <(find $fastqDir -name "*2P.fastq.gz" | sort) <(find $fastqDir -name "*1U.fastq.gz" -o -name "*r1_UP_*.fastq" | sort) <(find $fastqDir -name "*2U.fastq.gz" | sort) --delimiters '\t'  >> $outfile
elif [ ${mode} = "ncbi" ]
then
	# using find for ncbi names "_1", "_2"
#	echo -e "sample\tfq1\tfq2" > $outfile
#	paste <(find $fastqDir -name "*_1.fastq.gz" -exec basename {} \; | cut -d '_' -f 1 | sort | uniq) <(find $fastqDir -name "*_1.fastq.gz" -o -name "*_1.fastq" | sort) <(find $fastqDir -name "*_2.fastq.gz" -o -name "*_2.fastq" | sort) --delimiters '\t'  >> $outfile
    echo -e "sample\tfq1\tfq2" > ${outfile}
    for file in ${fastqDir}/*_1.fastq*; do
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
	for file in ${fastqDir}/*.fasta; do sample=`basename -s .fasta ${file}`; echo -e "$sample\t$file" >> ${outfile}; done
	#for file in $fastqDir/*.fna; do sample=`basename -s .fna $file`; echo -e "$sample\t$file" >> $outfile; done
else
	echo "mode not specified or recognized.Please specify mode in command line."
	echo "Choose from samples, illumina, trimmed, ncbi"
#	break
	exit 1

fi



# --------------
# remove undetermined 
# if [[ ! keep-undetermined = true && $dryrun != true ]]
if [[ ${dryrun} != true ]]
then
	sed -i '/Undetermined_/d' ${outfile}
fi


# check ---------------
if [[ ${check} == true ]]; then
    tail -n +2 ${outfile} | while read sample read1 read2; do
        if ! [[ ${read1} =~ "${sample}_" ]]; then
            echo "Error: $sample not contained in $read1"
            [[ "$0" = "$BASH_SOURCE" ]] && exit 1 || return 1
        fi
    done
    echo "Check successful"
fi


# ---------------


if [[ ${dryrun} != true ]]
then
	samplesfound=`tail -n +2 ${outfile} | wc -l `
	echo "Found $samplesfound samples"
fi



