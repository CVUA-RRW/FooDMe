# Installation

A detailled installation guide. 
To help new users, examples below will describe an installation in the `/home/user` folder.

## Requirements

Because of the tools used in the workflow only work in a UNIX environment,
Windows OS is *not* supported.

The ressources needed to run the pipeline vary according to the nucleotide database used.
In general be sure to have enought hardrive space to store the database and you raw data plus 
space to store the results. 

The workflow will remove tempory intermediate files on the fly to minimize the memory footprint.

Additionally the BLAST step will require an amount of virtual memory equivalent to the size of 
the database. This can be either hard drive space or RAM.

A minimal configuration for working with the BLAST nt database is therefore:

- 500 Gb Hard drive
- 8 Gb RAM

Increasing the number of cores will considerably speed up the workflow by taking advantage of 
parallelization.

The workflow can therefore run on a medium range laptop, even within a Virtual Machine
emulating Linux. 

## Conda

Snakemake makes intensive use of the environment manager [conda](https://docs.conda.io/en/latest/).
There are many different distributions of conda to choose from, each with their advantages or inconvenients.
For a new installation we recommend using the minimalistic distribution [miniconda](https://docs.conda.io/en/latest/miniconda.html).

The dependency solver of conda being notoriously slow and helpless in front 
of complex environment, it is advised to supplement the conda installation 
with a better solver called [mamba](https://github.com/mamba-org/mamba):

```bash
conda install -n base -c conda-forge mamba
```

## Download the repository

Install [git](https://git-scm.com/):

```bash
mamba install -n base -c conda-forge git
```

Then get a copy of the workflow:

```bash
cd ~ 
git clone https://github.com/CVUA-RRW/FooDMe.git
```

> **_NOTE:_** Without git, you can download the repository manually and upack the archive locally.

Having git installed will later allow you to get the latest version of the workflow
as well keep the older versions archived for reproducibility and traceability purposes.

For example update the workflow with:

```bash
cd ~/FooDMe
git pull origin
```

## Set up a conda environment to run the workflow

Running the workflow will require snakemake to be installed in the current environment.
To avoid future conflicts between software versions, it is recommended to create 
a new environment to execute snakemake:

```bash
mamba create -n snakemake -c conda-forge -c bioconda snakemake
```

The snakemake environment can then be toggled on and off with:

```bash 
conda activate snakemake
conda deactivate snakemake
```

## Test the installation

The repository comes with a minimal dataset allowing to run a quick test of the installation.
Running this example is also a good occasion to initialize all the software dependencies 
that will be reused on later runs.

In this examples we will store the workflow´s enviroments under `~/conda-envs`.

```bash
cd ~/FooDMe
conda activate snakemake
snakemake --use-conda --conda-prefix ~/conda-envs --cores 1
```

Snakemake will the start creating the nescessary conda environments (this can take a few minutes)
before analyzing the three examples samples in the `.tests` folder.

> **_NOTE:_** The `.tests` folder might be hidden on your file explorer. If you don´t see it,
> enable `view hidden file` in the options.

Feel free to explore the files produced by this first run in the `.tests` folder. 
More details about the output and the use of snakemake will be given in the later sections
of this guide.
