# Installation

Below is a detialled installation guide.
Beginners are encouraged to follow this guide step by step. 

!!! note 

    The examples below will describe an installation in the `/home/user` folder.
    It is assumed that raw sequencing data are stored in the `raw_data` folder in the user´s home directory.

## Requirements

Becausesome of the tools used by the workflow only work in a UNIX environment,
Windows OS is *not* supported.

The ressources needed to run the pipeline vary according to the nucleotide database used.
In general be sure to have enought hard-drive space to store the databases and the raw data plus 
space to store the results. 

The workflow will remove tempory intermediate files on the fly to minimize the memory footprint.

!!! note

    Depending on the analysis mehtod, sequencing depth, and samples complexity, the persistent output
    should be below 100MB per sample.

Additionally the BLAST step will require an amount of virtual memory equivalent to the size of 
the database. This can be either hard drive space or RAM.

A minimal configuration for working with the BLAST nt database is therefore:

- 500 Gb Hard drive
- 8 Gb RAM

Increasing the number of cores will considerably speed up the workflow by taking advantage of 
parallelization.

The workflow can therefore run on a medium range laptop, even within a Virtual Machine
emulating Linux. 

An internet connection will be nescessary for the first run of the pipeline. Successive runs can be 
performed without an internet connection

## Conda

Snakemake makes intensive use of the environment manager [conda](https://docs.conda.io/en/latest/).
There are many different distributions of conda to choose from, each with their advantages or inconvenients.
For a new installation we recommend using the minimalistic distribution [miniconda](https://docs.conda.io/en/latest/miniconda.html).

For an installation guide, see the [Bioconda documentation](https://bioconda.github.io/user/install.html#set-up-channels), specifically steps 1 and 2:

!!! quote

    1. Install conda

    Bioconda requires the conda package manager to be installed. If you have an Anaconda Python installation, you already have it. Otherwise, the best way to install it is with the Miniconda package. The Python 3 version is recommended.

    On MacOS, run:

    ```bash
    curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    sh Miniconda3-latest-MacOSX-x86_64.sh
    ```

    On Linux, run:

    curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    sh Miniconda3-latest-Linux-x86_64.sh

    Follow the instructions in the installer. If you encounter problems, refer to the Miniconda documentation. You can also join our Gitter channel to ask other users for help.

    2. Set up channels

    After installing conda you will need to add the bioconda channel as well as the other channels bioconda depends on. It is important to add them in this order so that the priority is set correctly (that is, conda-forge is highest priority).

    The conda-forge channel contains many general-purpose packages not already found in the defaults channel.

    ```bash
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    ```

The newer versions on Snakemake also require that conda be set in strict channel priority mode.
This can be done by changing the conda configuration:

```bash
conda config --set channel_priority strict
```

The dependency solver of conda being notoriously slow and helpless in front 
of complex environments, it is advised to supplement the conda installation 
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

!!! info 

    Without git, you can download the repository manually and unpack the archive locally.

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

!!! note 

    The created conda environment can be reused in later runs by specifying the argument
    `--conda-prefix ~/conda-envs` when executing snakemake. This will save the environment 
    creation time.

!!! warning 

    The `.tests` folder might be hidden on your file explorer. If you don´t see it,
    enable `view hidden file` in the options.

Feel free to explore the files produced by this first run in the `.tests` folder. 
More details about the [output](results.md) and the [use of snakemake](run.md) will be given in the later sections
of this guide.
