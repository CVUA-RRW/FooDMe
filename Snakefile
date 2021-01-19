import pandas as pd
import os 
import subprocess
import time
import shutil
from snakemake.utils import min_version


# Settings ------------------------------------------------------------------------------------------------------------------

min_version("5.10")


shell.executable("bash")


configfile: os.path.join(workflow.basedir, "tests", "config", "config.yaml")


workdir: config["workdir"]


sample_path = config["samples"]


samples = pd.read_csv(sample_path, index_col="sample", sep = "\t", engine="python")
samples.index = samples.index.astype('str', copy=False) # in case samples are integers, need to convert them to str


pipe_log = os.path.join(os.getcwd(), "PIPELINE_STATUS")


# Functions -----------------------------------------------------------------------------------------------------------------

def git_version():
    try:
        __version__ = subprocess.check_output(["git", "describe", "--always"], cwd= workflow.basedir).strip().decode("utf-8")
    except subprocess.CalledProcessError:
        __version__ = "Unknown"
    finally:
        return(__version__)

# Input rule ----------------------------------------------------------------------------------------------------------------
 
rule all:
    input: 
        expand("{sample}/reports/{sample}_report.html", sample = samples.index),
        "reports/report.html",

# Includes ------------------------------------------------------------------------------------------------------------------


include: "rules/trimming.smk"
include: "rules/vsearch.smk" if config["cluster"]["method"] == "otu" else "rules/dada2.smk" 
include: "rules/blast.smk"
include: "rules/reports.smk"


# Workflow ------------------------------------------------------------------------------------------------------------------

onstart:
    print("\nYou are using FooDMe version: {}".format(git_version()))
    with open(pipe_log, 'a') as f:
        f.write("[" + time.asctime(time.localtime(time.time())) + "]: Pipeline started\n")

onsuccess:
    try:
        for logfile in os.listdir(".snakemake/log/"):
            shutil.move(os.path.join(".snakemake/log", logfile), "logs")
        shutil.rmtree(".snakemake", ignore_errors=True)
    except:
        # if not executing .snakemake from workdir, the log file will be in execution directory
        # as far as I know, there is now way to access this form here
        pass
    
    print("\nWorkflow finished, no error")
    with open(pipe_log, 'a') as f:
        f.write("[" + time.asctime(time.localtime(time.time())) + "]: Pipeline succesfully finished\n")
    
onerror:
    for logfile in os.listdir(".snakemake/log/"):
        shutil.move(os.path.join(".snakemake/log", logfile), "logs")
    shutil.rmtree(".snakemake", ignore_errors=True)
    
    print("\nAn error occured, please consult the latest log file in {}".format(os.path.join(config["workdir"], ".snakemake", "log")))
    with open(pipe_log, 'a') as f:
        f.write("[" + time.asctime(time.localtime(time.time())) + "]: Pipeline stopped on error\n")
