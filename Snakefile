import pandas as pd
import os 
import subprocess
import time
import shutil


# Settings ------------------------------------------------------------------------------------------------------------------

shell.executable("bash")


workdir: config["workdir"]


pipe_log = os.path.join(config["workdir"], "PIPELINE_STATUS")

samples = pd.read_csv(config["samples"], index_col="sample", sep = "\t", engine="python")
samples.index = samples.index.astype('str', copy=False) # in case samples are integers, need to convert them to str

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


include: "rules/trimming.rule"
include: "rules/vsearch.rule" if config["cluster"]["method"] == "otu" else "rules/dada2.rule" 
include: "rules/blast.rule"
include: "rules/reports.rule"


# Workflow ------------------------------------------------------------------------------------------------------------------

onstart:
    print("\nYou are using FooDMe version: {}".format(git_version()))
    with open(pipe_log, 'a') as f:
        f.write("[" + time.asctime(time.localtime(time.time())) + "]: Pipeline started\n")

onsuccess:
    for logfile in os.listdir(".snakemake/log/"):
        shutil.move(os.path.join(".snakemake/log", logfile), "logs")
    shutil.rmtree(".snakemake", ignore_errors=True)
    
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
