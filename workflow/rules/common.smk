import pandas as pd
import os
import time
from snakemake.utils import validate


# Pipeline setup --------------------------------------
version = open(os.path.join(workflow.basedir, "..", "VERSION"), "r").read()
pipe_log = os.path.join(os.getcwd(), "PIPELINE_STATUS")


# Validating config ----------------------------------
validate(config, schema="../schema/config.schema.yaml")


# Loading and validating samples ---------------------
sample_path = config["samples"]
samples = pd.read_csv(sample_path, index_col="sample", sep="\t", engine="python")
validate(samples, schema="../schema/samples.schema.yaml")
samples.index = samples.index.astype("str", copy=False)


# General puprose functions --------------------------
def get_local_time():
    return time.asctime(time.localtime(time.time()))


# Input functions ------------------------------------
def get_fastq(wildcards, read_pair="fq1"):
    return samples.loc[(wildcards.sample), [read_pair]].dropna()[0]


def get_mask():
    if config["taxid_filter"] == "None":
        return "common/nomask"
    else:
        return "common/taxid_mask.txt"


def get_blocklist():
    if config["blocklist"] == "None":
        return "common/noblock"
    elif config["blocklist"] == "extinct":
        return os.path.join(workflow.basedir, "..", "data", "blocklist.txt")
    else:
        return config["blocklist"]


def get_acc_blocklist(wildcards):
    if config["seq_blocklist"] == "None":
        return f"{wildcards.sample}/taxonomy/{wildcards.sample}_blast_report.tsv"
    else:
        return f"{wildcards.sample}/taxonomy/{wildcards.sample}_blast_report_prefiltered.tsv"
