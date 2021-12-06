import pandas as pd
import os
import time
import subprocess
from snakemake.utils import validate


# Pipeline setup --------------------------------------
version = open(os.path.join(workflow.basedir, "..", "VERSION"), "r").read()
pipe_log = os.path.join(os.getcwd(), "PIPELINE_STATUS")


# Validating config ----------------------------------
validate(config, schema="../schema/config.schema.yaml")


# Loading and validation samples ---------------------
sample_path = config["samples"]
print(os.getcwd())
print(sample_path)
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
    if config["blast"]["taxid_filter"] == "None":
        return "common/nomask"
    else:
        return "common/taxid_mask.txt"


def get_blocklist():
    if config["blast"]["blocklist"] == "None":
        return "common/noblock"
    elif config["blast"]["blocklist"] == "extinct":
        return os.path.join(workflow.basedir, "..", "data", "blocklist.txt")
    else:
        return config["blast"]["blocklist"]


# Wrangling functions ------------------------------
def concatenate_uniq(entries):
    s = "; ".join(entries.to_list())
    df = pd.DataFrame(
        [e.rsplit(" (", 1) for e in s.split("; ")], columns=["name", "freq"]
    )  # parenthesis in names
    df["freq"] = df["freq"].str.replace(")", "", regex=False).astype(float)
    # Aggreagte, normalize, and sort
    tot = df["freq"].sum()
    df = df.groupby("name").apply(lambda x: x.sum() / tot)
    df = df.sort_values(by=["freq"], ascending=False)
    # Format as string
    uniq = df.to_dict()["freq"]
    uniq = [f"{name} ({round(freq, 2)})" for name, freq in uniq.items()]
    return "; ".join(uniq)
