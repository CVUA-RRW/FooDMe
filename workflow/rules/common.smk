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


# Loading and validationg benchmark reference --------
reference_path = config["benchmark"]["reference"]
reference = pd.read_csv(
    reference_path, index_col="sample", sep="\t", engine="python"
)
validate(reference, schema="../schema/reference.schema.yaml")
reference.index = reference.index.astype("str", copy=False)
# Get union of reference and samples to use for benchmarking
samples_set = set(samples.index)
reference_set = set(reference.index)
benchmark_index = list(samples_set.union(reference_set))


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
