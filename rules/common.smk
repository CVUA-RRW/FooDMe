import pandas as pd
import os
import time
import subprocess


sample_path = config["samples"]

samples = pd.read_csv(sample_path, 
                      index_col="sample", 
                      sep = "\t", 
                      engine="python")
samples.index = samples.index.astype('str', copy=False) 


pipe_log = os.path.join(os.getcwd(), "PIPELINE_STATUS")


def git_version():
    try:
        __version__ = subprocess.check_output(
                        ["git", "describe", "--always --tags"], 
                        cwd= workflow.basedir).strip().decode("utf-8")
    except subprocess.CalledProcessError:
        __version__ = "Unknown"
    finally:
        return(__version__)


def get_local_time():
    return time.asctime(time.localtime(time.time()))