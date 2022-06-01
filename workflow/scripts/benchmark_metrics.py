#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import pandas as pd
from sklearn.metrics import (
    precision_score, 
    recall_score, 
    f1_score, 
    average_precision_score,
)


sys.stderr = open(snakemake.log[0], "w")


def main(confmat, output, sample):
    conf_table = pd.read_csv(confmat, sep="\t")
    
    # get classification metrics
    precision = precision_score(conf_table['expected'], conf_table['predicted'])
    recall = recall_score(conf_table['expected'], conf_table['predicted'])
    fscore = f1_score(conf_table['expected'], conf_table['predicted'])
    prauc = average_precision_score(conf_table['expected'], conf_table['pred_ratio'])
    
    with open(output, "w") as fout:
        fout.write("\t".join(["Sample", "Precision", "Recall", "F1 score", "Average precision"]))
        fout.write("\n")
        fout.write("\t".join([sample, str(precision), str(recall), str(fscore), str(prauc)]))
        fout.write("\n")


if __name__ == '__main__':
    main(
        confmat=snakemake.input['confmat'],
        output=snakemake.output['metrics'],
        sample=snakemake.params['sample'],
    )