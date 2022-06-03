# Benchmark mode

FooDMe offers a benchmarking module that can be used to
automatically compare the anaylsis results to a theoritical 
sample composition and compute performance indicators.

This can be used directly for method validation or verification
or to compare the performance of the pipeline with different 
parameter sets.

## Benchmarking parameters

The benchmarking module requires three additional parameters:

- `reference`: Path to a table containing the sample composition information containing the followuing fields:
    - `sample`: Sample name
    - `taxid`: Taxonomic identifier
    - `proportion`: Expected fraction of the sample made up from this component, in the interval [0, 1]
- `threshold`: A minimal quantity for a component to be considered a 'true' result, given in percent.
- `target_rank`: A maximum taxonomic rank for a component to be considered a 'true' result

!!! warning
    
    Because of limitaitons of the NCBI taxonomiy classification, only
    Linnaean ranks are supported. This means you can choose only from 
    the following (case-sensitive) options:
    
    - species
    - genus
    - family
    - order
    - class
    - phylum
    - kingdom

## Usage

The basic usage only differs to the basic analysis by specifying the `benchmark` 
target directly after the `snakemake` call:

```bash
conda activate snakemake
snakemake benchmark \
  --use-conda --conda-prefix ~/conda-envs --cores 1 \
  --configfile ~/FooDMe/config/myconfig.yaml
```

This will either run the normal analysis followed by the benchmarking module,
or just the benchmarking module if the analysis already exists.

## Benchmarking results

The benchmarking results will be saved in the `benchmarking` subfolder, including 
tables and an interactive HTML report.
Sample specific metrics will also be saved in each sample's folder, in a specific subfolder.

### Yield

The yield is given as the amount of reads left in percent of total input reads for each of the three main 
analytical steps:

- Merged reads correspond to the yield of read pre-processing and merging
- Clustered reads corresponds to the yield of all steps up to the cluster (OTU or ASV) generation
- Assigned reads corresponds to the final yield of the analysis

!!! note
    
    For ASV, the merged reads yield includes the denoising steps

### Metrics

After filtering the analysis results with the given concentration threshold and
maximal rank, the results are compared to a "true" composition and several metrics are 
calculated to measure the reliability of the analysis.

Metrics are calculated for each sample and an aggregated metric is also calculatedthat considers 
the whole dataset as a single sample.

#### Classification metrics

Each result is categorised as either True-positive (TP), False-positive (FP), or False-Negative (FN)
depending on wether it is considered true and expected.

!!! note
    
    Because the number of true negative is roughly the number of taxid in the database
    and therefore several orders of magnitude higher than the number of true positives, 
    we do not calculate metrics reliying on the negative reuslt, such as specificity or ROC.

The following metrics are then calculated as:

- precision: part of true positives in all results predicted to be true. It is the reciprocal of the False-positve rate 

$$
P=\frac{TP}{TP+FP}
$$

- recall: part of real positives predicted as true. It is the reciprocal of the False-negative rate.

$$
R=\frac{TP}{TP+FN}
$$

- F1 score: harmonic mean of the precision and recall.

$$
F1=2\cdot\frac{P \cdot R}{P + R}
$$

- Average precision: summarizes the variation of precision and recall across a range of $n$ thresholds 
  (here concentration thresholds are used). It is an approximation of the area under the precision-recall curve.

$$
AP=\sum{n}(R_n - R_{n-1})P_n
$$

!!! info
    
    While the F1 score provides a direct measure of the analysis performance as it ran, 
    the average precision reflects the overall performance independently of the thrshold. 
    Combinations of both may be used to determine if the concentration threshold might be optimized.

#### Quantification metrics

Two quantification metrics are used here. The 'Distance' metric corresponds to the Euclidian distance
of the predicted and true results. It is a reflection of how far away are the results from the 
expected composition.
We also calculate the 'Mean Absolute Error' as the average value of the difference between predicted and expected concentrations for each samples:

$$
MAE= \sum_{n=1}^{N} |expected-predicted| \cdot \frac{1}{N}
$$

!!! note
    
    While false negative results account for 0% composition, false positive results 
    are only indirectly (through their contribution to the total composition) taken into
    account in the quantification metrics.

### Confusion matrix

The confusion table used to calculate the above mentionned metrics is given as a succint summary 
of the results with the following information:

- `Taxid`: the expected or determined taxonomic identifier
- `match_rank`: the rank to which a predicted result was matched to its expected value, given in the Linnaean taxonomical ranks.
- `predicted`/`expected`: `1`or `0` reprensenting true or false values
- `pred_ratio`: the predicted amount of this component in the sample
- `exp_ratio`: the expected amount of this component in the sample
