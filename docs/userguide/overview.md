# Overview

Welcome to the user manual of FooDme.
This documentation will guide you through your first steps with the pipeline.

## Aim

The aim of this workflow is the taxonomic assignement of paired-end Illumina amplicon sequencing reads.
It was developped with a focus of food and feed authenticity control.

## Workflow description

The workflow consists of three main steps that will be performed sequentially: reads pre-processing,
clustering, and taxonomic assignment.

### Input data

The pipeline expects paired-end Illumina reads, provided as paths in a tabular file.
See the [configuration help](configuration.md) for more details.

### Pre-processing

As a first analysis step, primers will be trimmed form the reads. By default primers are only matched on 
the 5' end of the reads. In some cases (e.g. Sequencing is longer than the amplicon length) one may want
to trim primers on the 3' end as well. This behaviour can be triggered in the [parameters](configuration.md). 

The reads will then be pre-processed for quality trimming on the 3' using a sliding window checking for minimal 
quality requirements.

### Clustering

FooDMe implements three different clustering strategies to choose from. Each has specific performances
that can be better suited to your specific needs.

#### Identity clustering

For this strategy, sequences are first dereplicated and ranked by abundance.
Each sequences is then comparaired to an itinially empty list of centroid and 
based on sequence similarity is either assigned as a new centroid or merged with 
the most similar existing centroid.

The degree of identity required for clustering the sequences can be freely set, 0.97 
being a commonly accepted value and 1.0 corresponding to a dereplication.

This approach is very good at smoothing out seqencing noise but can also results
in the clustering of highly similar sequences.

#### Dereplication

With the strategy, only identitcal reads will be clustered together.
This is effectively implemented as an identity clustering strategy with 100% identity.
This is more suitable for high sensitivity identification of amplicons but
the sequencing noise will not be filtered.

Due to the large number of clusters resulting form this strategy this is also the most 
computive intensive of the three.

#### Denoising

For the denoising strategy, the nucleotide substitution rate of the sequencing run 
will be modelled based on the available data. Using this model, reads can be corrected 
for sequencing-induced substitutions. This results in a typically low number of cluster, 
close to the biological reality of the sample.

### Taxonomic assignment

#### BLAST search and filtering

The representative sequences for each clusters are compaired to a user provided nucleotide
database using Basic Local Alignment Strategy (BLAST) and references satisfying specified 
similiraty critera are recovered.

Because this typically results in a large number of matching results (and taxa), the matches 
can be post-filtered based on their alignment scores.

#### Taxonomic consensus

As this process often results in a unclear mix of taxa, a consensus can be determined based 
on the underlying taxonomic hierachy and a minimal agreement level that can be freely set
between strict majority (0.51) or last-common ancestor (1.0).
