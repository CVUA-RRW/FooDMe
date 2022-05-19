# Welcome to the documentation of FooDMe

![CI](https://github.com/CVUA-RRW/FooDMe/workflows/CI/badge.svg?branch=master)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/CVUA-RRW/FooDMe)](https://github.com/CVUA-RRW/FooDMe/releases/latest)

FooDMe is a reproducible and scalable snakemake workflow for the analysis of 
DNA metabarcoding experiments, with a special focus on food and feed samples.

FooDMe processes paired-end Illumina reads in the following way:

- Primer trimming
- Quality filtering
- Sequence clustering based on identity clustering (OTUs), dereplication, or denoising (ASVs)
- Similiraty search cluster sequences in a user-provided dataabase
- Determine a taxonomic consensus for each cluster
- Outputs quality reports and results in human-readable formats

