# Welcome to the documentation of FooDMe
![CI](https://github.com/CVUA-RRW/FooDMe/workflows/CI/badge.svg?branch=master)
[![GitHub Clones](https://img.shields.io/badge/dynamic/json?color=success&label=Cloners&query=uniques&url=https://gist.githubusercontent.com/gregdenay/02b5545a991e1a51c423422e56f5500f/raw/clone.json&logo=github)](https://github.com/MShawon/github-clone-count-badge)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/CVUA-RRW/FooDMe)](https://github.com/CVUA-RRW/FooDMe/releases/latest)
[![DOI](https://zenodo.org/badge/296584559.svg)](https://zenodo.org/badge/latestdoi/296584559)

![logo](logo.png){ align=right }

FooDMe is a reproducible and scalable snakemake workflow for the analysis of 
DNA metabarcoding experiments, with a special focus on food and feed samples.

---

## Overview

FooDMe processes paired-end Illumina reads in the following way:

- Primer trimming
- Quality filtering
- Sequence clustering based on identity clustering (OTUs), dereplication, or denoising (ASVs)
- Similarity search cluster sequences in a user-provided database
- Determine a taxonomic consensus for each cluster
- Outputs quality reports and results in human-readable formats

---

## Quick links

- First time? Check the [User Guide](userguide/overview.md).
- For a history of the releases see [Changelog](about/changelog.md).
- A list of [Frequently Asked Questions](faqs.md).
- Get help [here](help.md)!

--- 

## Citation

If you use FooDMe fo reasearch please cite the project using the 
DOI at the top of this page.

---
## Connex tools

**BaRCoD:** A snakemake workflow aiming at recovering and analyzing barcodes for metabarcoding experiments.
Using a set of primers, finds possible amplicon in the database (or a taxonomic subset thereof) 
and performs pairwise alignements of barcodes. Available form [Github](https://github.com/CVUA-RRW/BaRCoD).

**TaxidTools:** A Python library to load, process and work with taxonomy definitions.
Documentation on the [homepage](https://cvua-rrw.github.io/taxidTools/).

