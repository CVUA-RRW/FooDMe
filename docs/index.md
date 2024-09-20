![logo](logo.png){ align=right width="200" }

# Welcome to the documentation of FooDMe

![CI](https://github.com/CVUA-RRW/FooDMe/workflows/CI/badge.svg?branch=master)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/CVUA-RRW/FooDMe)](https://github.com/CVUA-RRW/FooDMe/releases/latest)
[![DOI](https://zenodo.org/badge/296584559.svg)](https://zenodo.org/badge/latestdoi/296584559)
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)


**Support for FooDMe has ended ! Check out [FooDMe2](https://bio-raum.github.io/FooDMe2)**

FooDMe2 is a complete rewrite of FooDMe. We designed FooDMe 2 to be more flexible and take away some of the complexity encountered in FooDMe 1. This not only concerns the installation procedure, which is vastly streamlined now, but also the process of configuring and starting individual analysis runs. The new implementation also maes it easier to deploy, maintain, and to add additonal functionalities in the future.


---

## Overview

FooDMe processes paired-end Illumina reads in the following way:

- Primer trimming
- Quality filtering
- Sequence clustering based on identity clustering (OTUs), dereplication, or denoising (ASVs)
- Similarity-search cluster sequences in a user-provided database
- Determine a taxonomic consensus for each cluster sequence
- Outputs quality reports and results in human-readable formats

---

## Quick links

- First time? Check the [User Guide](userguide/overview.md).
- For a history of the releases see [Changelog](about/changelog.md).
- A list of [Frequently Asked Questions](faqs.md).
- Get help [here](help.md)!

--- 

## Citation

If you use FooDMe for research please cite us!

The source code is deposited on Zenodo and citable using the DOI at the top of this page (Last release)

We also published a benchmarking study for 16S Metabarcoding of meat products:

Denay, G.; Preckel, L.; Petersen, H.; Pietsch, K.; Wöhlke, A.; Brünen-Nieweler, C. 
Benchmarking and Validation of a Bioinformatics Workflow for Meat Species Identification Using 16S rDNA Metabarcoding. 
Foods 2023, 12, 968. https://doi.org/10.3390/foods12050968 

---

## Connex tools

**BaRCoD:** A snakemake workflow aiming at recovering and analyzing barcodes for metabarcoding experiments.
Using a set of primers, finds possible amplicon in the database (or a taxonomic subset thereof) 
and performs pairwise alignements of barcodes. Available from [Github](https://github.com/CVUA-RRW/BaRCoD).

**TaxidTools:** A Python library to load, process and work with taxonomy definitions.
Documentation on the [homepage](https://cvua-rrw.github.io/taxidTools/).

---

## About us 

FooDMe is developed by Grégoire Denay at the Chemical and Veterinary Investigation Office Rhein-Ruhr-Wupper, 
an official public laboratory in the field of consumer health protection.

More about [CVUA-RRW](https://www.cvua-rrw.de/).
