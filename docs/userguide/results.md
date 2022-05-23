# Viewing the results

The analysis produces a number of files containing quality reports, tables and graphs, 
as well as some of the processed data at different steps of the analysis.

!!! note

    By default, large files are removed during the analysis
    use snakemake´s `--notemp` to keep all files.

## Folder structure

The output folder follows the following structure:

```
workdir/
|- common/             \\ Taxonomy processing files
|- logs/
|   |- all/            \\ Logs for steps processing aggregated results
|   |- common/         \\ Logs for taxonomy produces
|   |- sample_name/    \\ Sample-wise log files
|- reports/            \\ All analysis reports
|- sample_name/        \\ Sample-wise processed data are located
    |- denoising/      \\    here, splitted in the relevant folders
    |- krona/
    |- taxonomy
    |- trimmmed
    |- reports         \\ Analysis reports for individual samples
```

## Report

All the reports are stored as tab-delimited text files in the `reports` folder.
Additionaly an HTML report is produced that contains the information of the individual files.
Tables in this report are interactive and can be filtered, columns moved or removed, and 
the tables can be exported in excel, csv, and pdf formats.

### Panel descriptions

The HTML report is divided in several sections reflecting the different processing steps.

Tables can be sorted using the arrows in each header cell, filtered using the fields below the header,
columns can be moved by drag and drop, and column visibility toggled using the buttoin above the table.
Additional buttons above the table allow to copy the entire table to the clipboard, print it, or export it 
in popular formats such as Excel or PDF.

#### Quality summary

This presents a succint summary of the most important quality statistics of the analyis.

- Sample: sample identifier 
- Q30 rate: proportion of bases above Q30 after quality processing
- Insert size peak: Estimated mean size of the insert
- Read number: total input reads
- Pseudo-reads: number of merged reads
- Reads in ASV/OTU: number of clustered reads
- ASV/OTU number: number of clusters
- Assigned reads: number of reads assigned to a taxa
- Rank-consensus: number of clusters with a consensus at the given rank
- No match: number of clusters that could not be successfully assigned

#### Trimming statistics

Summary of the primer and quality trimming processes.

- Sample: sample identifier
- Total raw reads: number of inputs reads (forward + reverse)
- Total reads after primer trimming: number of reads ´left after primer trimming. Reads where primers cannot be found are discarded.
- No primer found: percent of reads where primers could not be found
- Total bases before quality trim: number of nucleotide in the input reads
- Total reads after qualtiy trim: number of pre-processed reads available for clustering
- Total bases after quality trim: number of nucleotide after the pre-processing
- Q20 after: propoertion of bases above Q20 after quality processing
- Q30 rate: proportion of bases above Q30 after quality processing
- Duplication rate: number of duplicated reads
- Insert size peak: Estimated mean size of the insert
- links: Link to the fastp report

#### Read filtering statistics (OTU only)

Summary of the quality controls of reads prior to clustering.

- Sample: sample identifier
- Total reads: number of input reads for the error-correction process
- Pseudo-reads:  number of reads resulting form merging process
- Merging failure: proportion of reads that could not be merged
- Pseudo-reads PF: number of merged-reads passing size filters
- Discarded reads: proportion of merged-reads not passing size filter
- Unique sequences: Number and proportion of distinct sequences

#### Clustering statistics (OTU only)

Summary of the OTU clustering process.

- Sample: sample identifier
- Unique sequences: Number and proportion of distinct sequences
- Clusters: Number of clusters passing minimal population requirement
- Discarded clusters: proportion of reads and clusters not passing minimal population requirement
- Non-chimeric clusters: number of cluster not flagged as chimeric
- Chimeras: proportion of clusters and reads flagged as chimeric and removed
- Pseudo-reads clustered: number ands proportion of pseudo-reads successfully clustered

#### Denoising statistics (ASV only)

Summary of the denoising and read merging steps.

- Sample: sample identifier
- Total reads: number of input reads for the error-correction process
- Filtered reads: number of reads passing size filters
- Discarded reads: proportion of reads not passing size filter
- Denoised R1/2: number of successfully corrected reads in the respective orientation
- Merged: number of reads resulting form merging process
- Merging failure: proportion of reads that could not be merged
- ASV: number of clustered
- Size-filtered ASV: number of clusters passing minimal population requirements
- Discarded ASVs: proportion of clusters and reads not passing size requirements
- Non-chimeric ASV: number of cluster not flagged as chimeric
- Chimeras: proportion of clusters and reads flagged as chimeric and removed
- Reads in ASVs: final number and proportion of reads succesfully clustered

#### Taxonomic assignment

Summary of the taxonomic assignment of individual cluster.

- Sample: sample identifier
- Query: cluster identifier
- Count: cluster population
- Blast hits: number of BLAST results for this cluster
- Best bit-score: sequence alignment score of the best BLAST result
- Lowest bit-score: sequence alignment score of the worst BLAST result
- Bit-score threshold: minimal sequence alignment score for result to be accepted
- Saved Blast hits: number of sequences passing the sequence alignment threshold
- Consensus: scientific name of the determined taxonomic consensus
- Rank: rank of the taxonomic consensus
- Taxid: taxonomic identifier of the consensus
- Disambiguation: scientific names and proportions of all saved Blast hits used for taxonomic consensus
- blast_report: link to raw BLAST results
- filtered report: link to BLAST report after alignment quality filtering

#### Taxonomic assignment statistics

Sample-wise summarized assignemnt report.

- Sample: sample identifier
- Query: number of clusters
- Unknown sequences: number and proportion of sequnces that did not yield BLAST results
- Rank consensus: Number and proportion of clusters assigned at the given rank

#### Metabarcoding results

Sample composition summary.

- Sample: sample identifier
- Consensus: scientific name of the determined taxonomic consensus
- Rank: rank of the taxonomic consensus
- Taxid: taxonomic identifier of the consensus
- Count: number of reads assigned to the taxon
- Disambiguation: scientific names and proportions of all saved Blast hits used for taxonomic consensus
- Percent of total: proportion of reads assigned to this taxon
- Percent of assigned: proportion of the assigned reads (excluding those with no consensus) assigned to this taxon

#### Graphical results

A graphical overview of the sample compositions. 
Several options can be changed and a snapshot downloaded.

#### Versions

Summary of the tools and databases used and their version or modification dates.