# Workflow configuration

## Configuration

To configure the workflow, modifiy `config/config.yaml` according 
to your needs.

For specific details about the parameters 
check the documentation.

### Optimized configuration for 16S birds and mammals metabarcoding

A configuration file with optimized parameters (with no warranty!) for 
the analysis of 16S metabarcoding experiments of birds and mammals is 
included in this folder. 
**You will still need to modify the indicated paths to match your own folder structure.**

See for reference:

Denay, G.; Preckel, L.; Petersen, H.; Pietsch, K.; Wöhlke, A.; Brünen-Nieweler, C. 
Benchmarking and Validation of a Bioinformatics Workflow for Meat Species Identification Using 16S rDNA Metabarcoding. 
Foods 2023, 12, 968. https://doi.org/10.3390/foods12050968 

## Sample sheet

Add your samples to `config/samples.tsv` or use the 
the helper script from Bundes Insitut für Risikobewertung in 
`ressources/creat_sampleSheet.sh`.

## Reference table

If using the benchmark mode, sdd information for your reference samples in the `reference.tsv` template.
Note that proportions are given as fractions in the [0, 1] interval.

## Parameter space exploration configuration

To configure a parameter space exploration analysis, modify the file 
`config_paramspace.yaml` according to your needs.

For specific details about the parameters 
check the documentation.

