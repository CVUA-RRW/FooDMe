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
* Dobrovolny S, Blaschitz M, Weinmaier T, Pechatschek J, Cichna-Markl M, Indra A,  
Hufnagl P, Hochegger R. Development of a DNA metabarcoding method for the  
identification of fifteen mammalian and six poultry species in food. Food Chem.  
2019 Jan 30;272:354-361. doi: 10.1016/j.foodchem.2018.08.032.   
Epub 2018 Aug 9. PMID: 30309555.

* Preckel L, Bruenen-Nieweler C, Denay G, Petersen H, Cichna-Markl M, Dobrovolny S,  
Hochegger R. Identification of Mammalian and Poultry Species in Food and Pet Food  
Samples Using 16S rDNA Metabarcoding. Foods. 2021 Nov 20;10(11):2875.  
doi: 10.3390/foods10112875. PMID: 34829156; PMCID: PMC8620145.

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

