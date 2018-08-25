# Ongoing development

This directory contains the files and information used to develop the fendR package. This directory has been added to the [.Rbuildignore](../.Rbuildignore) file.

## Dataset formatting
The files in here are used to develop and test various aspects. We can create sub-directories as need be. So far this is limited to:
* [formatCcleData.R](formatCcleData.R): a script that pulls data from cBioPortal and pushes matrices to Synapse and to inst directory

## Preliminary analysis
We need to show that standard approaches can fail for drug identification
* ability to single out individual gene
* ability to target specific pathways
* ability to target specific proteins (via viper)

## Basic methodology
Describe here

## Tests of efficacy
These test break down into two facets:
* how well can we recover drugs that
  * we know behave differently across cell lines
  * _could_ behave differently across cell lines
* whether or not these drugs have been tested or not (i.e. whether or not we're
exploring the space of all drugs or just those that we can evaluate)

The results are captured in this table below:
|| Testing Known Drugs | Testing by Genotype|
|--- | --- | ---|
|Evaluating all drugs | [Completed](byDrugTestAllDrugs)| [In Process]()|
|Evaluating tested drugs | [In Process](byDrugOnTestedDrugs) | [In Process]()|

## Next steps
If the results seem to work well on a single dataset, does additional dataset
harmonization help?
