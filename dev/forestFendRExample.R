#this file creates examples that run fendR

library(fendR)

#load in network, files from system.file included in package
gene.data<-loadSampleData(system.file('CCLE_binary_mutation_matrix_ucscGenesFromCBioPortal.tsv',package='fendR'))
pheno.data<-loadPhenotypeData(system.file('CTRP_v20_AUC_vales_by_drug.tsv',package='fendR'))
network.data<-read.table('https://github.com/fraenkel-lab/OmicsIntegrator/raw/master/data/iref_mitab_miscore_2013_08_12_interactome.txt',sep='\t')

#create new forest class with data - both inheriting class info and additional
fObj <- forestFendR(network=network.data,
  featureData=gene.data,
  phenoData=pheno.data,
  forestPath='../../OmicsIntegrator')


fObj<-createNewFeaturesFromNetwork(fObj)

fObj<-buildModelFromEngineeredFeatures(fObj)

score<-scoreDataFromModel(fObj)

##we can add some generic fendR methods as well, such as plotting, statistics, loo, etc.
