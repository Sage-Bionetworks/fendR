#this file creates examples that run fendR

library(fendR)

#load in network, files from system.file included in package

##should we load the data or not? seems like a waste of time at this point
gene.file<-system.file('CCLE_binary_mutation_matrix_ucscGenesFromCBioPortal.tsv',package='fendR')
gene.data<-loadSampleData(gene.file)
pheno.file<-system.file('CTRP_v20_AUC_vales_by_drug.tsv',package='fendR')
pheno.data<-loadPhenotypeData(pheno.file)
network.file<-'https://github.com/fraenkel-lab/OmicsIntegrator/raw/master/data/iref_mitab_miscore_2013_08_12_interactome.txt'
#network.data<-read.table(network.file,sep='\t')

#create new forest class with data - both inheriting class info and additional
fObj <- forestFendR(network=network.file,
  featureData=gene.data,
  phenoData=pheno.data,
  forestPath='../../OmicsIntegrator')


graphs<-createNewFeaturesFromNetwork(fObj)

fObj<-buildModelFromEngineeredFeatures(fObj)

score<-scoreDataFromModel(fObj)

##we can add some generic fendR methods as well, such as plotting, statistics, loo, etc.
