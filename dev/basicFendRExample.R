#this file creates examples that run fendR

library(fendR)

#load in network, files from system.file included in package

##should we load the data or not? seems like a waste of time at this point
gene.file<-system.file('CCLE_binary_mutation_matrix_ucscGenesFromCBioPortal.tsv',package='fendR')
gene.data<-loadSampleData(gene.file)
pheno.file<-system.file('CTRP_v20_AUC_vales_by_drug.tsv',package='fendR')
pheno.data<-loadPhenotypeData(pheno.file)

target.file<-system.file('CTRP_v20_drug_target_vals.tsv',package='fendR')
target.data<-loadTargetData(target.file)

network.file<-'https://github.com/fraenkel-lab/OmicsIntegrator/raw/master/data/iref_mitab_miscore_2013_08_12_interactome.txt'




#create new basicFendR class with data - both inheriting class info and additional
fObj <- basicFendR(networkFile=network.file,
  featureData=gene.data,
  sampleOutcomeData=pheno.data,
  phenoFeatureData = target.data
 )

#sampling 10 drugs
testDrugs=unique(fObj$phenoFeatureData$Phenotype)

testDrugs<-sample(testDrugs,3)

#these are the four functions we need
fObj<-loadNetwork(fObj)
fObj <- createNewFeaturesFromNetwork(fObj,testDrugs)

#origMatrix<-originalResponseMatrix(fObj,phenotype=testDrugs)
#engMatrix<-engineeredResponseMatrix(fObj,phenotype=testDrugs)


##we can add some generic fendR methods as well, such as plotting, statistics, loo, etc.
res<-crossValidationCompare(fObj,
  modelCall='lm',
  modelArgs=list(),
  testPheno=testDrugs,
  sampleIndependent=TRUE)
plotModelResults(res)
write.table('fendRtestResults.tsv',sep='\t',header=T,row.names=F)
