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
weighted_network<-loadNetwork(network.file)


#create new forest class with data - both inheriting class info and additional
fObj <- basicFendR(network=weighted_network,
  featureData=gene.data,
  sampleOutcomeData=pheno.data,
  phenoFeatureData = target.data
 )

testDrugs=c('selumetinib',"sorafenib","vorinostat")

##TODO: leave out samples at a time, create features, then test.
fObj<-createNewFeaturesFromNetwork(fObj,testDrugs)

##i commented these so they won't run just yet, but should be run across all fendR objects
fObj<-buildModelFromEngineeredFeatures(fObj,testDrugs)

#score<-scoreDataFromModel(fObj)



##eventually move these to some other function.
##removes sample from feature data
removeSampleFromObject <- function(obj,sampleName){
  if(sampleName%in%unique(obj$sampleOutcomeData$Sample))
    obj$sampleOutcomeData<-subset(obj$sampleOutcomeData,Sample!=sampleName)
  else
    print(paste('Sample',sampleName,'not found in outcome data'))
  return(obj)
}

##performs loo cross validation by sample - should we move to other file?
looCrossValidation<-function(obj,testDrugs){
  #iterate over all cell lines
  all.samps<-intersect(obj$sampleOutcomeData$Sample,obj$featureData$Sample)

  vals<-sapply(all.samps,function(x){

    #subset out that data and re-assign original object
    test.data<-subset(obj$sampleOutcomeData,Sample==x)
    orig.test.features<-subset(obj$featureData,Sample==x)
    aug.test.features<-subset(obj$remappedFeatures,Sample==x)

    ##now remove from object
    newObj<-removeSampleFromObject(fObj,x)

    #create baseline model
    baselineModelObj<-buildModelFromOriginalFeatures(newObj,testDrugs)

    #get baseline predictions
    baselinePreds<-scoreDataFromModel(baselineModelObj,orig.test.features,test.data)

    #create new features
    augmentedObj<-createNewFeaturesFromNetwork(newObj,testDrugs)

    #create new model - this will replace the model list object in the class
    augmentedObj<-buildModelFromEngineeredFeatures(augmentedObj,testDrugs)

    updatedPreds<-scoreDataFromModel(baselineModelObj,aug.test.features,test.data)
    data.frame(select(baselinePreds,originalPred=Prediction,Actual),select(updatedPreds,augmentedPred=Prediction))


  })
  vals

}

##we can add some generic fendR methods as well, such as plotting, statistics, loo, etc.
