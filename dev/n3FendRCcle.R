## Run the n3 fendR class against CCLE data

suppressPackageStartupMessages(library(fendR))
suppressPackageStartupMessages(library(synapseClient))

synapseLogin()

## Read in CCLE mutation and drug response data.
gene.file <- system.file('CCLE_binary_mutation_matrix_ucscGenesFromCBioPortal.tsv', package='fendR')
gene.data <- loadSampleData(gene.file)
pheno.file <- system.file('CTRP_v20_AUC_vales_by_drug.tsv',package='fendR')
pheno.data <- loadPhenotypeData(pheno.file)

## Read in drug targets.
target.file <- system.file('CTRP_v20_drug_target_vals.tsv', package='fendR')
target.data <- loadTargetData(target.file)

## Download Yuanfang Guan's network (translated from mouse MGI ids to human Hugo symbols by
## aggregrating multiple columns to one using their mean)
network.file <- getFileLocation(synGet("syn8265096"))

target.genes <- unique(target.data$Gene)

# Create the n3 fendR object
fObj <- n3FendR(network = network.file,
  featureData = gene.data,
  sampleOutcomeData = pheno.data,
  phenoFeatureData = target.data,
  target.genes = target.genes
 )

testDrugs=c('selumetinib',"sorafenib","vorinostat")

res <- createNewFeaturesFromNetwork(fObj, testDrugs)

stop("stop")


#fObj<-buildModelFromEngineeredFeatures(fObj,testDrugs)




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
    print(paste('Removing sample',x,'to evaluate'))
    #subset out that data and re-assign original object
    test.data<-subset(obj$sampleOutcomeData,Sample==x)
    orig.test.features<-subset(obj$featureData,Sample==x)
    aug.test.features<-subset(obj$remappedFeatures,Sample==x)

    ##now remove from object
    newObj<-removeSampleFromObject(obj,x)

    #create baseline model
    baselineModelObj<-buildModelFromOriginalFeatures(newObj,testDrugs)

    #get baseline predictions
    baselinePreds<-scoreDataFromModel(baselineModelObj,orig.test.features,test.data)

    #create new features
    augmentedObj<-createNewFeaturesFromNetwork(newObj,testDrugs)

    #create new model - this will replace the model list object in the class
    augmentedObj<-buildModelFromEngineeredFeatures(augmentedObj,testDrugs)

    updatedPreds<-scoreDataFromModel(augmentedObj,aug.test.features,test.data)
    df<-data.frame(select(baselinePreds,originalPred=Prediction,Actual),select(updatedPreds,augmentedPred=Prediction))
    df$Drug<-rownames(df)
    df$Sample<-rep(x,nrow(df))
    df
  })
  vals

}

##we can add some generic fendR methods as well, such as plotting, statistics, loo, etc.
