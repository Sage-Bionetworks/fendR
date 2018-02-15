##various tests of Cross dataset functionality,
source("crossDatasetPrediction.R")


#test random forests
testCDPrf<-function(){
  ccle<-buildFendRFromDataset('CCLE',geneNorm=zscore,responseNorm=zscore)
  ccle$sampleOutcomeData$Phenotype<-tolower(ccle$sampleOutcomeData$Phenotype)
  ccle$phenoFeatureData$Phenotype<-tolower(sapply(ccle$phenoFeatureData$Phenotype,function(x) unlist(strsplit(x,split='_'))[1]))

  sanger<-buildFendRFromDataset('Sanger',geneNorm=zscore,responseNorm=zscore)
  sanger$sampleOutcomeData$Phenotype<-tolower(sapply(sanger$sampleOutcomeData$Phenotype,function(x) unlist(strsplit(x,split='_'))[1]))
  sanger$phenoFeatureData$Phenotype<-tolower(sapply(sanger$phenoFeatureData$Phenotype,function(x) unlist(strsplit(x,split='_'))[1]))

  ##now we have to clean up the sanger dataset here
  newsangervals<-mutate(sanger$sampleOutcomeData,SampPhen=paste(Sample,Phenotype,sep='_'))
  newresp<-sapply(unique(newsangervals$SampPhen),function(x){
    mv<-which(newsangervals$SampPhen==x)
    if(length(mv)>1)
      return(mean(newsangervals$Response[mv],na.rm=T))
    else
      return(newsangervals$Response[mv])
  })

  newsanger<-data.frame(SampPhen=names(newresp),Response=newresp)%>%separate(SampPhen,c("Sample","Phenotype"),sep='_')
  sanger$sampleOutcomeData<-newsanger

  library(randomForest)
  res<-predict(ccle,sanger,numDrugs=10,modelCall='randomForest',modelArgs=list(na.action=na.omit))
  pl<-plotModelResults(res,'ccleToSanger10drugsrf')

  library(glmnet)
  res<-predict(ccle,sanger,numDrugs=10,modelCall='glmnet',predictArgs=list(s=0.001,type='response'))#,modelArgs=list(na.action=na.omit))
  pl<-plotModelResults(res,'ccleToSanger10drugsglmnet')

  res<-predict(ccle,sanger,numDrugs=10,modelCall='glmnet',modelArgs=list(alpha=0.5),predictArgs=list(s=0.001,type='response'))#,modelArgs=list(na.action=na.omit))
  pl<-plotModelResults(res,'ccleToSanger10drugsglmnetAlpha5')

  res<-predict(ccle,sanger,numDrugs=10,modelCall='glmnet',modelArgs=list(alpha=0.8),predictArgs=list(s=0.001,type='response'))#,modelArgs=list(na.action=na.omit))
  pl<-plotModelResults(res,'ccleToSanger10drugsglmnetAlpha8')

}
testCDP<-function(){

  ccle<-buildFendRFromDataset('CCLE',geneNorm=zscore,responseNorm=zscore)
  ccle$sampleOutcomeData$Phenotype<-tolower(ccle$sampleOutcomeData$Phenotype)
  ccle$phenoFeatureData$Phenotype<-tolower(sapply(ccle$phenoFeatureData$Phenotype,function(x) unlist(strsplit(x,split='_'))[1]))

  sanger<-buildFendRFromDataset('Sanger',geneNorm=zscore,responseNorm=zscore)
  sanger$sampleOutcomeData$Phenotype<-tolower(sapply(sanger$sampleOutcomeData$Phenotype,function(x) unlist(strsplit(x,split='_'))[1]))
  sanger$phenoFeatureData$Phenotype<-tolower(sapply(sanger$phenoFeatureData$Phenotype,function(x) unlist(strsplit(x,split='_'))[1]))

  ##now we have to clean up the sanger dataset here
  newsangervals<-mutate(sanger$sampleOutcomeData,SampPhen=paste(Sample,Phenotype,sep='_'))
  newresp<-sapply(unique(newsangervals$SampPhen),function(x){
    mv<-which(newsangervals$SampPhen==x)
    if(length(mv)>1)
      return(mean(newsangervals$Response[mv],na.rm=T))
    else
      return(newsangervals$Response[mv])
  })

  newsanger<-data.frame(SampPhen=names(newresp),Response=newresp)%>%separate(SampPhen,c("Sample","Phenotype"),sep='_')
  sanger$sampleOutcomeData<-newsanger

  res<-predict(ccle,sanger,numDrugs=10)
  pl<-plotModelResults(res,'ccleToSanger10drugslm')

  res2<-predict(sanger,ccle,numDrugs=10)
  pl2<-plotModelResults(res2,'SangerToCcle10drugslm')

  ntap<-buildFendRFromDataset('NTAP',geneNorm=zscore,responseNorm=zscore)

  ntap$sampleOutcomeData$Phenotype<-tolower(ntap$sampleOutcomeData$Phenotype)
  res<-predict(sanger,ntap,numDrugs=10)

  pl<-plotModelResults(res,'SangerToNTAP10drugslm')
}


