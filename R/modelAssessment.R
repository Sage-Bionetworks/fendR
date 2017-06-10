##
## This file contains the tools to compare models after feature engineering
##

#' Calculate the Area Under the ROC Curve
#' @description Calculate the Area Under the ROC Curve
#' @param response  A vector of the response data
#' @param predicted.response  A vector of the predicted response.
#' @import pROC
#' @return  The Area Under the ROC Curve
calculate.auc <- function(response, predicted.response, ...) {
  suppressPackageStartupMessages(library("pROC"))
  pROC::auc(response=response, predictor=predicted.response, ...)
}

#' Calculate the ROC Curve
#' @description Calculate the ROC Curve
#' @param response  A vector of the response data
#' @param predicted.response  A vector of the predicted response.
#' @import pROC
#' @return  A list of class "roc"
calculate.roc <- function(response, predicted.response, ...) {
  suppressPackageStartupMessages(library("pROC"))
  pROC::roc(response=response, predictor=predicted.response, ...)
}


#' Assess fendR and modeling using cross validation
#' @description Use leave-one-out cross validation to assess the feature mapping and
#' modeling approach
#' @param fendRObj The object class
#' @param modelCall A String to call the model and predict
#' @param modelArgs a list of extra arguments to pass into do.call with model
#' @param testPheno a list of phenotypes to limit the test
#' @param sampleIndependent set to TRUE if your engineered features are independent of samples
#' @keywords INCOMPLETE
#' @export
#' @return Data frame with 5 columsn: Phenotype, Sample, TrueValue,OriginalPred,EngineeredPred
crossValidationCompare <- function(fendRObj,
  modelCall='lm',
  modelArgs=list(),
  testPheno=c(),
  k=10,
  sampleIndependent=TRUE){


  if(length(testPheno)==0)
    testPheno <- unique(fendRObj$sampleOutcomeData$Phenotype)

  if(is.null(fendRObj$graph))
    fendRObj <- loadNetwork(fendRObj) ###only need to load graph once

  if(sampleIndependent)##if the samples are independent we can generate this once
    fendRObj<-createNewFeaturesFromNetwork(fendRObj,testPheno)


  ##for now we assume that all models are assume independence between samples AND Drugs
  origMatrix<-originalResponseMatrix(fendRObj)
  print(paste("Dimensions of original matrix",paste(dim(origMatrix),collapse=',')))
  engMatrix<-engineeredResponseMatrix(fendRObj)
  print(paste("Dimensions of engineered matrix",paste(dim(engMatrix),collapse=',')))

  ##get a list of all samples

  all.samps <- intersect(origMatrix$Sample,engMatrix$Sample)
  testPheno<-testPheno[testPheno%in%fendRObj$phenoFeatureData$Phenotype]

  #doMC::registerDoMC()
  #for each sample, leave one out
  require(caret)
  kfold<-sapply(caret::createFolds(all.samps,k=k),function(x) all.samps[x])
#  kfold<-split(all.samps,base::sample(1:k,size=length(all.samps),replace=TRUE))
  vals<-plyr::llply(kfold,function(x){

    #subset out that data and re-assign original object
    test.df<-dplyr::filter(fendRObj$sampleOutcomeData,Sample%in%x)%>%filter(Phenotype%in%testPheno)
    test.data<-test.df$Response
    names(test.data)<-test.df$Phenotype
    orig.test.features<-subset(fendRObj$featureData,Sample%in%x)
    #artificially expand to do acast

    otf<-reshape2::acast(orig.test.features,Sample~Gene)
    aug.test.features<-subset(fendRObj$remappedFeatures,Sample%in%x)

    #build original model
    orig.pred<-data.frame(sapply(testPheno,function(p){
      mod.dat<-dplyr::filter(origMatrix,!Sample%in%x)%>%dplyr::filter(Phenotype==p)%>%dplyr::select(-Phenotype)
      rownames(mod.dat)<-mod.dat$Sample
      mod.dat<-dplyr::select(mod.dat,-Sample)
      orig.mod<-do.call(modelCall,args=list(formula='Response~.',data=mod.dat))
      predict(orig.mod,newdata=data.frame(otf))
    }))

  orig.pred$Sample=rownames(orig.pred)
  orig.df<-gather(orig.pred,"Phenotype","Response",1:(ncol(orig.pred)-1))
  orig.df$PredType=rep("OriginalPrediction",nrow(orig.df))
  ##build engineered model
    eng.pred<-data.frame(sapply(testPheno,function(p){
      mod.dat<-dplyr::filter(engMatrix,Sample!=x)%>%dplyr::filter(Phenotype==p)%>%dplyr::select(-Phenotype)
      rownames(mod.dat)<-mod.dat$Sample
      mod.dat<-dplyr::select(mod.dat,-Sample)

      atf<-reshape2::acast(select(dplyr::filter(aug.test.features,Phenotype==p),Gene,Sample,Value),Sample~Gene,value.var='Value')

      eng.mod<-do.call(modelCall,args=list(formula='Response~.',data=mod.dat))
      predict(eng.mod,newdata=data.frame(atf))

    }))

    eng.pred$Sample=rownames(eng.pred)
    eng.df<-gather(eng.pred,"Phenotype","Response",1:(ncol(eng.pred)-1))
    eng.df$PredType=rep('EngineeredPrediction',nrow(eng.df))

    test.df$PredType=rep('TrueValues',nrow(test.df))
    full.df<-rbind(test.df,orig.df,eng.df)
    full.df
      },.parallel = TRUE)

    all.res<-do.call('rbind',vals)

    return(all.res)

}

#' Visualirze fendR and modeling using cross validation
#' @description plots results of LOO CV code above
#' @param modelingDataFrame output of \code{crossValidationCompare}
#' @keywords
#' @import ggplot2 tidyr
#' @export
#' @return image

plotModelResults <- function(modelingDataFrame,prefix=''){
  ##data frame is result from LOO
  require(tidyr)
  require(ggplot2)
  res.df<-tidyr::spread(modelingDataFrame,PredType,Response)

  origCor=cor(res.df$OriginalPrediction,res.df$TrueValues,use='pairwise.complete.obs')
  engCor=cor(res.df$EngineeredPrediction,res.df$TrueValues,use='pairwise.complete.obs')

  byDrug<-res.df%>%dplyr::group_by(Phenotype)%>%summarise(original=cor(OriginalPrediction,TrueValues,use='pairwise.complete.obs'),engineered=cor(EngineeredPrediction,TrueValues,use='pairwise.complete.obs'))

  bySample<-res.df%>%dplyr::group_by(Sample)%>%summarise(original=cor(OriginalPrediction,TrueValues,use='pairwise.complete.obs'),engineered=cor(EngineeredPrediction,TrueValues,use='pairwise.complete.obs'))

  corDrugDf<-byDrug%>%tidyr::gather(Features,Correlation,2:3)
  corSampDf<-bySample%>%tidyr::gather(Features,Correlation,2:3)
  pdf(paste0(prefix,'modelResults.pdf',sep=''))
  p<-  ggplot2::ggplot(corDrugDf)+ggplot2::geom_point(ggplot2::aes(x=Phenotype,y=Correlation,col=Features))+ggtitle('Correlation across samples by drug')+ theme(axis.text.x = element_text(angle = 90, hjust = 1))

  print(p)

  p<-  ggplot2::ggplot(corSampDf)+ggplot2::geom_point(ggplot2::aes(x=Sample,y=Correlation,col=Features))+ggtitle('Correlation across drugs by sample')+
    theme(axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())
  print(p)

  ##now print the correlations themselves
  p<-ggplot(data=res.df)+geom_point(aes(x=TrueValues,y=EngineeredPrediction,col=Phenotype))
  print(p)
  p<-ggplot(data=res.df)+geom_point(aes(x=TrueValues,y=OriginalPrediction,col=Phenotype))
  print(p)
  dev.off()

}

#' Assess fendR and modeling using random weights
#' @description Compare performance of fendR with random weights added to network/features
#' @param fendRObj The object class
#' @param testPheno a list of phenotypes to limit the test
#' @keywords INCOMPLETE
#' @export
#' @return Not sure yet...
compareModelToRandom <-function(fendRObj,testPheno){

}


