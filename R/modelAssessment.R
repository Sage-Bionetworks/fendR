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
#' @import plyr dplyr doMC
#' @return Not sure yet...
crossValidationCompare <- function(fendRObj,
  modelCall='lm',
  modelArgs=list(),
  testPheno=c(),
  sampleIndependent=TRUE){

  ##get a list of all samples
  all.samps<-intersect(fendRObj$sampleOutcomeData$Sample,fendRObj$featureData$Sample)

  if(length(testPheno)==0)
    testPheno <- unique(fendRObj$sampleOutcomeData$Phenotype)

  if(is.null(fendRObj$graph))
    fendRObj <- loadNetwork(fendRObj) ###only need to load graph once

  if(sampleIndependent)##if the samples are independent we can generate this once
    fendRObj<-createNewFeaturesFromNetwork(fendRObj,testPheno)

  ##FIX
  ##for now we assume that all models are assume independence between samples AND Drugs
  origMatrix<-originalResponseMatrix(fendRObj)
  engMatrix<-engineeredResponseMatrix(fendRObj)


  doMC::registerDoMC(2)
  #for each sample, leave one out
  vals<-plyr::laply(all.samps,function(x){

    #subset out that data and re-assign original object
    test.df<-dplyr::filter(fendRObj$sampleOutcomeData,Sample==x)%>%filter(Phenotype%in%testPheno)
    test.data<-test.df$Response
    names(test.data)<-test.df$Phenotype
    orig.test.features<-subset(fendRObj$featureData,Sample==x)
    #artificially expand to do acast
    otf<-orig.test.features$Value
    names(otf)<-orig.test.features$Gene

    aug.test.features<-subset(fendRObj$remappedFeatures,Sample==x)
    atf<-reshape2::acast(select(aug.test.features,Phenotype,Gene,Value),Phenotype~Gene)


    #build original and updated model
    orig.pred<-sapply(testPheno,function(p){
      mod.dat<-dplyr::filter(origMatrix,Sample!=x)%>%filter(Phenotype==p)%>%dplyr::select(-Phenotype,-Sample)
      orig.mod<-do.call(modelCall,args=list(formula='Response~.',data=mod.dat))
      predict(orig.mod,newdata=data.frame(t(otf)))[[1]]
    })

    mod.dat<-dplyr::filter(engMatrix,Sample!=x)%>%dplyr::select(-Phenotype,-Sample)
    eng.mod<-do.call(modelCall,args=list(formula='Response~.',data=mod.dat))
    eng.pred<-predict(eng.mod,data.frame(atf))


      df<-data.frame(OriginalPred=orig.pred,
          EngineeredPred=eng.pred,Phenotype=testPheno,
          Sample=rep(x,length(testPheno)),TrueValue=test.data)

      df
      },.parallel = TRUE)

    all.res<-do.call('rbind',vals)

    return(all.res)

}


plotModelResults <- function(modelingDataFrame){
  ##data frame has a

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


