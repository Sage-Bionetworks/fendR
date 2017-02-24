##
## This file contains the tools to compare models after feature engineering
##


#' Assess fendR and modeling using cross validation
#' @description Use leave-one-out cross validation to assess the feature mapping and
#' modeling approach
#' @param fendRObj The object class
#' @param modelCall A String(?) to call the model and predict
#' @param testPheno a list of phenotypes to limit the test
#' @param featuresIndependent set to TRUE if your engineered features are independent of samples
#' @keywords INCOMPLETE
#' @export
#' @return Not sure yet...
crossValidationCompare <- function(fendRObj,
  modelCall='lm',
  testPheno=c(),
  samplesIndependent=TRUE){

  ##get a list of all samples
  all.samps<-intersect(fendRObj$sampleOutcomeData$Sample,fendRObj$featureData$Sample)

  if(length(testPheno)==0)
    testPheno <- unique(fendRObj$sampleOutcomeData$Phenotype)

  if(is.null(fendRObj$graph))
    fendRObj <- loadNetwork(fendRObj) ###only need to load graph once

  if(samplesIndependent)##if the samples are independent we can generate this once
    fendRObj<-createNewFeaturesFromNetwork(fendRObj,testPheno)

  ##FIX
  ##for now we assume that all models are assume independence between samples AND Drugs
  origMatrix<-originalResponseMatrix(fendRObj)
  engMatrix<-engineeredResponseMatrix(fendRObj)


  #for each sample, leave one out
  vals<-lapply(all.samps,function(x){


    #subset out that data and re-assign original object
    test.df<-filter(fendRObj$sampleOutcomeData,Sample==x)%>%filter(Phenotype%in%testPheno)
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
      mod.dat<-filter(origMatrix,Sample!=x)%>%filter(Phenotype==p)%>%select(-Phenotype,-Sample)
 #   zvars<-which(apply(mod.dat,2,var)==0)
#    if(length(zvars)>0)
#      mod.dat<-mod.dat[,-zvars]
      orig.mod<-do.call(modelCall,list(formula='Response~.',data=mod.dat))
      predict(orig.mod,newdata=data.frame(t(otf)))[[1]]
    })

    mod.dat<-filter(engMatrix,Sample!=x)%>%select(-Phenotype,-Sample)
 #   zvars<-which(apply(mod.dat,2,var)==0)
#    if(length(zvars)>0)
#      mod.dat<-mod.dat[,-zvars]
    eng.mod<-do.call(modelCall,list(formula='Response~.',data=mod.dat))
    eng.pred<-predict(eng.mod,data.frame(atf))


      df<-data.frame(OriginalPred=orig.pred,
          EngineeredPred=eng.pred,
          Sample=rep(x,length(testPheno)),TrueValue=test.data)

      df
      })

    all.res<-do.call('rbind',vals)
  vals

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


