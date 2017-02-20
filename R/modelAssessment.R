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
    test.data<-subset(fendRObj$sampleOutcomeData,Sample==x)
    orig.test.features<-subset(fendRObj$featureData,Sample==x)
    otf<-orig.test.features$Value
    names(otf)<-orig.test.features$Gene

    aug.test.features<-subset(engMatrix,Sample==x)
    atf<-aug.test.features$Value
    names(atf)<-aug.test.features$Gene

    #build original and updated model
    orig.mod<-do.call(modelCall,list(formula='Response~.',data=origMatrixList[[p]]))
    orig.pred<-predict(orig.mod,data=otf)

    eng.mod<-do.call(modelCall,list(formula='Response~.',data=engMatrixList[[p]]))
    eng.pred<-predict(eng.mod,data=atf)


      df<-data.frame(OriginalPred=orig.pred,EngineeredPred=eng.pred)
      df$Pheno=rep(p,nrow(df))
      df$LeftOutSample=rep(x,nrow(df))
      df$Vals<-test.data
      df
      })
    all.res<-do.call('rbind',vals)



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


