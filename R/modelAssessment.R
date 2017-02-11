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
crossValidationCompare <- function(fendRObj,modelCall='lm',testPheno=c(),featuresIndependent=TRUE){

  ##get a list of all samples
  all.samps<-intersect(obj$sampleOutcomeData$Sample,obj$featureData$Sample)

  if(is.null(fendRObj$graph))
    fendRObj <- loadNetwork(fendRObj) ###only need to load graph once

  if(featuresIndependent)
    fendRObj<-createNewFeaturesFromNetwork(fendRObj)

  #for each sample, leave one
  vals<-sapply(all.samps,function(x){

    print(paste('Removing sample',x,'to evaluate'))


    #subset out that data and re-assign original object
    test.data<-subset(fendRObj$sampleOutcomeData,Sample==x)
    orig.test.features<-subset(fendRObj$featureData,Sample==x)
    aug.test.features<-subset(fendRObj$remappedFeatures,Sample==x)

    testObj<-removeSampleFromObject(fendRObj,x)

    ##if features depend on samples we need to re-engineer?
    if(!featuresIndependent){
      testObj<-createNewFeaturesFromNetwork(testObj,testDrugs)
    }

    origMatrix<-originalResponseMatrix(testObj)
    engMatrixList<-engineeredResponseMatrix(testObj)

    #build original and updated model
    orig.mod<-do.call(modelCall,list(formula='Response~.',data=origMatrix))
    orig.pred<-predict(orig.mod,data=test.data)

    eng.preds<-lapply(engMatrixList,function(engMatrix){
      eng.mod<-do.call(modelCall,list(formula='Response~.',data=engMatrix))
      predict(eng.mod,data=aug.test.data)

      })

    #predictions on test data
    updatedPreds<-scoreDataFromModel(augmentedObj,aug.test.features,test.data)
    df<-data.frame(select(baselinePreds,originalPred=Prediction,Actual),select(updatedPreds,augmentedPred=Prediction))
    df$Drug<-rownames(df)
    df$Sample<-rep(x,nrow(df))
    df
  })
  vals

}

#' Helper function to remove feature data
#' @description Un-exported for now
#' @param fendRObj The object class
#' @param sampleName
#' @keywords leaveOneOut
#' @return updated object
removeSampleFromObject <- function(obj,sampleName){
  if(sampleName%in%unique(obj$sampleOutcomeData$Sample))
    obj$sampleOutcomeData<-subset(obj$sampleOutcomeData,Sample!=sampleName)
  else
    print(paste('Sample',sampleName,'not found in outcome data'))
  return(obj)
}

