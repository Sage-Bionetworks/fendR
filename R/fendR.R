## This contains the abstract functionality of the fendR
## might become semi-abstract if we move dataProcessing tools in this

######################################################################
# Create the fendR class

#' An S3 class to represent a fendR predictive network algorithm
#'
#' @param network A network or path to network, depending on the underlying class
#' @param featureData a data.frame that contains rows representing genes and columns representing samples
#' @param sampleOutcomeData a data.frame representing at least one column of phenotype and rows representing samples
#' @param phenoFeatureData a data.frame where rows represent genes and columns represent a relationship between phenotype and gene
#' @export
#' @return a fendR object
fendR<-function(networkFile, featureData, phenoFeatureData, sampleOutcomeData){
  me <- list(network=networkFile,
    featureData=featureData,
    phenoFeatureData=phenoFeatureData,
    sampleOutcomeData=sampleOutcomeData,
    graph = NULL, #to be filled in using selectFeaturesFromNetwork function
    remappedFeatures=NULL # to be filled in using engineerFeaturesFromNetwork function
  )
  class(me) <- append(class(me),"fendR")
  return(me)
}

######################################################################

######################################################################
# these are the generic function to be implemented

#' Load network from file or path
#' @description Takes the network handle and loads relevant features from network into iGraph object
#' @param object A fendR object
#' @param testPheno A set of phenotypes to limit the analysis scope
#' @keywords network
#' @export
#' @return A fendR object with a data frame with three columns: Gene, Phenotype, NetworkScore
loadNetwork <- function(object){
  UseMethod('loadNetwork',object)
}


#' Engineer Features from Network
#' @description Takes the igraph populated by \code{loadNetwork}
#' and propagates scores to remapped features
#' @param object fendR object after running \code{loadNetwork}
#' @keywords network propagation
#' @export
#' @return fendR object with the remappedFeature data.frame populated with four columns: Gene, Sample, Phenotype, Value
createNewFeaturesFromNetwork<-function(object,testDrugs){
  UseMethod('createNewFeaturesFromNetwork',object)
}


#' Get model-ready matrix from original features
#' @description return the original feature data in a response matrix for modeling
#' @param object a fendrObject
#' @keywords
#' @export
#' @return A response matrix to use for modeling with the formula 'Response~.'
originalResponseMatrix <- function(object, phenotype=c(), ...){
  UseMethod('originalResponseMatrix',object)
}

#' Get model-ready matrix from original features for fendR class
#' @description return the original feature data in a response matrix for modeling
#' @param object a fendrObject
#' @keywords
#' @export
#' @return A response matrix to use for modeling with the formula 'Response~.'
originalResponseMatrix.fendR <- function(object,phenotype=c()){
  if(length(phenotype)==0)
    phenotype <- unique(object$sampleOutcomeData$Phenotype)

  fres<-dplyr::inner_join(object$featureData,subset(object$sampleOutcomeData,Phenotype%in%phenotype),by='Sample')%>%dplyr::select(Sample,Gene,Value,Response,Phenotype)
  dres<-tidyr::spread(fres,Gene,value=Value,drop=TRUE)

  return(dres)

}

#' Get re-engineered feature matrix
#' @description Grabs data from the re-engineered features and creates a model-ready response matrix
#' @keywords model fendR
#' @export
#' @return A response matrix to use for modeling with the formula Response~.
engineeredResponseMatrix <- function(object,phenotype=c()){
  UseMethod('engineeredResponseMatrix',object)
}


#' Helper function to remove feature data or sample data from any fendR object
#' @description Un-exported for now
#' @param fendRObj The object class
#' @param sampleName
#' @keywords leaveOneOut
#' @return updated object
removeSampleOrPhenoFromObject <- function(obj,sampleName='',phenoName=''){
  if(sampleName!=''&&sampleName%in%unique(obj$sampleOutcomeData$Sample)){
    obj$sampleOutcomeData<-subset(obj$sampleOutcomeData,Sample!=sampleName)
#    obj$featureData <- subset(obj$featureData,Sample!=sampleName)
  }

  if(phenoName!='' && phenoName%in%unique(obj$sampleOutcomeData$Phenotype)){
    obj$sampleOutcomeData<-subset(obj$sampleOutcomeData,Phenotype!=phenoName)
  #  obj$phenoFeatureData<-subset(obj$phenoFeatureData,Phenotype!=phenoName)
  }
  #  print(paste('Sample',sampleName,'not found in outcome data'))
  return(obj)
}
