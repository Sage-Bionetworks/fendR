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
fendR<-function(networkFile, featureData, phenoFeatureData, sampleOutcomeData,targetGenes=NULL,geneNorm=NA,responseNorm=NA){

  ##reduce genes by targetGenes if possible  - moved this from n3 fendR
  full.gene.set<-unique(featureData$Gene)
  num.genes.in.full.feature.space <- length(full.gene.set)
  reduced.gene.set<-full.gene.set
  if(!is.na(targetGenes) && !is.null(targetGenes)) {
    if(!any(targetGenes %in% full.gene.set)) {
      stop("None of target genes are in the featureData")
    }
    reduced.gene.set <- full.gene.set[full.gene.set %in% targetGenes]
    cat(paste0("Reducing feature space from ", num.genes.in.full.feature.space, " to ", length(reduced.gene.set), "\n"))
  } else {
    cat(paste0("Using all ", num.genes.in.full.feature.space, " genes in featureData.\n"))
    cat("Not sparsifying data or network\n")
  }


  ##remove phenotype data for which there are not two tests of data
  ##figure out which phenotypes have both feature data and outcome data
  phenos<-intersect(sampleOutcomeData$Phenotype,phenoFeatureData$Phenotype)
  all.phenos<-union(sampleOutcomeData$Phenotype,phenoFeatureData$Phenotype)
  print(paste("Found",length(phenos),'phenotypes that have feature data and outcome data out of',length(all.phenos)))

  ##figure out which samples have associated response data
  com.samps<-intersect(sampleOutcomeData$Sample,featureData$Sample)
  all.samps<-union(sampleOutcomeData$Sample,featureData$Sample)
  print(paste("Found",length(com.samps),'samples that have feature data and outcome data out of',length(all.samps)))

  #first create base object, then use helper functions to update
  me <- list(network=networkFile,
    featureData=featureData,
    phenoFeatureData=phenoFeatureData,
    sampleOutcomeData=sampleOutcomeData,
    targetGenes=reduced.gene.set,
    graph = NULL, #to be filled in using selectFeaturesFromNetwork function
    remappedFeatures=NULL # to be filled in using engineerFeaturesFromNetwork function
  )
  class(me) <- append(class(me),"fendR")

  #first normalize if there are functions populated
  if(!is.na(geneNorm))
      me<-normalizeFeatureValues(me,geneNorm)
  if(!is.na(responseNorm)){
    me<-normalizeResponse(me,responseNorm)
  }
  ##now do some reduction
  me<-reduceGenes(me,reduced.gene.set)
  me<-reducePhenos(me,phenos)
  me<-reduceSamples(me,com.samps)

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


#' Normalize response to drug
#' @description Applies `normFunction` to each cell line across the drug screening library
#' @param object A fendR object
#' @param normFunction a function that will normalize the vector
#' @keywrods normalization
#' @export
#' @return A fendR object with normalized response
normalizeResponse <-function(object,normFunction){
  UseMethod('normalizeResponse',object)

}

#' Normalize gene expressoin data
#' @description Applies `normFunction` to each cell line across the features/genes
#' @param object A fendR object
#' @param normFunction a function that will normalize the vector
#' @keywrods normalization
#' @export
#' @return A fendR object with normalized gene expression dataset
normalizeFeatureValues <-function(object,normFunction){
  UseMethod('normalizeFeatureValues',object)

}

#' Engineer Features from Network
#' @description Takes the igraph populated by \code{loadNetwork}
#' and propagates scores to remapped features
#' @param object fendR object after running \code{loadNetwork}
#' @keywords network propagation
#' @export
#' @return fendR object with the remappedFeature data.frame populated with four columns: Gene, Sample, Phenotype, Value
createNewFeaturesFromNetwork<-function(object, ...){
  UseMethod('createNewFeaturesFromNetwork',object)
}

#'
#'@export
reduceGenes <-function(object,gene.list){
  UseMethod('reduceGenes',object)
}

#' Reduce model to only consider a subset of genes
#' @export
reduceGenes.fendR<-function(object,gene.list){
  object$featureData=subset(object$featureData,Gene%in%gene.list)
  if(!is.null(object$remappedFeatures))
    object$remappedFeatures<-subset(object$remappedFeatures,Gene%in%gene.list)
  return(object)
}

reduceSamples<-function(object,samp.list){
  UseMethod("reduceSamples",object)
}
#' Reduce mode to only consider a set of samples
reduceSamples.fendR <- function(object,samp.list){
  object$featureData=subset(object$featureData,Sample%in%samp.list)
  object$sampleOutcomeData=subset(object$sampleOutcomeData,Sample%in%samp.list)
  return(object)
}

reducePhenos <-function(object,phenos){
  UseMethod('reducePhenos',object)
}
#' Reduce model to only consider a set of phenotypes
reducePhenos.fendR <-function(object,phenos){
  object$sampleOutcomeData=subset(object$sampleOutcomeData,Phenotype%in%phenos)
  object$phenoFeatureData=subset(object$phenoFeatureData,Phenotype%in%phenos)
  return(object)
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

originalSparseResponseMatrix <- function(object, phenotype=c(), ...){
  UseMethod('originalSparseResponseMatrix',object)
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
  dres<-tidyr::spread(unique(fres),Gene,value=Value,drop=TRUE)

  return(dres)

}


#' Normalize response to drug
#' @description Applies `normFunction` to each cell line across the drug screening library
#' @param object A fendR object
#' @param normFunction a function that will normalize the vector
#' @keywrods normalization
#' @export
#' @return A fendR object with normalized response
normalizeResponse.fendR <-function(object,normFunction){

  updated.data<-group_by(object$sampleOutcomeData,Sample)%>%mutate(normResponse=normFunction(Response))
  updated.data<-data.frame(updated.data)%>%select(-Response)%>%select(Sample,Phenotype,Response=normResponse)

  object$sampleOutcomeData <-updated.data
  object
}

#' Normalize gene expressoin data
#' @description Applies `normFunction` to each cell line across the features/genes
#' @param object A fendR object
#' @param normFunction a function that will normalize the vector
#' @keywrods normalization
#' @export
#' @return A fendR object with normalized gene expression dataset
normalizeFeatureValues.fendR <-function(object,normFunction){

  updated.data<-group_by(object$featureData,Sample)%>%mutate(normValue=normFunction(Value))
  updated.data<-data.frame(updated.data)%>%select(-Value)%>%select(Gene,Sample,Value=normValue)

  object$featureData <-updated.data
  object

}

#' Get re-engineered feature matrix
#' @description Grabs data from the re-engineered features and creates a model-ready response matrix
#' @keywords model fendR
#' @export
#' @return A response matrix to use for modeling with the formula Response~.
engineeredResponseMatrix <- function(object,phenotype=c()){
  UseMethod('engineeredResponseMatrix',object)
}

engineeredSparseResponseMatrix <- function(object, phenotype=c(), ...){
  UseMethod('engineeredSparseResponseMatrix',object)
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
