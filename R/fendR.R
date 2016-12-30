## This contains the abstract functionality of the fendR
## might become semi-abstract if we move dataProcessing tools in this

######################################################################
# Create the fendR class

#' An S3 class to represent a fendR predictive network algorithm
#'
#' @param network A matrix (exact format TBD)
#' @param featureData a data.frame that contains rows representing genes and columns representing samples
#' @param phenoData a data.frame representing at least one column of phenotype and rows representing samples
#' @export
#' @return a fendR object
fendR<-function(network, featureData, phenoData){
  me <- list(network=network,
    featureData=featureData,
    phenoData=phenoData,
    remappedFeatures=featureData, ##default to original data for testing
    featureModel=NULL)
  class(me) <- append(class(me),"fendR")
  return(me)
}

######################################################################

######################################################################
# these are the generic methods to be implemented

#' Engineer Features from Network
#' \code{createNewFeaturesFromNetwork} takes the gene-based measurements and alters their #' score using a network
#' @param object That contains a data frame and network
#' @keywords
#' @export
#' @return data.frame representing new gene by sample matrix that is augmented
createNewFeaturesFromNetwork<-function(object){
  UseMethod('createNewFeaturesFromNetwork',object)
}

createNewFeaturesFromNetwork.fendR <- function(object){
  print("This method cannot be called on generic fendR class")
  return(object)
}

#' Builds predictive model from network-augmented feature
#' \code{buildModelFromEngineeredFeatures} takes the engineered features and creates a model based on an underlying
#' @param object That contains a phenotypic data
#' @param newFeatureSet from \code{createNewFeaturesFromNetwork}
#' @keywords
#' @export
#' @return an object of class \code{lm}
buildModelFromEngineeredFeatures <- function(object){
  UseMethod('buildModelFromEngineeredFeatures',object)
}

##now create warnin
#' Builds predictive model from network-augmented feature
#' \code{buildModelFromEngineeredFeatures} takes the engineered features and creates a model based on an underlying
#' @param object That contains a phenotypic data
#' @param newFeatureSet from \code{createNewFeaturesFromNetwork}
#' @keywords
#' @export
#' @return an object of class \code{lm}
buildModelFromEngineeredFeatures.fendR <- function(object){
  print('Building linear model from remapped features')
  over<-intersect(rownames(object$phenoData),colnames(object$remappedFeatures))
  all.mods<-apply(object$phenoData[over,],2,function(x){
    df<-data.frame(drug=x,t(object$remappedFeatures[,over]))
    mod<-lm(drug~.,data=df)
  })
  names(all.mods)<-colnames(object$phenoData)
  object$model <- all.mods
  return(object)


}

#' \code{scoreDataFromModel} takes the new model and predicts a phenotype from an input set
#' @param model
#' @param unseenFeatures
#' @keywords
#' @export
#' @return list of scores for each of the columns of the unseen feature data frame
scoreDataFromModel <- function(object, unseenFeatures){
  UseMethod('scoreDataFromModel',object)
}

scoreDataFromModel.fendR <- function(object, unseenFeatures){
  print("This method cannot be called on generic fendR class")
  return(rep(0.0),ncol(unseenFeatures))
}
