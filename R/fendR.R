## This contains the abstract functionality of the fendR
## might become semi-abstract if we move dataProcessing tools in this

######################################################################
# Create the forestFendR class

#' An S4 class to represent a fendR predictive network algorithm
#'
#' @slot network A matrix (exact format TBD)
#' @slot featureData a data.frame that contains rows representing genes and columns representing samples
#' @slot phenoData a data.frame representing at least one column of phenotype and rows representing samples
fendR <- setClass(
  "fendR",

  slots = c(
    network = "matrix",
    featureData = "data.frame",
    phenoData = "data.frame"
  ), ##TODO: add slots for additional feature sets

  #prototype = c() ##TODO set this

  validity = function(object)
  {
    #two checks (so far)
    #1 make sure that the overlap between features and networks is greater than some threshold (e.g. 100)
    #2 make sure dimensions of feature data match the dimensions of the phenoData
  }
)
######################################################################

######################################################################
# Set generic methods that must be implemented by final class

#' Engineer Features from Network
#' \code{createNewFeaturesFromNetwork} takes the gene-based measurements and alters their #' score using a network
#' @param fendrObject That contains a data frame and network
#' @keywords
#' @export
#' @return data.frame representing new gene by sample matrix that is augmented
setGeneric(name="createNewFeaturesFromNetwork",
  def=function(fendrObject){
    standardGeneric("createNewFeaturesFromNetwork")
  })

#' Builds predictive model from network-augmented feature
#' \code{buildModelFromEngineeredFeatures} takes the engineered features and creates a model based on an underlying
#' @param fendrObject That contains a phenotypic data
#' @param newFeatureSet from \code{createNewFeaturesFromNetwork}
#' @keywords
#' @export
#' @return an object of class \code{lm}
setGeneric(name="buildModelFromEngineeredFeatures",
  def=function(fendrObject,newFeatureSet){
    standardGeneric("buildModelFromEngineeredFeatures")
    }
  )

#' \code{scoreDataFromModel} takes the new model and predicts a phenotype from an input set
#' @param model
#' @param unseenFeatures
#' @keywords
#' @export
#' @return list of scores for each of the columns of the unseen feature data frame
setGeneric(name="scoreDataFromModel",
  def=function(model,unseenFeatures){
    standardGeneric("scoreDataFromModel")
    }
  )
