
##this is a basic implementation of the fendR class to demonstrate how to create a class
##TODO: discuss moving some of these features to the parent class.

######################################################################
# Create the basicFendR class
#
# This is used to represent a basicimplementation of the fendR framework

#' An S3 class to represent a basic implementation of the fendR predictive
#' network algorithm
#' @inheritParams fendR
#' @export
#' @return basicFendR object
basicFendR<-function(network, featureData, phenoFeatureData,sampleOutcomeData){
 me <-fendR(network, featureData, phenoFeatureData,sampleOutcomeData)
 class(me) <- append(class(me),'basicFendR')
 return(me)
}

######################################################################
# Set methods implemented by basicFendR



#' Engineer Features from Network
#' \code{createNewFeaturesFromNetwork} takes the phenotype scores and updates the gene scores based on phenotype
#' @param object That contains a data frame and network
#' @keywords
#' @export
#' @return list of gene features for each phenotype/drug response
createNewFeaturesFromNetwork.basicFendR<-function(object){
    ##figure out which phenotypes have both feature data and outcome data

    ##for each phenotype, for each gene
    ##update from featureData the score by shortest weighted path to target genes

    ##Reduction strategy: remove any genes that don't change across drugs?
}


#' \code{scoreDataFromModel} takes the new model and predicts a phenotype from an input set
#' @param model
#' @param unseenFeatures
#' @keywords
#' @export
#' @return list of scores for each of the columns of the unseen feature data frame
scoreDataFromModel.basicFendR<-function(object){

}


