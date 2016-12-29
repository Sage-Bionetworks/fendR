
##this is the forest implemetnation of the fendR class

######################################################################
# Create the forestFendR class
#
# This is used to represent a forest-based implementation of the fendR framework

#' An S3 class to represent a Forest-based implementation of the fendR predictive
#' network algorithm
#' @inheritParams fendR
#' @param forestPath Path to python forest scripts
#' @export
#' @return forestFendR object
forestFendR<-function(network, featureData, phenoData, forestPath){
 me <-fendR(network,featureData,phenoData)
 me$forestPath <- forestPath
 class(me) <- append(class(me),'forestFendR')
 return(me)
}

######################################################################
# Set methods implemented by forestFendR

#' Engineer Features from Network
#' \code{createNewFeaturesFromNetwork} takes the gene-based measurements and alters their #' score using a network
#' @param object That contains a data frame and network
#' @keywords
#' @export
#' @return data.frame representing new gene by sample matrix that is augmented
createNewFeaturesFromNetwork.forestFendR<-function(object){
    ##summarize phenotypic features across samples
    geneSums<-rowSums(object$featureData)
    geneSums<-subset(geneSums,geneSums>0)
    #these become our nodeWeights for forest

    #write network, weight, config file to working directory
    prizeFileName=''
    netFileName=''
    confFile=''

    #then call forest.py with arguments
    cmd=paste('python ',object$forestPath,'/scripts/forest.py -p ',prizeFileName,' -e ',netFileName,
      ' -c ',confFile,sep='')

    #optimize arguments for number of trees and size.
    mu.range<-seq(0.001,0.004,by=0.001)
    beta.range<-seq()
    w.range<-seq()

    #iterate over all possible parameters

    ##read in files and compare number of nodes in graph

    #assign new features based on shortest path, using igraph


}

#internal helper function to load file to igraph
loadForestNetworkToIgraph<-function(fname){

}


#' Builds predictive model from network-augmented feature
#' \code{buildModelFromEngineeredFeatures} takes the engineered features and creates a model based on an underlying
#' @param object That contains a phenotypic data
#' @param newFeatureSet from \code{createNewFeaturesFromNetwork}
#' @keywords
#' @export
#' @return an object of class \code{lm}
buildModelFromEngineeredFeatures.forestFendR<-function(object){
  #start with a linear model? haven't tested this yet

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
scoreDataFromModel.forestFendR<-function(object){

}


