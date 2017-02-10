## This contains the abstract functionality of the fendR
## might become semi-abstract if we move dataProcessing tools in this

######################################################################
# Create the fendR class

#' An S3 class to represent a fendR predictive network algorithm
#'
#' @param network An iGraph network
#' @param featureData a data.frame that contains rows representing genes and columns representing samples
#' @param sampleOutcomeData a data.frame representing at least one column of phenotype and rows representing samples
#' @param phenoFeatureData a data.frame where rows represent genes and columns represent a relationship between phenotype and gene
#' @export
#' @return a fendR object
fendR<-function(network, featureData, phenoFeatureData, sampleOutcomeData){
  me <- list(network=network,
    featureData=featureData,
    phenoFeatureData=phenoFeatureData,
    sampleOutcomeData=sampleOutcomeData,
    remappedFeatures=featureData, ##default to original data for testing
    featureModel=NULL)
  class(me) <- append(class(me),"fendR")
  return(me)
}

######################################################################

######################################################################
# these are the generic methods to be implemented


#' Engineer Features from Network
#' \code{createNewFeaturesFromNetwork} takes the gene-based measurements and alters their
#' score using the list of networks
#' @param object That contains a data frame and network
#' @keywords
#' @export
#' @return data.frame representing new gene by sample matrix that is augmented
createNewFeaturesFromNetwork<-function(object,testDrugs){
  UseMethod('createNewFeaturesFromNetwork',object)
}

#createNewFeaturesFromNetwork.fendR <- function(object,testDrugs){
#  print("This method cannot be called on generic fendR class")
#  return(object)
#}

#' Builds predictive model from network-augmented feature
#' \code{buildModelFromEngineeredFeatures} takes the engineered features and creates a model based on an underlying
#' @param object That contains a phenotypic data
#' @param newFeatureSet from \code{createNewFeaturesFromNetwork}
#' @keywords
#' @export
#' @return an object of class \code{lm}
buildModelFromEngineeredFeatures <- function(object,testDrugs){
  UseMethod('buildModelFromEngineeredFeatures',object)
}


#' Builds predictive model from network-augmented feature
#' \code{buildModelFromEngineeredFeatures} takes the engineered features and creates a model based on an underlying
#' @param object That contains a phenotypic data
#' @param testDrug name of drug to select from dataset to test in smaller context
#' @keywords
#' @export
#' @import dplyr
#' @return an object of class \code{fendR} with a list of \code{lm} objects
buildModelFromEngineeredFeatures.fendR <- function(object,testDrugs=NULL){
  suppressPackageStartupMessages(library(dplyr))
  print('Building linear model from remapped features')


   over<-intersect(object$sampleOutcomeData$Phenotype,object$remappedFeatures$Phenotype)
  if(!is.null(testDrugs))
    over<-testDrugs

  all.mods<-lapply(over,function(x){
      print(paste('Creating model for',x))
      mod.df<-dplyr::inner_join(subset(object$sampleOutcomeData,Phenotype==x)%>%dplyr::select(Sample,Response),subset(object$remappedFeatures,Phenotype==x),by="Sample")%>%dplyr::select(Sample,Gene,Value,Response)


        res<-tidyr::spread(mod.df,Gene,value=Value)
        res<-res[,-which(colnames(res)=='Sample')]
        mod<-lm(Response~.,data=res)
      mod
  })
  names(all.mods)<-over
  object$featureModel <- all.mods
  return(object)

}

#' Builds predictive model from original features
#' \code{buildModelFromOriginalFeatures} takes the original features and creates a model based on an underlying
#' @param object That contains a phenotypic data
#' @param testDrug name of drug to select from dataset to test in smaller context
#' @keywords
#' @export
#' @import dplyr
#' @return an fendR object with a list of \code{lm} objects
buildModelFromOriginalFeatures.fendR <- function(object,testDrugs=NULL){
  library(dplyr)
  print('Building linear model from original features')


  over<-unique(object$sampleOutcomeData$Phenotype)
  if(!is.null(testDrugs))
    over<-testDrugs

  #let's do a basic reduction here to test for genes who are not differentially mutated across samples
  gene.var<-object$featureData%>%dplyr::group_by(Gene)%>%dplyr::summarize(Variance=var(Value))
 # print(head(gene.var))

  #nonzero genes
  nz.genes<-gene.var$Gene[which(gene.var$Variance!=0)]
  print(paste('Reducing gene dataset to',length(nz.genes),'features that vary across all samples'))

  red.df<-subset(object$featureData,Gene%in%nz.genes)
  all.mods<-lapply(over,function(x){
    print(paste('Creating model for',x))

    mod.df<-dplyr::inner_join(subset(object$sampleOutcomeData,Phenotype==x),red.df,by="Sample")%>%select(Sample,Gene,Value,Response)
    res<-tidyr::spread(mod.df,Gene,value=Value)
    res<-res[,-which(colnames(res)=='Sample')]
    mod<-lm(Response~.,data=res)
    mod
  })
  names(all.mods)<-over
  object$featureModel <- all.mods
  return(object)
}


#' Builds predictive model from original features
#' \code{buildModelFromOriginalFeatures} takes the original features and creates a model based on an underlying
#' @param object That contains a phenotypic data
#' @keywords
#' @export
#' @return an fendR object with a list of \code{lm} objects
buildModelFromOriginalFeatures <- function(object,testDrugs){
  UseMethod('buildModelFromOriginalFeatures',object)

}

#' \code{scoreDataFromModel} takes the new model and predicts a phenotype from an input set
#' @param model
#' @param unseenFeatures
#' @param unseenResponse
#' @keywords
#' @export
#' @return list of scores for each of the columns of the unseen feature data frame
scoreDataFromModel <- function(object, unseenFeatures,unseenResponse){
  UseMethod('scoreDataFromModel',object)
}

#' \code{scoreDataFromModel} takes the new model and predicts a phenotype from an input set
#' @param model
#' @param unseenFeatures
#' @param unseenResponse
#' @keywords
#' @export
#' @return list of scores for each of the columns of the unseen feature data frame
scoreDataFromModel.fendR <- function(object, unseenFeatures,unseenResponse){
  #new.df<-left_join(unseenFeatures,object$featureData,by='Sample')
  model.preds<-sapply(names(object$featureModel),function(p){
    #subset for phenotype if we have different ones for each drug
    if('Phenotype'%in%names(unseenFeatures))
      ddf<-subset(unseenFeatures,Phenotype==p)[,-which(colnames(unseenFeatures)=='Phenotype')]
    else
      ddf<-unseenFeatures

    newvec<-dplyr::select(ddf,Gene,Value)%>%tidyr::spread(Gene,value=Value)
    predict(object$featureModel[[p]],newdata=newvec)[[1]]
  })

  model.preds<-data.frame(Prediction=model.preds,Actual=unseenResponse$Response[match(names(object$featureModel),unseenResponse$Phenotype)])

  return(model.preds)
}
