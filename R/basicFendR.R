
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
#' @import dplyr
#' @return list of gene features for each phenotype/drug response
createNewFeaturesFromNetwork.basicFendR<-function(object,testDrugs=NA){
    library(dplyr)
    ##figure out which phenotypes have both feature data and outcome data
    phenos<-intersect(object$sampleOutcomeData$Phenotype,object$phenoFeatureData$Phenotype)
    all.phenos<-union(object$sampleOutcomeData$Phenotype,object$phenoFeatureData$Phenotype)
    print(paste("Found",length(phenos),'phenotypes that have feature data and outcome data out of',length(all.phenos)))

    if(!is.na(testDrugs)&&length(intersect(testDrugs,phenos))>1){
      print(paste("Reducing scope to only focus on",paste(testDrugs,collapse=',')))
      phenos=intersect(phenos,testDrugs)
    }


    #TODO: investigate how these can be done with dplyr/mutate?

    ##for each phenotype, update the gene value by the shortest path to the gene target
    pheno.updates<-lapply(phenos,function(p){
      dt<-as.character(subset(object$phenoFeatureData,Phenotype==p)$Gene)
      print(paste('Calculating shortest path to',p,'target(s):',paste(dt,collapse=',')))
       #calculate shortest path between all drug targets and genes in feature set (that are in network)
        gd<-distances(object$network,intersect(dt,names(V(object$network))),
            intersect(object$featureData$Gene,names(V(object$network))))
        #get minimum across all drug targets
        min.to.targ<-apply(gd,2,min)
        #remove Inf values
        min.to.targ<-min.to.targ[which(is.finite(min.to.targ))]
        min.to.targ
    })

    ##update from featureData the score by shortest weighted path to target genes
    ##this is ridiculously time-consuming
    pheno.features<-lapply(pheno.updates,function(x){
      ##find out features with graph data
      nzFeatures<-intersect(names(x),object$featureData$Gene)

      zFeat<-setdiff(object$featureData$Gene,names(x))

      #create new data frame with features
      ddf<-data.frame(Gene=as.character(c(nzFeatures,zFeat)),
        FracDistance=c(1/x[nzFeatures],rep(0,length(zFeat))))
      ddf$FracDistance[!is.finite(ddf$FracDistance)]<-0

      new.fd<-left_join(object$featureData,ddf,by="Gene")%>%mutate(NetworkValue=Value+FracDistance)
      return(new.fd)
    })
    #move to data frame
    newdf<-do.call('rbind',pheno.features)
    #now add back phenotype information
    phen<-c()
    for(i in 1:length(phenos))
      phen<-c(phen,rep(phenos[i],nrow(pheno.features[[i]])))
    newdf$Phenotype<-phen


    ##Reduction strategy:
    #if we have multiple drugs: remove any genes that don't change across drugs.
    #eventually do something more complicated
    if(is.na(testDrugs)||length(testDrugs>1)){
      gene.var<-newdf%>%group_by(Gene)%>%summarize(Variance=var(NetworkValue))
      nzvars<-which(gene.var$Variance>0)
      genes<-gene.var$Gene[nzvars]
      print(paste('Keeping',length(nzvars),'gene values that change across drug treatments out of',length(gene.var$Gene)))
    }else{
       nzvars<-which(mutate(newdf,Diff=Value-NetworkValue)$Diff!=0)
       genes<-newdf$Gene[nzvars]
      print(paste('Keeping',length(nzvars),'gene values are altered by the network out of',nrow(newdf)))
    }
    newdf<-subset(newdf,Gene%in%genes)

    object$remappedFeatures<-newdf%>%select(Gene,Sample,Phenotype,
      Value=NetworkValue)

  return(object)

}


#' \code{scoreDataFromModel} takes the new model and predicts a phenotype from an input set
#' @param model
#' @param unseenData
#' @keywords
#' @export
#' @return list of scores for each of the columns of the unseen feature data frame
scoreDataFromModel.basicFendR<-function(object,unseenData){

}


