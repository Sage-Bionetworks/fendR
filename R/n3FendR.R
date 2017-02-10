
## This is an implementation of the nearest-neighbor network (n3) fendR class.
## It is based on the methods used by Yuanfang Guan in DREAM challenges.

######################################################################
# Create the n3FendR class
#
# This is used to represent the nearest-neighbor network implementation of the fendR framework.

#' An S3 class to represent the nearest-neighbor network implementation of the fendR predictive
#' network algorithm.
#' @param network the file name of a feather or big.matrix object representing a dense matrix
#' @param featureData a data.frame that contains rows representing genes and columns representing samples
#' @param sampleOutcomeData a data.frame representing at least one column of phenotype and rows representing samples
#' @param phenoFeatureData a data.frame where rows represent genes and columns represent a relationship between phenotype and gene
#' @param target.genes a vector of target gene names (a subset of those in featureData) that will be used to sparsify the network
#' @param network.type the on-disk representation of the network (either "feather" or "big.memory")
#' @inheritParams fendR
#' @export
#' @return n3FendR object
n3FendR<-function(network, featureData, phenoFeatureData,sampleOutcomeData, target.genes, network.type){
 me <-fendR(network, featureData, phenoFeatureData,sampleOutcomeData)
 me <- append(me, list(target.genes = target.genes, network.type = network.type))
 class(me) <- append(class(me),'n3FendR')
 return(me)
}

######################################################################
# Set methods implemented by n3FendR



#' Engineer Features from Network
#' \code{createNewFeaturesFromNetwork} takes the phenotype scores and updates the gene scores based on phenotype
#' @param object That contains a data frame and network
#' @keywords
#' @export
#' @import dplyr
#' @return list of gene features for each phenotype/drug response
createNewFeaturesFromNetwork.n3FendR <- function(object, testDrugs = NA){
    library(dplyr)
    ##figure out which phenotypes have both feature data and outcome data
    phenos <- intersect(object$sampleOutcomeData$Phenotype,object$phenoFeatureData$Phenotype)
    all.phenos <- union(object$sampleOutcomeData$Phenotype,object$phenoFeatureData$Phenotype)
    cat(paste0("Found ", length(phenos), "phenotypes that have feature data and outcome data out of", length(all.phenos), "\n"))

    if(!is.na(testDrugs) && any(testDrug%in%phenos)) {
      cat(paste0("Reducing scope to only focus on ", paste(testDrugs, collapse=','), " drugs\n"))
      phenos <- phenos[phenos %in% testDrug]
    }

    ## Sparsify the data (featureData and network) by only considering a subset of
    ## curated target.genes -- these are likely drug targets and/or "cancer genes," etc.
    ## NB: by doing this up front/here we do not propagate mutations in non target genes
    ## to the target genes.
    full.gene.set <- unique(object$featureData$gene)
    num.genes.in.full.feature.space <- length(full.gene.set)
    reduced.gene.set <- full.gene.set
    if(!is.na(object$target.genes) && !is.null(object$target.genes)) {
      if(!any(object$target.genes %in% full.gene.set)) {
        stop("None of target genes are in the featureData")
      }
      reduced.gene.set <- full.gene.set[full.gene.set %in% object$target.genes]
      cat(paste0("Reducing feature space from ", num.genes.in.full.feature.space, " to ", length(reduced.gene.set), "\n"))
    } else {
      cat(paste0("Using all ", num.genes.in.full.feature.space, " genes in featureData.\n"))
      cat("Not sparsifying data or network\n")
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
    #newdf$Phenotype<-unlist(sapply(phenos,rep,nrow(object$featureData)))

#  pf<-gather(data.frame(pheno.features),"Phenotype","NetworkDistance",1:ncol(pheno.features))

    ##Reduction strategy:
    #if we have multiple drugs: remove any genes that don't change across drugs.
    #eventually do something more complicated
    if(is.na(testDrug)||length(testDrug>1)){
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
