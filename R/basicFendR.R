
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
basicFendR<-function(networkFile, featureData, phenoFeatureData,sampleOutcomeData){
 me <-fendR(networkFile, featureData, phenoFeatureData,sampleOutcomeData)
 class(me) <- append(class(me),'basicFendR')
 return(me)
}

######################################################################
# Set methods implemented by basicFendR


#' Load Network data
#'
#' \code{loadNetwork} takes a file path and formats it as network
#' @param Path to network file
#' @keywords network feather
#' @export
#' @return a fendR object with the graph parameter populated with an iGraph object where edge weights represent distance between nodes (smaller means *more* association)
#' @examples
loadNetwork.basicFendR <- function(fObj){
  tab<-read.table(fObj$network,stringsAsFactors =FALSE)
  net<-igraph::graph_from_data_frame(tab,directed=F)
  E(net)$weight<-1-min(tab[,3],1)
  net<-delete_vertices(net,'UBC')
  fObj$graph <- net
  return (fObj)
}


#' Engineer Features from Network
#' \code{createNewFeaturesFromNetwork} takes the phenotype scores and updates the gene scores based on phenotype
#' @param object That contains a data frame and network
#' @keywords
#' @export
#' @return list of gene features for each phenotype/drug response
createNewFeaturesFromNetwork.basicFendR<-function(object,testDrugs=NA){
    doMC::registerDoMC()

    ##figure out which phenotypes have both feature data and outcome data
    phenos<-intersect(object$sampleOutcomeData$Phenotype,object$phenoFeatureData$Phenotype)
    all.phenos<-union(object$sampleOutcomeData$Phenotype,object$phenoFeatureData$Phenotype)
    print(paste("Found",length(phenos),'phenotypes that have feature data and outcome data out of',length(all.phenos)))


    if(!is.na(testDrugs)&&any(testDrugs%in%phenos)){
      print(paste("Reducing scope to only focus on",paste(testDrugs,collapse=',')))
      phenos=phenos[phenos %in% testDrugs]
    }

    ##for each phenotype, update the gene value by the shortest path to the gene target
    pheno.updates<-plyr::ldply(phenos,function(p){
      dt<-as.character(subset(object$phenoFeatureData,Phenotype==p)$Gene)
      print(paste('Calculating shortest path to',p,'target(s):',paste(dt,collapse=',')))
       #calculate shortest path between all drug targets and genes in feature set (that are in network)
        gd<-distances(object$graph,intersect(dt,names(V(object$graph))),
            intersect(object$featureData$Gene,names(V(object$graph))))

        #get minimum across all drug targets
        min.to.targ<-apply(gd,2,min)
        #remove Inf values
        min.to.targ<-min.to.targ[which(is.finite(min.to.targ))]
        min.to.targ
    },.parallel = TRUE)
    pheno.updates<-t(pheno.updates)
    colnames(pheno.updates)<-phenos
    #pheno.updates<-data.frame(pheno.updates)
    #rownames(pheno.updates)<0
    ##update from featureData the score by shortest weighted path to target genes
    ##this is ridiculously time-consuming
    pheno.features<-plyr::ldply(phenos,function(y){
      x<-pheno.updates[,y]
      ##find out features with graph data
      nzFeatures<-intersect(names(x),object$featureData$Gene)

      zFeat<-setdiff(object$featureData$Gene,names(x))

      #create new data frame with features
      ddf<-data.frame(Gene=as.character(c(nzFeatures,zFeat)),
        FracDistance=c(1/x[nzFeatures],rep(0,length(zFeat))), stringsAsFactors=FALSE)

      ddf$FracDistance[!is.finite(ddf$FracDistance)]<-0
      new.fd<-left_join(object$featureData,ddf,by="Gene")#%>%mutate(NetworkValue=Value+FracDistance)
      new.fd$NetworkValue<-apply(select(new.fd,Value,FracDistance),1,sum)

      new.fd$Phenotype<-rep(y,nrow(new.fd))
      return(new.fd)

    },.parallel =FALSE)
    newdf<-pheno.features

    ##Reduction strategy:
    #if we have multiple drugs: remove any genes that don't change across drugs.
    #eventually do something more complicated
    if(is.na(testDrugs)||length(testDrugs>1)){
      gene.var<-newdf%>%dplyr::group_by(Gene)%>%dplyr::summarize(Variance=var(NetworkValue))
      var.thresh=0.99
      nzvars<-which(gene.var$Variance>quantile(gene.var$Variance,na.rm=T,var.thresh)[[paste(var.thresh*100,'%',sep='')]])

      #nzvars<-which(gene.var$Variance>0)
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

#' Get Engineered Features as model matrix
#' @description Gets a \code{list} of response matrices for a phenotype
#' @export
engineeredResponseMatrix.basicFendR<-function(fObj,phenotype=c()){
  if(length(phenotype)==0)
    phenotype <- unique(fObj$remappedFeatures$Phenotype)


  out.dat<-subset(fObj$sampleOutcomeData,Phenotype%in%phenotype)
  in.dat<-subset(fObj$remappedFeatures,Phenotype%in%phenotype)

  mod.df<-dplyr::inner_join(out.dat,dplyr::select(in.dat,Sample,Gene,Value),by="Sample")%>%dplyr::select(Sample,Gene,Value,Phenotype,Response)

#  mod.df<-filter(mod.df,!Sample%in%sampsToOmit)
  mod.df<-dplyr::mutate(mod.df,SamplePheno=paste(Sample,Phenotype,sep='_'))
  dupes<-which(duplicated(select(mod.df,Gene,SamplePheno)))
  res<-tidyr::spread(select(mod.df[-dupes,],Gene,Value,SamplePheno,Response,Phenotype,Sample),Gene,Value)

 # rownames(res)<-res$SamplePheno
  res<-res[,-which(colnames(res)%in%c('Gene','SamplePheno'))]

  avar<-apply(res,2,var)
  #var.thresh=0.99
 # zvar<-which(avar<quantile(avar,na.rm=T,var.thresh)[[paste(var.thresh*100,'%',sep='')]])#
  zvar<-which(avar==0)

  if(length(zvar)>0){
  print(paste('Removing',length(zvar),'features from matrix out of',length(avar)))
  res<-res[,-zvar]
  }
  res



}
