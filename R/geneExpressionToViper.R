#------------------
#
#process differential expression data to extract proteins of interest
#
# Current goal is to get a differential expression dataset into viper
#------------------

library(viper)
library(tidyverse)

#' \code{loadEset} takes the rna seq and pheno test files and builds expression set
#' @param rna.seq.data Tidied rna seq file
#' @param pheno.file Tidied drug response
#' @keywords
#' @export
#' @examples
#' @return expression set
loadEset<-function(rna.data,pheno.data,useEntrez=TRUE){
  library(org.Hs.eg.db)

  ##work with list here
  if(is.list(pheno.data))
    samps<-unique(unlist(lapply(pheno.data,function(x) x$Sample)))
  else
    samps<-pheno.data$Sample

  samples<-intersect(rna.data$Sample,samps)

  #create rna matrix
  rna.mat<-NULL

  rna.mat<-unique(rna.data)%>%tidyr::spread(key=Sample,value=Value)
  #if(is.null(rna.mat))
  #  rna.mat<-unique(rna.data)%>%dplyr::mutate(id=1:n())%>%spread(key=Sample,value=Value,fill=0)%>%dplyr::select(-id)
  rownames(rna.mat)<-rna.mat$Gene
  rna.mat<-rna.mat[,which(colnames(rna.mat)%in%samples)]

  #map to ENTREZ for viper
  if(useEntrez){
    x <- org.Hs.egSYMBOL2EG
    # Get the entrez gene identifiers that are mapped to a gene symbol
    mapped_genes <- AnnotationDbi::mappedkeys(x)
    # Convert to a list
    xx <- AnnotationDbi::as.list(x[mapped_genes])

    inds<-match(rownames(rna.mat),AnnotationDbi::mappedkeys(x))
    red.mat<-rna.mat[!is.na(inds),]
    rownames(red.mat)<-xx[which(!is.na(inds))]
    rna.mat<-red.mat

  }

  #maybe we can use each row as input into viper
  if(is.list(pheno.data)){
    phen.class<-Reduce(full_join,lapply(pheno.data,function(pd) spread(pd,key=Phenotype,value=Response,convert=TRUE)%>%as.data.frame()))
  }else{
    phen.class<-spread(pheno.data,key=Phenotype,value=Response,convert=TRUE)%>%as.data.frame()
  }
  rownames(phen.class)<-phen.class$Sample
  colnames(phen.class)<-tolower(colnames(phen.class))

  phen.class<-phen.class%>%dplyr::rename(sampleID="sample")


  #now i just need a regulon!!!
  eset<-Biobase::ExpressionSet(as.matrix(rna.mat))#,phenoData=Biobase::AnnotatedDataFrame(phen.class[which(rownames(phen.class)%in%samples),]))

  Biobase::phenoData(eset)<-Biobase::AnnotatedDataFrame(phen.class[which(rownames(phen.class)%in%samples),])
  return(eset)

}

#' \code{addResponseClass} takes an expression set and thresholds data to call differential expression based on drug
#' @param eset is expression set
#' @param thresholds are 2 values to use as quantiles
#' @keywords
#' @export
#' @examples
#' @return expression set
addResponseClass<-function(eset,drug,thresholds=c(0.25,0.75)){
  library(Biobase)
  #extract pheno data to ascribe high/low
  pheno.data <-unique(pData(eset)%>%gather(key="Phenotype",value="Response",-sampleID))
  ##ascribe high/low values based on quantile data
  with.class<-pheno.data%>%group_by(Phenotype)%>%mutate(ResponseType=ifelse(Response<quantile(Response,thresholds,na.rm=T)[1],"Low",ifelse(Response>quantile(Response,thresholds,na.rm=T)[2],"High","Mid")))

  phen.class<-spread(dplyr::select(data.frame(with.class),-Response),key=Phenotype,value=ResponseType)%>%as.data.frame()

  rownames(phen.class)<-phen.class$sampleID
  phen.class<-phen.class%>%dplyr::rename(Sample='sampleID')

  Biobase::phenoData(eset)<-Biobase::AnnotatedDataFrame(phen.class)

  return(eset)
}

#' \code{addGenotypeClass} takes a WT and KO set of features and creates a single class
#' @param eset is expression set
#' @param conditions is a list of conditions to ascribe to the response,
#' @param genotype is a particular genotype to apply the ocnditions to
#' @keywords
#' @export
#' @examples
#' @return expression set
addGenotypeClass<-function(eset,conditions=list(homozygous=list(WT="+/+",KO="-/-")),genotype='nf1 genotype'){
  library(Biobase)
  #extract pheno data to ascribe high/low
  pheno.data <-unique(pData(eset)%>%gather(key="Phenotype",value="Response",-sampleID))
  ##ascribes genotype data base3d on NF1 status, assuming NF1 is only gene.
  with.class<-lapply(conditions,function(x) pheno.data%>%mutate(ResponseType=ifelse(Response%in%x$WT,"WT",ifelse(Response%in%x$KO,"KO","NA"))))

  with.cond<-do.call('rbind',lapply(names(conditions),function(x) cbind(with.class[[x]],Condition=rep(x,nrow(with.class[[x]])))))
  ##we could do this if we end up having combinations of genotypes...
  #with.cond$ResponseType[which(!with.cond$Phenotype%in%genotype)]<-'NA'

  ##or just subset
  with.cond<-subset(with.cond,Phenotype%in%genotype)
  phen.class<-spread(dplyr::select(data.frame(with.cond),-Response),key=Condition,value=ResponseType)%>%as.data.frame()

  rownames(phen.class)<-phen.class$sampleID
  phen.class<-phen.class%>%dplyr::rename(Sample='sampleID')

  Biobase::phenoData(eset)<-Biobase::AnnotatedDataFrame(phen.class)

  return(eset)
}



#' \code{runViperOnDset} takes an expression set and runs viper with all aracne
#' # networks
#' @param eset is expression set
#' @keywords
#' @export
#' @examples
#' @return viper object
runViperOnDset <- function(eset){
  library(aracne.networks)
  #get aracne networks
  net.names <- data(package="aracne.networks")$results[, "Item"]
  all.networks <- lapply(net.names,function(x) get(x))
  names(all.networks) <- net.names

  res <- viper(Biobase::exprs(eset), all.networks)
  return(res)
}
#
#' \code{getViperForDrug} takes an expression set and a list of drugs and identifies a viper list of proteins that explain differential expression
#' @param drugs is a list of drug names
#' @keywords
#' @export
#' @examples
#' @return viper object
getViperForDrug <- function(v.res,high,low,pvalthresh=0.05,useEntrez=TRUE,p.corr=TRUE){

  #TODO: increase number of permuations! this is too low!!

   sig <-viper::rowTtest(v.res[,high],v.res[,low])$statistic
   pval<-viper::rowTtest(v.res[,high],v.res[,low])$p.value
   if(p.corr){
     pval <- p.adjust(pval)
   }
  sig.ps<-which(pval<pvalthresh)
  ret<-sig[sig.ps]
  names(ret)<-rownames(v.res)[sig.ps]

  if(useEntrez){
    library(org.Hs.eg.db)

      #we have to map back to HUGO
      x <- org.Hs.egSYMBOL2EG
      # Get the entrez gene identifiers that are mapped to a gene symbol
      mapped_genes <- AnnotationDbi::mappedkeys(x)
      # Convert to a list
      xx <- AnnotationDbi::as.list(x[mapped_genes])
      inds=match(names(ret),xx)
      names(ret)<-names(xx)[inds]


  }

  return(ret)
}


