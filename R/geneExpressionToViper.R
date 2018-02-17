#------------------
#
# convert differential gene expression dataset to viper
#
# Goal is to get a differential expression dataset into viper and then
#------------------

library(viper)
library(aracne.networks)
library(tidyverse)

#' \code{loadEset} takes the rna seq and pheno test files and builds expression set
#' @param rna.seq.data Tidied rna seq file
#' @param pheno.file Tidied drug response
#' @keywords
#' @export
#' @examples
#' @return expression set
loadEset<-function(rna.seq.data,pheno.file,useEntrez=TRUE){
  library(org.Hs.eg.db)

  #get data files
  rna.data<-loadSampleData(rna.seq.data)

  pheno.data<-loadPhenotypeData(pheno.file)

  samples<-intersect(rna.data$Sample,pheno.data$Sample)

  #create rna matrix
  rna.mat<-unique(rna.data)%>%spread(key=Sample,value=Value)
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
  phen.class<-spread(pheno.data,key=Phenotype,value=Response)%>%as.data.frame()
  rownames(phen.class)<-phen.class$Sample
  phen.class<-phen.class%>%dplyr::rename(sampleID="Sample")


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
  #extract pheno data to ascribe high/low
  pheno.data <-pData(eset)%>%gather(key="Phenotype",value="Response",-sampleID)
  ##ascribe high/low values based on quantile data
  with.class<-pheno.data%>%group_by(Phenotype)%>%mutate(ResponseType=ifelse(Response<quantile(Response,thresholds,na.rm=T)[1],"Low",ifelse(Response>quantile(Response,thresholds,na.rm=T)[2],"High","Mid")))

  phen.class<-spread(dplyr::select(data.frame(with.class),-Response),key=Phenotype,value=ResponseType)%>%as.data.frame()

  rownames(phen.class)<-phen.class$sampleID
  phen.class<-phen.class%>%dplyr::rename(Sample='sampleID')

  Biobase::phenoData(eset)<-Biobase::AnnotatedDataFrame(phen.class)

  return(eset)
}

#' \code{runViperOnDrugs} takes an expression set and a list of drugs and identifies a viper list of proteins that explain differential expression
#' @param eset is expression set
#' @param drugs is a list of drug names
#' @keywords
#' @export
#' @examples
#' @return viper object

runViperOnDrugs <- function(eset, drugs = c("BRD-K11533227","dacarbazine")){

  drugs<-c("BRD-K11533227","dacarbazine")##two drugs for which we have diff ex genes
  all.sigs<-sapply(drugs,
      function(drug) rowTtest(eset, pheno=drug,group1='High',group2='Low')$statistic)
  rownames(all.sigs)<-rownames(exprs(eset))


  #TODO: increase number of permuations! this is too low!!
  null.sigs<-lapply(drugs,function(drug) ttestNull(eset, pheno=drug,group1='High',group2='Low',per=100))

  #TODO: get aracne networks

  #TODO: run viper
  res <- msviper(all.sigs[,1], reg,null.sigs[[1]])

  return(res)
}


