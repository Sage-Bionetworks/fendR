##
# testKnownDrugs
#
# This script iterates through a set of high-throughput screens with
# matching expression data and compares the ability of this approach
# to identify known drugs that affect cell viability
##
library(fendR)
library(plyr)


#' \code{findDrugsWithTargetsAndGenes} Identifies drugs in a
#' @param eset.file Expression set with expression and phenotype data
#' @param viper.file Viper file with networks for all phenotypes
#' @param drug.name
#' @param thresholds
#' @param w
#' @param b
#' @param mu
#' @keywords
#' @export
#' @examples
#' @return list of network result objects
#'
findDrugsWithTargetsAndGenes <-function(eset.file,
    viper.file,
  drug.name,
    thresholds=c(0.25,0.75),
  w=2,
  b=1,
  mu=5e-04){

  library(synapser)
  synLogin()

  require(doParallel)
  require(Biobase)
  cl <- makeCluster(12)

  registerDoParallel(cl)
    #require(parallel)
  #load eset **from SYNAPSE**
 # eset<-loadEset(synGet(rna.seq.data)$path,synGet(pheno.file)$path,useEntrez=TRUE)
 #

  eset<-readRDS(synGet(eset.file)$path)
  pset<-fendR::addResponseClass(eset,thresholds)

  #get drugs that have target ids
  matched.ids <- getDrugIds(varLabels(pset))
  tested.drugs <- matched.ids$ids
  #print(matched.ids)

  if(!missing(drug.name)){
    inds <- which(tolower(matched.ids$drugs)%in%tolower(drug.name))
    if(length(inds)>0)
      matched.ids<-matched.ids[inds,]
  }
  library(viper)

  #get those with significantly differentially expressed genes
  all.pvals<-sapply(tolower(matched.ids$drugs),function(drug) viper::rowTtest(pset, pheno=drug,group1='High',group2='Low')$p.value)

  sig.genes<-apply(all.pvals,2,function(x) length(which(p.adjust(x)<0.05)))

  nz.sig<-names(which(sig.genes>5))
  print(paste("found",length(nz.sig),'drugs at least 5 differentially expressed genes'))

  #build network
  drug.graph <- fendR::loadDrugGraph()
  combined.graph <-fendR::buildNetwork(drug.graph)
  all.drugs <- fendR::getDrugsFromGraph(drug.graph)

  v.obj <- readRDS(synGet(viper.file)$path)

  fname=paste(paste(eset.file,viper.file,w,b,mu,sep='_'),'.rds',sep='')

  #TODO: make this multi-core, possibly break into smaller functions
  all.res <- llply(nz.sig,.fun=function(drug,pset,all.drugs,w,b,mu,fname){
    #create viper signature from high vs. low
    print(drug)
     high = which(Biobase::pData(pset)[[drug]] =='High')
     low = which(Biobase::pData(pset)[[drug]]=='Low')
  #print(high)

  newf=paste(drug,fname,sep='_')
  v.res<-fendR::getViperForDrug(v.obj,high,low,0.1,TRUE)

  if(file.exists(newf)){
    pcsf.res<-readRDS(newf)
  } else{
   # print(v.res)
    pcsf.res <-fendR::runPcsfWithParams(ppi=combined.graph,terminals=abs(v.res),dummies=all.drugs,w=w,b=b,mu=mu,doRand=TRUE)
    saveRDS(pcsf.res,file=newf)

    }
  drug.res <- igraph::V(pcsf.res)$name[which(igraph::V(pcsf.res)$type=='Compound')]
  print(paste("Selected",length(drug.res),'drugs in the graph'))

    ##collect stats, store in synapse table
    list(network=pcsf.res,
      drugs=drug.res,
      w=w,
      b=b,
      mu=mu,
      viperProts=names(v.res),
      inputDrug=drug,
      file=newf)

  },pset,all.drugs=tested.drugs,w=w,b=b,mu=mu,fname,.parallel=TRUE,.paropts = list(.export=ls(.GlobalEnv)))

  stopCluster(cl)
  #TODO: evaluate all graphs with reference to network
  all.res

}


#'
#'plotGenesByDrug
#'Ranks cell lines by drug efficacy and then plots expression of
#'genes in gene list
#'@param eset
#'@param protMat
#'@param geneList
plotGenesByDrug<-function(eset,
    protMat,
    geneList,
    drug,
    drug.vals,
    genesOrProteins=c('genes','proteins')){

    require(pheatmap)
    require(viridis)
    rows=featureNames(eset)[which(p.adjust(all.pvals)<0.05)]
    cols=order(pData(eset)[,drug])
    cols=cols[which(!is.na(pData(eset)[cols,drug]))]
    fname=paste(drug,'response',length(rows),genesOrProteins,'heatmap.png',sep='')

    pheatmap(exprs(eset)[rows,cols],
    cluster_cols=F,
      color=viridis(20),
      show_rownames=F,show_colnames=F,clustering_distance_rows='correlation',
      annotation_col = data.frame(Response=drug.vals,AUC=pData(eset)[,drug]),
      main=paste('Response to',drug),filename=fname)

}



#'
#'trackNetworkStats takes a list of results from the drug test and shares them on synapse
#'@param pcsf.res.list
#'@param synTableId
#'@param esetFileId
#'@param viperFileId
#'@param thresholds
#'
trackNetworkStats<-function(pcsf.res.list,synTableId='syn12000477',esetFileId,viperFileId,thresholds){
  require(synapser)

  pcsf.parent='syn12000478'
  this.script='https://github.com/Sage-Bionetworks/fendR/blob/master/dev/testKnownDrugs.R'
  #decouple pcsf.res.list into data frame

  require(doParallel)
  cl <- makeCluster(12)

  registerDoParallel(cl)


  fin<-llply(pcsf.res.list,.fun=function(x,thresholds){
    #first store network
    network=x[['network']]
    drug=x[['inputDrug']]
    w=x[['w']]
    b=x[['b']]
    mu=x[['mu']]
    fname=x[['file']]

    res=synapser::synStore(File(fname,parentId=pcsf.parent),used=c(esetFileId,viperFileId),executed=this.script)
     upl<-data.frame(`Input Drug`=drug,w=w,beta=b,mu=mu,
                     `Viper Proteins`=paste(sort(x$viperProts),collapse=','),
                     `Output Drugs`=paste(sort(x$drugs),collapse=','),
                     `Original eSet`=esetFileId,`Original metaViper`=viperFileId,
                     `mean TMD`=0,`PCSF Result`=res$properties$id,
                     `Mean Jaccard Distance`=0,
                     Quantiles=paste(thresholds,collapse=','),
                     check.names=F)

     tres<-synapser::synStore(Table(synTableId,upl))
  },thresholds,.parallel=TRUE,.paropts = list(.export=ls(.GlobalEnv)))
  stopCluster(cl)
  #store as synapse table

}


all.res<-findDrugsWithTargetsAndGenes(eset.file='syn11912257',
  viper.file='syn11910413',
  thresholds=c(0.25,0.75),
  drug.name=c('parthenolide','gefitinib','selumetinib'))


all.res<-findDrugsWithTargetsAndGenes(eset.file=eset.file,
                                      viper.file=viper.file,
                                      thresholds=thresholds)#,
  #                                    drug.name=c('parthenolide','gefitinib','selumetinib'))

trackNetworkStats(all.res,esetFileId=eset.file,viperFileId=viper.file,thresholds=thresholds)

