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

  require(parallel)
  require(Biobase)
#  cl <- makeCluster(nnodes=8)


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
  v.obj <- readRDS(synapser::synGet(viper.file)$path)

  matched.drugs <- which(sapply(toupper(varLabels(pset)),function(x) unlist(strsplit(x,split='_'))[1])%in%matched.ids$drugs)

  #get those with significantly differentially expressed genes
  all.vprots<-sapply(varLabels(pset)[matched.drugs],function(drug){
    high = which(Biobase::pData(pset)[[drug]] =='High')
    low = which(Biobase::pData(pset)[[drug]]=='Low')
   # print(paste("found",length(high),'high and',length(low),'low samples for',drug,sep=' '))
    res<-fendR::getViperForDrug(v.obj,high,low,0.01,TRUE,FALSE)
    print(paste("Found ",paste(names(res),collapse=','),' for drug ',drug))
    return(res)
  })
  names(all.vprots)<-tolower(varLabels(pset)[matched.drugs])

 # all.pvals<-sapply(tolower(matched.ids$drugs),function(drug) viper::rowTtest(pset, pheno=drug,group1='High',group2='Low')$p.value)
#  sig.genes<-apply(all.pvals,2,function(x) length(which(p.adjust(x)<0.05)))

  nz.sig<-which(sapply(all.vprots,length)>5)
  print(paste("found",length(nz.sig),'drugs at least 5 differentially expressed prots'))

  #build network
  drug.graph <- fendR::loadDrugGraph()
  combined.graph <-fendR::buildNetwork(drug.graph)
  all.drugs <- fendR::getDrugsFromGraph(drug.graph)

  fname=paste(paste(eset.file,viper.file,w,b,mu,paste(thresholds,collapse='_'),sep='_'),'.rds',sep='')
  #print(names(all.vprots)[nz.sig])
  #TODO: make this multi-core, possibly break into smaller functions
  all.res <- mclapply(names(all.vprots)[nz.sig],function(drug,all.vprots,tested.drugs,w,b,mu,fname){
    #create viper signature from high vs. low
    cat(drug)
  #print(high)
  v.res=all.vprots[[drug]]
  newf=paste(drug,fname,sep='_')

  if(file.exists(newf)){
    pcsf.res<-readRDS(newf)
  } else{
   # print(v.res)
    pcsf.res.id <-fendR::runPcsfWithParams(ppi=combined.graph,terminals=abs(v.res),dummies=all.drugs,w=w,b=b,mu=mu,doRand=TRUE)
    pcsf.res <-fendR::renameDrugIds(pcsf.res.id,all.drugs)


    saveRDS(pcsf.res,file=newf)

    }
  drug.res <- igraph::V(pcsf.res)$name[which(igraph::V(pcsf.res)$type=='Compound')]
  cat(paste("Selected",length(drug.res),'drugs in the graph'))

    ##collect stats, store in synapse table
    list(network=pcsf.res,
      drugs=drug.res,
      w=w,
      b=b,
      mu=mu,
      viperProts=names(v.res),
      inputDrug=unlist(strsplit(drug,split='_'))[1],
      file=newf)

  },all.drugs=tested.drugs,all.vprots,w=w,b=b,mu=mu,fname,mc.cores=28)#.parallel=TRUE,.paropts = list(.export=ls(.GlobalEnv)))

  names(all.res)<-names(all.vprots)[nz.sig]
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
trackNetworkStats<-function(pcsf.res.list,synTableId='syn12209124',esetFileId,viperFileId,thresholds){
  require(synapser)

  pcsf.parent='syn12209125'
  this.script='https://github.com/Sage-Bionetworks/fendR/blob/master/dev/testKnownDrugs_NTAP.R'
  #decouple pcsf.res.list into data frame

#  require(doMC)
#  cl <- makeCluster(nnodes=8)
  require(parallel)
 # registerDoMC(cores=28)


  fin<-mclapply(pcsf.res.list,function(x,thresholds){
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
  },thresholds,mc.cores=28)#.parallel=TRUE,.paropts = list(.export=ls(.GlobalEnv)))
#  stopCluster(cl)
  #store as synapse table

}

####ntap files
eset.file='syn11912279'
viper.file='syn11912228'
#thresholds=c(0.25,0.75)
for(w in c(2,3,4,5)){
 for(b in c(1,2,5,10)){
  for(mu in c(5e-05,5e-04,5e-03,5e-02)){
  thresholds=c(0.05,0.95)
  all.res<-findDrugsWithTargetsAndGenes(eset.file=eset.file,
                                        viper.file=viper.file,
                                        thresholds=thresholds,
                                        w=w,b=b,mu=mu)#,
  #                                    drug.name=c('parthenolide','gefitinib','selumetinib'))
  trackNetworkStats(all.res,esetFileId=eset.file,viperFileId=viper.file,thresholds=thresholds)



  thresholds=c(0.25,0.75)
  all.res<-findDrugsWithTargetsAndGenes(eset.file=eset.file,
                                        viper.file=viper.file,
                                        thresholds=thresholds,
                                        w=w,b=b,mu=mu)#,
  #                                    drug.name=c('parthenolide','gefitinib','selumetinib'))
  trackNetworkStats(all.res,esetFileId=eset.file,viperFileId=viper.file,thresholds=thresholds)

}}}
