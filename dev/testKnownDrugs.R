##
# testKnownDrugs
#
# This script iterates through a set of high-throughput screens with
# matching expression data and compares the ability of this approach
# to identify known drugs that affect cell viability
##
library(fendR)

eset.file='syn11912257'
viper.file='syn11910413'
thresholds=c(0.25,0.75)

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

  #load eset **from SYNAPSE**
 # eset<-loadEset(synGet(rna.seq.data)$path,synGet(pheno.file)$path,useEntrez=TRUE)
 #

  eset<-readRDS(synGet(eset.file)$path)
  pset<-addResponseClass(eset,thresholds)

  #get drugs that have target ids
  matched.ids <- getDrugIds(varLabels(pset))
  tested.drugs <- matched.ids$ids

  if(!missing(drug.name)){
    inds <- which(tolower(matched.ids$drugs)%in%tolower(drug.name))
    if(length(inds)>0)
      matched.ids<-matched.ids[inds,]
  }
  library(viper)

  #get those with significantly differentially expressed genes
  all.pvals<-sapply(tolower(matched.ids$drugs),function(drug) rowTtest(pset, pheno=drug,group1='High',group2='Low')$p.value)

  sig.genes<-apply(all.pvals,2,function(x) length(which(p.adjust(x)<0.05)))

  nz.sig<-names(which(sig.genes>5))
  print(paste("found",length(nz.sig),'drugs at least 5 differentially expressed genes'))

  #build network
  drug.graph <- loadDrugGraph()
  combined.graph <-buildNetwork(drug.graph)
  all.drugs <- getDrugsFromGraph(drug.graph)

  v.obj <- readRDS(synGet(viper.file)$path)


  #TODO: make this multi-core, possibly break into smaller functions
  all.res <- lapply(nz.sig,function(drug,pset,all.drugs,w,b,mu){
    #create viper signature from high vs. low
     high = which(pData(pset)[[drug]] =='High')
     low = which(pData(pset)[[drug]]=='Low')
    v.res<-getViperForDrug(v.obj,high,low,0.1,TRUE)

    pcsf.res <-runPcsfWithParams(ppi=combined.graph,terminals=abs(v.res),dummies=all.drugs,w=w,b=b,mu=mu,doRand=TRUE)

   drug.res <- V(pcsf.res)$name[which(V(pcsf.res)$type=='Compound')]
    #now get average tamimoto distance between that drug and drug of interest
    #print(paste("Selected",length(drug.res),'drugs in the graph'))

    ##collect stats, store in synapse table
    list(network=pcsf.res,
      drugs=drug.res,
      w=w,
      b=b,
      mu=mu,
      viperProts=names(v.res),
      inputDrug=drug)

  },pset,all.drugs=tested.drugs,w=w,b=b,mu=mu)

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



trackNetworkStats<-function(pcsf.res.list,synTableId='syn12000477',esetFileId,viperFileId){
  require(synapser)

  pcsf.parent='syn12000478'
  this.script='https://github.com/Sage-Bionetworks/fendR/blob/master/dev/testKnownDrugs.R'
  #decouple pcsf.res.list into data frame
  fin<-lapply(pcsf.res.list,function(x){
    #first store network
    network=x[['network']]
    drug=x[['inputDrug']]
    w=x[['w']]
    b=x[['b']]
    mu=x[['mu']]
    fname=paste(paste(esetFileId,viperFileId,drug,w,b,mu,sep='_'),'.rds',sep='')
    saveRDS(network,file=fname)
    res=synStore(File(fname,parentId=pcsf.parent),used=c(esetFileId,viperFileId),executed=this.script)
     upl<-data.frame(`Input Drug`=drug,w=w,beta=b,mu=mu,`Viper Proteins`=paste(sort(x$viperProts),collapse=','),`Output Drugs`=paste(sort(x$drugs),collapse=','),`Original eSet`=esetFileId,`Original metaViper`=viperFileId,`mean TMD`=0,`PCSF Result`=res$properties$id,check.names=F)

     tres<-synStore(Table(synTableId,upl))
  })

  #store as synapse table

}


all.res<-findDrugsWithTargetsAndGenes(eset.file='syn11912257',
  viper.file='syn11910413',
  thresholds=c(0.25,0.75),
  drug.name=c('parthenolide','gefitinib','selumetinib'))
