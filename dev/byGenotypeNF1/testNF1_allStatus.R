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
    genotype='nf1 genotype',
    conditions=list(homozygous=list(WT="+/+",KO="-/-"),
                    KOvsHets=list(WT=c("+/+","+/-"),KO="-/-"),
                    InclHets=list(WT="+/+",KO=c("+/-","-/-"))),
  w=2,
  b=1,
  mu=5e-04){

  library(synapser)
  synLogin()

  require(parallel)
  require(Biobase)
#  cl <- makeCluster(nnodes=8)

  eset<-readRDS(synGet(eset.file)$path)
  pset<-fendR::addGenotypeClass(eset,conditions,genotype)

  #get drugs that have target ids
#  matched.ids <- getDrugIds(varLabels(pset))
#  tested.drugs <- matched.ids$ids
  #print(matched.ids)

 # if(!missing(drug.name)){
#    inds <- which(tolower(matched.ids$drugs)%in%tolower(drug.name))
#    if(length(inds)>0)
#      matched.ids<-matched.ids[inds,]
#  }

  library(viper)
  v.obj <- readRDS(synapser::synGet(viper.file)$path)

 # matched.drugs <- which(sapply(toupper(varLabels(pset)),function(x) unlist(strsplit(x,split='_'))[1])%in%matched.ids$drugs)

  #get those with significantly differentially expressed genes
  all.vprots<-lapply(names(conditions),function(cond){
    wt = which(Biobase::pData(pset)[[cond]] =='WT')
    ko= which(Biobase::pData(pset)[[cond]]=='KO')

    wt.names=intersect(colnames(v.obj),Biobase::pData(pset)$Sample[wt])
    ko.names=intersect(colnames(v.obj),Biobase::pData(pset)$Sample[ko])
   # print(paste("found",length(high),'high and',length(low),'low samples for',drug,sep=' '))
    res<-fendR::getViperForDrug(v.obj,wt.names,ko.names,0.001,TRUE,FALSE)
    print(paste("Found ",paste(names(res),collapse=','),' for condition ',cond))
    return(res)
  })
  names(all.vprots)<-names(conditions)

 # all.pvals<-sapply(tolower(matched.ids$drugs),function(drug) viper::rowTtest(pset, pheno=drug,group1='High',group2='Low')$p.value)
#  sig.genes<-apply(all.pvals,2,function(x) length(which(p.adjust(x)<0.05)))

  nz.sig<-which(sapply(all.vprots,length)>5)
  print(paste("found",length(nz.sig),'drugs at least 5 differentially expressed prots'))

  #build network
  drug.graph <- fendR::loadDrugGraph()
  combined.graph <-fendR::buildNetwork(drug.graph)
  all.drugs <- fendR::getDrugsFromGraph(drug.graph)
  dids<-as.character(getDrugIds(names(pData(eset)),split='_')[,1])
  fname=paste(paste(eset.file,viper.file,w,b,mu,sep='_'),'.rds',sep='')
  #print(names(all.vprots)[nz.sig])
  #TODO: make this multi-core, possibly break into smaller functions
  all.res <- lapply(names(all.vprots)[nz.sig],function(cond,all.vprots,w,b,mu,fname,conditions){
    #create viper signature from high vs. low
    cat(cond)
  #print(high)
  v.res=all.vprots[[cond]]
  newf=paste(cond,fname,sep='_')

  if(file.exists(newf)){
    pcsf.res<-readRDS(newf)
  } else{
   # print(v.res)
    pcsf.res.id <-fendR::runPcsfWithParams(ppi=combined.graph,terminals=abs(v.res),dummies=dids,w=w,b=b,mu=mu,doRand=TRUE)
    pcsf.res <-fendR::renameDrugIds(pcsf.res.id,dids)
    saveRDS(pcsf.res,file=newf)

  }

  drug.res <- igraph::V(pcsf.res)$name[which(igraph::V(pcsf.res)$type=='Compound')]
  cat(paste("Selected",length(drug.res),'drugs in the graph'))

  pvalsAndFigs=plotDrugs(eset.file,drug.res,genotype)
  tab<-do.call(cbind,vertex.attributes(pcsf.res))
  ttab<-as.data.frame(tab)%>%filter(type=='Compound')%>%mutate(Drug=tolower(name))%>%dplyr::select(c(Drug,prize))
  #now we need to store all of these in the updated table.
  res=left_join(pvalsAndFigs,ttab,by='Drug')


    ##collect stats, store in synapse table
    list(network=pcsf.res,
      drugs=drug.res,
      w=w,
      b=b,
      mu=mu,
      ko=paste(conditions[[cond]]$KO,collapse=','),
      wt=paste(conditions[[cond]]$WT,collapse=','),
      viperProts=names(v.res),
    #  inputDrug=unlist(strsplit(drug,split='_'))[1],
      file=newf,
      compoundStats=res)

  },all.vprots,w=w,b=b,mu=mu,fname,conditions)#,mc.cores=28)#.parallel=TRUE,.paropts = list(.export=ls(.GlobalEnv)))

  names(all.res)<-names(all.vprots)[nz.sig]
  all.res

}


#'
#'trackNetworkStats takes a list of results from the drug test and shares them on synapse
#'@param pcsf.res.list
#'@param synTableId
#'@param esetFileId
#'@param viperFileId
#'
trackNetworkStats<-function(pcsf.res.list,synTableId='syn16780706',esetFileId,viperFileId,dsetName='',  pcsf.parent='syn15734434',  plot.parent='syn15734433'){
  require(synapser)


  this.script='https://github.com/Sage-Bionetworks/fendR/blob/master/dev/byGenotypeNF1/testNF1_allStatus.R'
  #decouple pcsf.res.list into data frame

#  require(doMC)
#  cl <- makeCluster(nnodes=8)
  require(parallel)
 # registerDoMC(cores=28)


  fin<-lapply(pcsf.res.list,function(x){
    #first store network
    network=x[['network']]
    w=x[['w']]
    b=x[['b']]
    mu=x[['mu']]
    fname=x[['file']]
    ko=x[['ko']]
    wt=x[['wt']]
    ds=x[['compoundStats']]%>%rename(Drug='Selected Drug',p.value='Drug Wilcoxon P-value')%>%mutate('Drug Prize Value'=as.numeric(prize))%>%ungroup()
    ds$`Drug Boxplot`=sapply(ds$figFile,function(y) synStore(File(y,parentId=plot.parent))$properties$id)
    res=synapser::synStore(File(fname,parentId=pcsf.parent),used=c(esetFileId,viperFileId),executed=this.script)

    ds=ds%>%dplyr::select(-figFile,-prize)
    #store image file
    upl<-data.frame(`NF1 KO`=ko,`NF1 WT`=wt,w=w,beta=b,mu=mu,
      `Viper Proteins`=paste(sort(x$viperProts),collapse=','),
      `Original eSet`=esetFileId,`Original metaViper`=viperFileId,
      `PCSF Result`=res$properties$id,`Dataset name`=dsetName,check.names=F)#,
    #                     check.names=F)

  upl2=merge(ds,upl)

     tres<-synapser::synStore(Table(synTableId,upl2))
  })#,mc.cores=28)
  #.parallel=TRUE,.paropts = list(.export=ls(.GlobalEnv)))
#  stopCluster(cl)
  #store as synapse table

}

####ntap files
#synIds<-list(NTAP=list(results='syn12333924',eset.file='syn12333863',viper.file='syn12333867',tableId='syn12334021'),
  synIds=list(
              CCLE=list(eset.file='syn12549491',viper.file='syn12549589'),
            Sanger=list(eset.file='syn12549635',viper.file='syn12549806'))

wvals=c(2,3,4,5)
bvals=c(1,2,5,10)
muvals=c(5e-05,5e-04,5e-03,5e-02)

#for(w in c(2,3,4,5)){
# for(b in c(1,2,5,10)){
#  for(mu in c()){

all.params=expand.grid(w=wvals,b=bvals,mu=muvals,dname=names(synIds))

#all.params=all.params[1:10,]

fr=mdply(.data=all.params,.fun=function(w,b,mu,dname){

  x=synIds[[dname]]

  all.res<-findDrugsWithTargetsAndGenes(eset.file=x$eset.file,
                                        viper.file=x$viper.file,
                                        genotype='nf1',
                                        conditions=list(KOvsWT=list(KO=1,WT=0)),
                                        w=w,b=b,mu=mu)


 trackNetworkStats(all.res,esetFileId=x$eset.file,viperFileId=x$viper.file, dsetName=dname)
},.parallel=TRUE)

