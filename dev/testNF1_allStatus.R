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
   # print(paste("found",length(high),'high and',length(low),'low samples for',drug,sep=' '))
    res<-fendR::getViperForDrug(v.obj,wt,ko,0.001,TRUE,FALSE)
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
  dids<-getDrugIds(names(pData(eset)),split='_')[,1]
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
      ko=paste(conditions[[cond]]$KO,collapse=','),
      wt=paste(conditions[[cond]]$WT,collapse=','),
      viperProts=names(v.res),
    #  inputDrug=unlist(strsplit(drug,split='_'))[1],
      file=newf)

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
trackNetworkStats<-function(pcsf.res.list,synTableId='syn12334021',esetFileId,viperFileId){
  require(synapser)

  pcsf.parent='syn12333924'
  this.script='https://github.com/Sage-Bionetworks/fendR/blob/master/dev/testNF_Status.R'
  #decouple pcsf.res.list into data frame

#  require(doMC)
#  cl <- makeCluster(nnodes=8)
  require(parallel)
 # registerDoMC(cores=28)


  fin<-mclapply(pcsf.res.list,function(x){
    #first store network
    network=x[['network']]
    w=x[['w']]
    b=x[['b']]
    mu=x[['mu']]
    fname=x[['file']]
    ko=x[['ko']]
    wt=x[['wt']]

    res=synapser::synStore(File(fname,parentId=pcsf.parent),used=c(esetFileId,viperFileId),executed=this.script)
    upl<-data.frame(`NF1 KO`=ko,`NF1 WT`=wt,w=w,beta=b,mu=mu,
                     `Viper Proteins`=paste(sort(x$viperProts),collapse=','),
                     `Output Drugs`=paste(sort(x$drugs),collapse=','),
                     `Original eSet`=esetFileId,`Original metaViper`=viperFileId,
                     `mean TMD`=0,`PCSF Result`=res$properties$id,
                     `Mean Jaccard Distance`=0,
                     check.names=F)

     tres<-synapser::synStore(Table(synTableId,upl))
  },mc.cores=28)
  #.parallel=TRUE,.paropts = list(.export=ls(.GlobalEnv)))
#  stopCluster(cl)
  #store as synapse table

}

####ntap files
eset.file='syn12333863'
viper.file='syn12333867'

##ccle files
eset.file='syn12549491'

viper.file='syn12549589'
##sanger files

eset.file='syn12549635'
viper.file='syn12549806'

#for(w in c(2,3,4,5)){
# for(b in c(1,2,5,10)){
#  for(mu in c(5e-05,5e-04,5e-03,5e-02)){
   mu=0.005
   b=1
   w=2

  all.res<-findDrugsWithTargetsAndGenes(eset.file=eset.file,
                                        viper.file=viper.file,
                                        genotype='nf1',
                                        conditions=list(KOvsWT=list(KO=1,WT=0)),
                                        w=w,b=b,mu=mu)

  plotDrugs(eset,all.res$KOvsWT$drugs,'nf1')
  tab<-do.call(cbind,vertex.attributes(all.res$KOvsWT$network))
 # trackNetworkStats(all.res,esetFileId=eset.file,viperFileId=viper.file)




#}}}

