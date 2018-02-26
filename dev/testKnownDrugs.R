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
#' @param rna.seq.data Tidied rna seq file
#' @param pheno.file Tidied drug response
#' @keywords
#' @export
#' @examples
#' @return
#'
findDrugsWithTargetsAndGenes <-function(eset.file,
    viper.file,
    thresholds=c(0.25,0.75)){

  library(synapser)
  synLogin()

  #load eset **from SYNAPSE*
 # eset<-loadEset(synGet(rna.seq.data)$path,synGet(pheno.file)$path,useEntrez=TRUE)
 #

  eset<-readRDS(synGet(eset.file)$path)
  pset<-addResponseClass(eset,thresholds)
  #get drugs that have target ids
  matched.ids <- getDrugIds(varLabels(pset))

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
  tested.drugs <- matched.ids$ids

  v.obj <- readRDS(synGet(viper.file)$path)

  #TODO: make this multi-core, possibly break into smaller functions
  all.res <- lapply(nz.sig,function(drug,pset,all.drugs){
    #create viper signature from high vs. low
     high = which(pData(pset)[[drug]] =='High')
     low = which(pData(pset)[[drug]]=='Low')
    v.res<-getViperForDrug(v.obj,high,low,0.1,TRUE)

    pcsf.res <-runPcsfWithParams(combined.graph,abs(v.res),all.drugs,w=5,b=5,mu=5e-02)
    drug.res <- intersect(all.drugs,V(pcsf.res)$name)
    #now get average tamimoto distance between that drug and drug of interest
    print(paste("Selected",length(drug.res),'drugs in the graph'))

    ##collect stats, store in synapse table

  })

  #TODO: evaluate all graphs with reference to network


}


trackNetworkStats<-function(pcsf.res,w,mu,b){


}
