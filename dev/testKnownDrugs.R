##
# testKnownDrugs
#
# This script iterates through a set of high-throughput screens with
# matching expression data and compares the ability of this approach
# to identify known drugs that affect cell viability
##
library(fendR)

rna.seq.data<-system.file('SANGER_brainarray_rma_expr.tsv.gz', package='fendR')
pheno.file<-system.file('SANGER_dose_response_AUC.tsv',package='fendR')


rna.seq.data<-system.file('CCLE_medianZscore_rnaSeq_ucscGenesFromCbioPortal.tsv.gz', package='fendR')
pheno.file<-system.file('CTRP_v20_AUC_vales_by_drug.tsv',package='fendR')


#' \code{findDrugsWithTargetsAndGenes} Identifies drugs in a
#' @param rna.seq.data Tidied rna seq file
#' @param pheno.file Tidied drug response
#' @keywords
#' @export
#' @examples
#' @return
#'
findDrugsWithTargetsAndGenes <-function(rna.seq.data,
    pheno.file,
    thresholds=c(0.25,0.75)){

  #load eset
  eset<-loadEset(rna.seq.data,pheno.file,useEntrez=TRUE)
  pset<-addResponseClass(eset,thresholds)

  #get drugs that have target ids
  matched.ids <- getDrugIds(varLabels(pset))

  #get those with significantly differentially expressed genes
  all.pvals<-sapply(tolower(matched.ids$drugs),function(drug) rowTtest(pset, pheno=drug,group1='High',group2='Low')$p.value)

  sig.genes<-apply(all.pvals,2,function(x) length(which(p.adjust(x)<0.05)))

  nz.sig<-names(which(sig.genes>5))
  print(paste("found",length(nz.sig),'drugs at least 5 differentially expressed genes'))

  #build network
  drug.graph <- loadDrugGraph()
  combined.graph <-buildNetwork()
  all.drugs <- getDrugsFromGraph(drug.graph)

  #TODO: make this multi-core
  all.res <- lapply(nz.sig,function(drug,pset,all.drugs){
    v.res <- runViperOnDrug(pset,drug)


  })

}
