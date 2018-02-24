#' Store viper networks on synapse
# this could take a while with all the networks

library(synapser)
synLogin()

library(fendR)

storeViperOnSynapse <- function(rna.seq.data,pheno.file,parentId,datasetName){
  #load eset **from SYNAPSE*
  eset<-loadEset(synGet(rna.seq.data)$path,synGet(pheno.file)$path,useEntrez=TRUE)
  viper.res<-runViperOnDset(eset)
  fname=paste(datasetName,'withViper.rds',sep='')
  saveRDS(viper.res,file=fname)
  this.script=''
  synStore(File(fname,description='Viper run on gene expression data',parent=parentId),Activity('Viper on expression data',used=c(rna.seq.data,pheno.file),executed=this.script))

}
