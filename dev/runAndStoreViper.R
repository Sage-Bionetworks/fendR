#' Store viper networks on synapse
# this could take a while with all the networks

library(synapser)
synLogin()

library(fendR)

storeViperOnSynapse <- function(rna.seq.data,pheno.file,used=c(),parentId,esetParent,datasetName){
  #load eset **from SYNAPSE*
  eset<-loadEset(rna.seq.data,pheno.file,useEntrez=TRUE)
  this.script='https://github.com/Sage-Bionetworks/fendR/pull/49/commits/d5aabd9260e0a3d15c26122c06ff963a798a87f6'

  #store the eset
  fname=paste(datasetName,'expressionSet.rds',sep='')
  saveRDS(eset,file=fname)
  synStore(File(fname,description='Expression set with drug info',parent=esetParent),activity=Activity('Expression data',used=used,executed=this.script))

 # viper.res<-runViperOnDset(eset)
  fname=paste(datasetName,'withViper.rds',sep='')
#  saveRDS(viper.res,file=fname)
#  synStore(File(fname,description='Viper run on gene expression data',parent=parentId),activity=Activity('Viper on expression data',used=used,executed=this.script))


}

file.list<-list(CCLE=list(rna='syn11902828',pheno='syn7466611'),Sanger=list(rna='syn9987858',pheno='syn9987866'),pNFCell=list(rna='syn8304627',pheno='syn8304620'))

sapply(names(file.list),function(x){
  rf<-synGet(file.list[[x]]$rna)$path
  pf<-synGet(file.list[[x]]$pheno)$path
  if(x%in%c("Sanger","CCLE")){
    #get data files
    rna.data<-loadSampleData(rf)
    pheno.data<-loadPhenotypeData(pf)
  }else{
    rna.data<-read.table(rf,header=T)
    pheno.data<-read.table(pf,header=T)
  }
  storeViperOnSynapse(rna.data,pheno.data,used=c(file.list[[x]]$rna,file.list[[x]]$pheno),'syn11889550','syn11912239',x)
})
