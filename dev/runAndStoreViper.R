#' Store viper networks on synapse
# this could take a while with all the networks

library(synapser)
synLogin()

library(fendR)

storeViperOnSynapse <- function(rna.seq.data,pheno.file,used=c(),parentId,esetParent,datasetName,useEntrez=useEntrez,runViper=TRUE){
  #load eset **from SYNAPSE*
  eset<-loadEset(rna.seq.data,pheno.file,useEntrez=useEntrez)
  this.script='https://raw.githubusercontent.com/Sage-Bionetworks/fendR/master/dev/runAndStoreViper.R?token=ABwyOvseVJgofvKwfoMxe2Zx0zYUQaTXks5bLRZzwA%3D%3D'

  #store the eset
  fname=paste(datasetName,'expressionSet.rds',sep='')
  saveRDS(eset,file=fname)
  synStore(File(fname,description='Expression set with drug info',parent=esetParent),activity=Activity('Expression data',used=used,executed=this.script))

  if(runViper){
    viper.res<-runViperOnDset(eset)
    fname=paste(datasetName,'withViper.rds',sep='')
    saveRDS(viper.res,file=fname)
    synStore(File(fname,description='Viper run on gene expression data',parent=parentId),activity=Activity('Viper on expression data',used=used,executed=this.script))
  }

}

#file.list<-list(CCLE=list(rna='syn11902828',pheno='syn7466611'),Sanger=list(rna='syn9987858',pheno='syn9987866'),pNFCell=list(rna='syn8304627',pheno='syn8304620'))

#now adding genotype updates to all cell lines!
file.list<-list(updatedPNF=list(rna='syn12333637',pheno=c('syn12333638',"syn8304620")),
  updatedCCLE=list(rna='syn11902828',pheno=c('syn7466611','syn7466552')),
  updatedSanger=list(rna='syn9987858',pheno=c('syn9987866','syn9988097')))

sapply(names(file.list),function(x){
  print(x)
  rf<-synGet(file.list[[x]]$rna)$path
  pf<-sapply(file.list[[x]]$pheno,function(y) synGet(y)$path)
  if(x%in%c("Sanger","CCLE",'updatedCCLE','updatedSanger')){
    #get data files
    rna.data<-loadSampleData(rf)
    pheno.data<-loadPhenotypeData(pf)
  }else{
    rna.data<-read.table(rf,header=T)
    pheno.data<-lapply(pf,function(y) read.table(y,header=T))
  }
 # if('Value'%in%colnames(pheno.data))
#    colnames(pheno.data)[which(colnames(pheno.data)=='Value')]<-'Response'

  storeViperOnSynapse(rna.data,pheno.data,used=c(file.list[[x]]$rna,file.list[[x]]$pheno),'syn11889550','syn12550830',x,useEntrez=FALSE,runViper=FALSE)
})
