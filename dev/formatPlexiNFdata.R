##format NTAP data into fendR format


library(dplyr)

#here is where we will store the files
par.id<-'syn8282028'


#drug sensitivity
ncat.file="https://raw.githubusercontent.com/sgosline/pnfCellLines/master/bin/ncatsSingleAgentScreens.R"
source(ncat.file)
this.file=''

targs<-ncatsDrugTargets()
colnames(targs)<-c("Phenotype","Gene")
write.table(targs,file='ncatsDrugTargetTidied.tsv',sep='\t',row.names=F,col.names=T)
synStore(File("ncatsDrugTargetTidied.tsv",parentId=par.id),executed=list(list(url=ncat.file),list(url=this.file)))

aucs<-getValueForAllCells('TAUC')
lac50<-getValueForAllCells('LAC50')
maxr<-getValueForAllCells("MAXR")

#RNA-Seq
rna.seq.data<-read.table(synGet('syn7124098')@filePath,header=T,as.is=T)

#mutation data
mut.data<-read.table(synGet('syn6174638')@filePath,sep='\t',header=T,as.is=T)
