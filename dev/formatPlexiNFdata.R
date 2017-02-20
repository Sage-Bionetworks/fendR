##format NTAP data into fendR format


library(dplyr)

#here is where we will store the files
par.id<-'syn8282028'


#drug sensitivity
source("https://raw.githubusercontent.com/sgosline/pnfCellLines/master/bin/ncatsSingleAgentScreens.R")

targs<-ncatsDrugTargets()
aucs<-getValueForAllCells('TAUC')
lac50<-getValueForAllCells('LAC50')
maxr<-getValueForAllCells("MAXR")

#RNA-Seq
rna.seq.data<-read.table(synGet('syn7124098')@filePath,header=T,as.is=T)

#mutation data
mut.data<-read.table(synGet('syn6174638')@filePath,sep='\t',header=T,as.is=T)
