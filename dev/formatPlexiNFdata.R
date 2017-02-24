##format NTAP data into fendR format


library(dplyr)
require(githubr)

#here is where we will store the files
par.id<-'syn8282028'


#drug sensitivity
ncat.file="https://raw.githubusercontent.com/sgosline/pnfCellLines/master/bin/ncatsSingleAgentScreens.R"
source(ncat.file)
this.file='https://raw.githubusercontent.com/Sage-Bionetworks/fendR/master/dev/formatPlexiNFdata.R?token=ABwyOt9lkMgjvEtDEf408VHcaDvjCCUXks5Yt1v2wA%3D%3D'

targs<-ncatsDrugTargets()
colnames(targs)<-c("Phenotype","Gene")
write.table(targs,file='ncatsDrugTargetTidied.tsv',sep='\t',row.names=F,col.names=T)


aucs<-getValueForAllCells('TAUC')
aucs<-as.data.frame(aucs)
aucs$Phenotype<-rownames(aucs)
auc.vals<-tidyr::gather(aucs,Sample,Response,1:(ncol(aucs)-1))

lac50<-getValueForAllCells('LAC50')
lac50<-as.data.frame(lac50)
lac50$Phenotype<-rownames(lac50)
lac.vals<-tidyr::gather(lac50,Sample,Response,1:(ncol(lac50)-1))

maxr<-getValueForAllCells("MAXR")
maxr<-as.data.frame(maxr)
maxr$Phenotype<-rownames(maxr)
maxr.vals<-tidyr::gather(maxr,Sample,Response,1:(ncol(maxr)-1))

write.table(auc.vals,'ncatsNtapTaucTidied.tsv',sep='\t',row.names=F,col.names=T)
write.table(lac50,'ncatsNtapLac50Tidied.tsv',sep='\t',row.names=F,col.names=T)
write.table(maxr,'ncatsNtapMaxrTidied.tsv',sep='\t',row.names=F,col.names=T)
#now write all to directory and upload to synapse

#RNA-Seq
rna.seq.data<-read.table(synGet('syn7124098')@filePath,header=T,as.is=T)
rna.seq.data<-data.frame(rna.seq.data)
rna.seq.data$Gene<-rownames(rna.seq.data)
rna.seq<-tidyr::gather(rna.seq.data,Sample,Value,1:(ncol(rna.seq.data)-1))
rna.seq$Sample<-sapply(as.character(rna.seq$Sample),function(x) gsub("..mixed.clone."," (mixed clone)",gsub("..single.clone."," (single clone)",x,fixed=T),fixed=T))
write.table(rna.seq,'ntapRnaSeqTpmTidied.tsv',sep='\t',row.names=F,col.names=T)

#mutation data
mut.data<-read.table(synGet('syn6174638')@filePath,sep='\t',header=T,as.is=T)
red.mut.data<-subset(mut.data,KnownGeneExonFunction%in%c("nonsynonymous SNV","stoploss SNV",'stopgainSNV','frameshift deletion','frameshiftinsertion'))
md.mat<-reshape2::acast(red.mut.data,KnownGeneGeneName~CellLine)
md.mat[which(md.mat>0,arr.ind=T)]<-1
md.mat<-as.data.frame(md.mat)
md.mat$Gene<-rownames(md.mat)

md<-tidyr::gather(md.mat,Sample,Value,1:(ncol(md.mat)-1))

write.table(md,'ntapMutDataTidied.tsv',sep='\t',row.names=F,col.names=T)

file.list=c("ncatsDrugTargetTidied.tsv",'ncatsNtapTaucTidied.tsv','ncatsNtapLac50Tidied.tsv','ncatsNtapMaxrTidied.tsv','ntapRnaSeqTpmTidied.tsv','ntapMutDataTidied.tsv')
for (file in file.list)
  synStore(File(file,parentId=par.id),executed=list(list(url=ncat.file),list(url=this.file)))
