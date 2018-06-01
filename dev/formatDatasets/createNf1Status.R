#create NF1 file for CCLE data

require(synapser)
require(tidyverse)
synLogin()

mut.file<-read.table(synGet('syn7466552')$path,header=T)

tab<-gather(mut.file['NF1',],key=Sample)
tab$Phenotype=rep("NF1 Genotype",nrow(tab))
tab$Value=rep("+/+",nrow(tab))
tab$Value[which(tab$value==1)]<-'+/-'

head(tab)
write.table(tab%>%select(-value),file='NF1GenotypeCCLEcellLines.tsv',sep='\t',row.names=F,col.names=T)

synStore(File('NF1GenotypeCCLEcellLines.tsv',parentId='syn7465504'),used='syn7466552',executed=this.script)
