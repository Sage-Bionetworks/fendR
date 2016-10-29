##create test matrix from MAF file - Gene by sample

library(synapseClient)
library(dplyr)
library(cgdsr)
library(data.table)

synapseLogin()

##this is just a UCSC list of gene names
all.genes<<-unique(fread('https://raw.githubusercontent.com/sgosline/RASPathwaySig/master/data/ucsc_kgXref_hg19_2015_10_29.csv')$geneSymbol)

##CCLE Broad data
mycancerstudy<-'cellline_ccle_broad'
mycgds = CGDS("http://www.cbioportal.org/public-portal/")
profile<-"cellline_ccle_broad_mutations" ##think about adding CNA data

##now go through the process of pulling from cbioportal
caseLists<-getCaseLists(mycgds,mycancerstudy)
#print('Got caselists')
mutSamps<-caseLists$case_list_id[grep("sequenced",caseLists[,1])]
#print(paste('Collecting CCLE mutation data for',tiss,'tissue'))
gene.groups=split(all.genes, ceiling(seq_along(all.genes)/500))
dat<-lapply(gene.groups,function(g) getProfileData(mycgds,g,profile,mutSamps))

ddat<-matrix()
for(i in which(sapply(dat,nrow)!=0)){
  ddat<-cbind(ddat,dat[[i]])
}
nans<-which(apply(ddat,2,function(x) all(is.nan(x)||is.na(x))))
# nas<-which(apply(ddat,2,function(x) all(is.na(x))))

ddat<-ddat[,-nans]
##now set to binary matrix
dfdat<-apply(ddat,1,function(x){
  sapply(unlist(x),function(y) !is.na(y) && y!='NaN')
})

bin.mat<-dfdat*1
fname='CCLE_binary_mutation_matrix_ucscGenesFromCBioPortal.tsv'
write.table(bin.mat,file=fname,row.names=T,col.names=T,sep='\t',quote=F)

##now we can re-rupload to new synapse project
fendRDatDir='syn7465504'
synStore(File(fname,parentId=fendRDatDir))

#might as well push to the CCLE repo
synStore(File(fname,parentId='syn5706496'))

###now get the CTRP data


