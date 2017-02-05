##create test matrix from MAF file - Gene by sample

library(synapseClient)
library(dplyr)
library(cgdsr)
library(data.table)

synapseLogin()
##########MUTATION PROFILE#################

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
colnames(bin.mat)[which(colnames(bin.mat)=='TT_THYROID')]<-'TTTHYROID'
colnames(bin.mat)<-sapply(colnames(bin.mat),function(x) unlist(strsplit(x,split='_'))[1])
fname='../inst/CCLE_binary_mutation_matrix_ucscGenesFromCBioPortal.tsv'
write.table(bin.mat,file=fname,row.names=T,col.names=T,sep='\t',quote=F)

##now we can re-rupload to new synapse project
this.script='https://raw.githubusercontent.com/Sage-Bionetworks/fendR/sara/testDataPrep/formatCcleData.R'
fendRDatDir='syn7465504'
synStore(File(fname,parentId=fendRDatDir),used=list(list(url=this.script)))

#might as well push to the CCLE repo
synStore(File(fname,parentId='syn5706496'),used=this.script)

########DRUG RESPONSE DATA #########################
###now get the CTRP data
auc_data<-synGet('syn5622711')@filePath
auc_dat=read.table(auc_data,sep='\t',header=T,quote='"')

ccl_data=synGet('syn5632194')@filePath
ccl_dat=read.table(ccl_data,header=T,sep='\t')
ccl_metadata<-synGet('syn5632192')@filePath
ccl_metadat<-read.table(ccl_metadata,sep='\t',header=T,as.is=T)

drug_data=synGet("syn5632193")@filePath
drug_dat=read.table(drug_data,sep='\t',header=T,fill=T,quote='"')


#now match it all up into a single data frame
new.df<-data.frame(AUC=as.numeric(auc_dat$area_under_curve),
                   Drug=drug_dat$cpd_name[match(auc_dat$master_cpd_id,drug_dat$master_cpd_id)],
                   CCL_id=ccl_dat$master_ccl_id[match(auc_dat$experiment_id,ccl_dat$experiment_id)])
new.df$CCL=ccl_metadat$ccl_name[match(new.df$CCL_id,ccl_metadat$master_ccl_id)]

require(reshape2)
drugmat<-acast(new.df,CCL~Drug,value.var="AUC",fun.aggregate=function(x) mean(x,na.rm=T))

fname='../inst/CTRP_v20_AUC_vales_by_drug.tsv'
write.table(drugmat,file=fname,sep='\t',row.names=T,col.names=T)

##upload to Synapse
synStore(File(fname,parentId=fendRDatDir),used=list(list(url=this.script)))

##lastly get drug target data
pheno.data<-read.table('../inst/CTRP_v20_AUC_vales_by_drug.tsv')
drug.target<-read.table(synGet('syn5632193')@filePath,header=T,sep='\t')
targs<-drug.target$gene_symbol_of_protein_target[match(colnames(pheno.data),sapply(drug.target$cpd_name,function(x) gsub(' |-','.',x)))]
names(targs)<-colnames(pheno.data)
targs<-targs[-which(is.na(targs))]

sing.targs<-sapply(as.character(targs),function(x) unlist(strsplit(x,split=';')))
#names(sing.targs)<-names(targs)

drugs.to.targs<-NULL
for(i in 1:length(sing.targs))
  drugs.to.targs<-rbind(drugs.to.targs,cbind(Drug=rep(names(targs)[i],length(sing.targs[[i]])),Target=sing.targs[[i]]))

    #data.frame(Drug=names(unlist(sing.targs)),Target=unlist(sing.targs))
write.table(drugs.to.targs,row.names=F,col.names=T,file='../inst/CTRP_v20_drug_target_vals.tsv',sep='\t')
synStore(File('../inst/CTRP_v20_drug_target_vals.tsv',parentId='syn7465504'),used=list(list(entity='syn5632193'),list(url=this.script)))
#drugs.to.targs=data.frame(Drug=colnames(pheno.data),Target=targs)
#drugs.to.targs<-drugs.to.targs[-which(is.na(drugs.to.targs$Target)),]
#remove nas

##copy to inst directory



