##create test matrix from Sanger data

library(synapseClient)
#library(cgdsr)
#library(data.table)
require(plyr)
require(dplyr)


synapseLogin()
synId='syn7466648'
this.script='https://raw.githubusercontent.com/Sage-Bionetworks/fendR/master/dev/formatSangerData.R?token=ABwyOrOUhZANLj8CaNeuV5HGP-lhjdsFks5ZRZj7wA%3D%3D'
an<-'Formatting of SANGER dataset'
##this is just a UCSC list of gene names
#all.genes<<-unique(read.table('https://raw.githubusercontent.com/sgosline/RASPathwaySig/master/data/ucsc_kgXref_hg19_2015_10_29.csv')$geneSymbol)

require(readxl)
drug.targs<-readxl::read_excel('~/Downloads/Screened_Compounds.xlsx')%>%select(`Drug ID`,`Drug Name`,`Target`)
names(drug.targs)<-c("DRUG_ID","DRUG_NAME","Gene")
drug.targs<-dplyr::mutate(drug.targs,Phenotype=paste(DRUG_NAME,DRUG_ID,sep='_'))
all.targs<-tidyr::separate_rows(drug.targs,Gene,sep=', ')

write.table(all.targs,'../inst/SANGER_drugGeneTargets.tsv',sep='\t',quote=F)
synStore(File('../inst/SANGER_drugGeneTargets.tsv',parentId=synId),executed=list(list(url=this.script)),activityName=an)

##need to update column names and move ensg identifiers to HUGO
mrna.norm<-read.table(gzfile('~/Downloads/sanger1018_brainarray_ensemblgene_rma.txt.gz'),header=T,row.names=1)
colnames(mrna.norm)<-sapply(colnames(mrna.norm),function(x) gsub('X','',x))
library(biomaRt)
ensembl=useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host='www.ensembl.org')
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
entrez='ensembl_gene_id'
egene='hgnc_symbol'
gene.mapping<-getBM(attributes=c(egene,entrez),filters=c(entrez),values=rownames(mrna.norm),mart=ensembl)

mrna.norm$ensembl_gene_id<-rownames(mrna.norm)

full.ids<-left_join(gene.mapping,mrna.norm,key='ensembl_gene_id')
blanks=which(full.ids$hgnc_symbol=="")
full.ids<-full.ids[-blanks,]
rownames(full.ids)<-full.ids$hgnc_symbol
full.ids<-dplyr::select(full.ids,-hgnc_symbol,-ensembl_gene_id)
#write.table(round(full.ids,digits=4),file='../inst/SANGER_brainarray_rma_expr.tsv',quote=F,sep='\t')
gz<-gzfile('../inst/SANGER_brainarray_rma_expr.tsv.gz','w')
write.table(round(full.ids,digits=4),gz,quote=F,sep='\t')
close(gz)

synStore(File('../inst/SANGER_brainarray_rma_expr.tsv.gz',parentId=synId),executed=list(list(url=this.script)),activityName=an)

##need to figure out these files and reformat to  a single file
unzip('~/Downloads/CellLines_CG_BEMs.zip')
bin.dir='CellLines_CG_BEMs'

tab<-read.table(paste(bin.dir,'PANCAN_SEQ_BEM.txt',sep='/'),header=T)
  ftab<-tidyr::gather(tab,'Sample','Mutation',2:ncol(tab))
  ftab$Sample<-sapply(ftab$Sample,function(x) gsub("X","",x))
  colnames(ftab)[1]<-'Gene'

  bem_mat<-tidyr::spread(ftab,Sample,Mutation)
  rownames(bem_mat)<-bem_mat$Gene
  write.table(bem_mat%>%dplyr::select(-Gene),file='../inst/SANGER_binaryEventMatrix.tsv',sep='\t',quote=F)

  synStore(File('../inst/SANGER_binaryEventMatrix.tsv',parentId=synId),executed=list(list(url=this.script)),activityName=an)
##read in AUC values
drug.resp<-readxl::read_excel('~/Downloads/v17_fitted_dose_response.xlsx')%>%dplyr::select(COSMIC_ID,DRUG_ID,AUC)
full.drug<-dplyr::left_join(drug.resp,select(drug.targs,DRUG_ID,Phenotype),by='DRUG_ID')
drug.mat<-as.data.frame(tidyr::spread(full.drug%>%select(COSMIC_ID,AUC,Phenotype),Phenotype,AUC))
rownames(drug.mat)<-drug.mat$COSMIC_ID

write.table(select(drug.mat,-COSMIC_ID),file='../inst/SANGER_dose_response_AUC.tsv',sep='\t')


synStore(File('../inst/SANGER_dose_response_AUC.tsv',parentId=synId),executed=list(list(url=this.script)),activityName=an)

##lastly we want to have the CTRP/SANGER mapping file which is here
ccle.map<-readxl::read_excel('~/Downloads/GDSC-CCLE-CTRP_conversion.xlsx')%>%select(COSMIC_ID=`GDSC1000 cosmic id`,`CCLE name`)

cell.line.dis<-readxl::read_excel('~/Downloads/Cell_Lines_Details.xlsx')
min.info<-cell.line.dis%>%select(`Sample Name`,COSMIC_ID=`COSMIC identifier`,`Cancer Type\r\n(matching TCGA label)`)%>%left_join(ccle.map,by="COSMIC_ID")
write.table(min.info,file='../inst/SANGER_CCLE_map.tsv',sep='\t',quote=F,row.names=F)

synStore(File('../inst/SANGER_CCLE_map.tsv',parentId=synId),executed=list(list(url=this.script)),activityName=an)
