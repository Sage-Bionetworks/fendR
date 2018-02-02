##
## Generate gene lists of differentially expressed genes between most and least
## responsive cell lines for a drug of interest
library(fendR)
rna.seq.data<-system.file('SANGER_brainarray_rma_expr.tsv.gz', package='fendR')
pheno.file<-system.file('SANGER_dose_response_AUC.tsv',package='fendR')


rna.seq.data<-system.file('CCLE_medianZscore_rnaSeq_ucscGenesFromCbioPortal.tsv.gz', package='fendR')
pheno.file<-system.file('CTRP_v20_AUC_vales_by_drug.tsv',package='fendR')



#' collect gene lists
#'
#' \code{loadSampleData} takes a file path and loads it into a data frame
#' @param Path to sample file
#' @keywords
#' @export
#' @examples
#' @return tidied data
collectGeneLists<-function(rna.seq.data,pheno.file,thresholds=c(0.25,0.75),parentId=''){

  #get data files
  rna.data<-loadSampleData(rna.seq.data)

  pheno.data<-loadPhenotypeData(pheno.file)

  samples<-intersect(rna.data$Sample,pheno.data$Sample)
  ##ascribe high/low values based on quantile data
  with.class<-pheno.data%>%group_by(Phenotype)%>%mutate(ResponseType=ifelse(Response<quantile(Response,thresholds,na.rm=T)[1],"Low",ifelse(Response>quantile(Response,thresholds,na.rm=T)[2],"High","Mid")))

  #create rna matrix
  rna.mat<-unique(rna.data)%>%spread(key=Sample,value=Value)
  rownames(rna.mat)<-rna.mat$Gene
  rna.mat<-select(rna.mat,-Gene)%>%select(samples)

  ##iterate over every drug

  #maybe we can use each row as input into viper
  phen.class<-spread(select(with.class,-Response),key=Phenotype,value=ResponseType)%>%as.data.frame()
  rownames(phen.class)<-phen.class$Sample
  phen.class<-phen.class%>%rename(sampleID=Sample)



#now i just need a regulon!!!
  eset<-Biobase::ExpressionSet(as.matrix(rna.mat))#,phenoData=Biobase::AnnotatedDataFrame(phen.class[which(rownames(phen.class)%in%samples),]))

  Biobase::phenoData(eset)<-Biobase::AnnotatedDataFrame(phen.class[which(rownames(phen.class)%in%samples),])

 # drug<-'selumetinib'
#  sig <- rowTtest(eset, pheno=drug,group1='High',group2='Low')$statistic

  all.sigs<-sapply(unique(with.class$Phenotype),function(drug) rowTtest(eset, pheno=drug,group1='High',group2='Low')$statistic)

  all.pvals<-sapply(unique(with.class$Phenotype),function(drug) rowTtest(eset, pheno=drug,group1='High',group2='Low')$p.value)

  corrected.pvals<-apply(all.pvals,2,p.adjust,method='BH')

  rownames(all.sigs)<-rownames(rna.mat)
  rownames(all.pvals)<-rownames(rna.mat)

  sig.genes<-apply(corrected.pvals,2,function(x) rownames(all.sigs)[which(x<0.1)])
  sapply(names(sig.genes)[which(sapply(sig.genes,length)>5)],function(x) write.csv(sig.genes[[x]],file=paste(gsub("[( ) /]","",x),'diffExGenesqLessThan0.1.txt',sep=''),row.names=F,col.names=F,quote=F))


  }

