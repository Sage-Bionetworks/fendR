##
## Generate gene lists of differentially expressed genes between most and least
## responsive cell lines for a drug of interest
library(fendR)
library(viper)
rna.seq.data<-system.file('SANGER_brainarray_rma_expr.tsv.gz', package='fendR')
pheno.file<-system.file('SANGER_dose_response_AUC.tsv',package='fendR')


rna.seq.data<-system.file('CCLE_medianZscore_rnaSeq_ucscGenesFromCbioPortal.tsv.gz', package='fendR')
pheno.file<-system.file('CTRP_v20_AUC_vales_by_drug.tsv',package='fendR')



loadEset<-function(rna.seq.data,pheno.file){
  #get data files
  rna.data<-loadSampleData(rna.seq.data)

  pheno.data<-loadPhenotypeData(pheno.file)

  samples<-intersect(rna.data$Sample,pheno.data$Sample)

  #create rna matrix
  rna.mat<-unique(rna.data)%>%spread(key=Sample,value=Value)
  rownames(rna.mat)<-rna.mat$Gene
  rna.mat<-select(rna.mat,-Gene)%>%select(samples)

  #maybe we can use each row as input into viper
  phen.class<-spread(pheno.data,key=Phenotype,value=Response)%>%as.data.frame()
  rownames(phen.class)<-phen.class$Sample
  phen.class<-phen.class%>%rename(sampleID=Sample)


  #now i just need a regulon!!!
  eset<-Biobase::ExpressionSet(as.matrix(rna.mat))#,phenoData=Biobase::AnnotatedDataFrame(phen.class[which(rownames(phen.class)%in%samples),]))

  Biobase::phenoData(eset)<-Biobase::AnnotatedDataFrame(phen.class[which(rownames(phen.class)%in%samples),])
  return(eset)

}

addResponseClass<-function(eset,thresholds=c(0.25,0.75)){
  #extract pheno data to ascribe high/low

  pheno.data <-pData(eset)%>%gather(key="Phenotype",value="Response",-sampleID)
  ##ascribe high/low values based on quantile data
  with.class<-pheno.data%>%group_by(Phenotype)%>%mutate(ResponseType=ifelse(Response<quantile(Response,thresholds,na.rm=T)[1],"Low",ifelse(Response>quantile(Response,thresholds,na.rm=T)[2],"High","Mid")))

  phen.class<-spread(select(data.frame(with.class),-Response),key=Phenotype,value=ResponseType)%>%as.data.frame()

  rownames(phen.class)<-phen.class$sampleID
  phen.class<-phen.class%>%rename(Sample='sampleID')

  Biobase::phenoData(eset)<-Biobase::AnnotatedDataFrame(phen.class)

  return(eset)
}

#' collect gene lists
#'
#' \code{loadSampleData} takes a file path and loads it into a data frame
#' @param Path to sample file
#' @keywords
#' @export
#' @examples
#' @return tidied data
collectGeneLists<-function(eset,parentId=''){

  pheno.data <-pData(eset)%>%gather(key="Phenotype",value="Response",-sampleID)

  all.sigs<-sapply(unique(pheno.data$Phenotype),function(drug) rowTtest(eset, pheno=drug,group1='High',group2='Low')$statistic)

  all.pvals<-sapply(unique(pheno.data$Phenotype),function(drug) rowTtest(eset, pheno=drug,group1='High',group2='Low')$p.value)

  corrected.pvals<-apply(all.pvals,2,p.adjust,method='BH')

  rownames(all.sigs)<-rownames(rna.mat)
  rownames(all.pvals)<-rownames(rna.mat)

  sig.genes<-apply(corrected.pvals,2,function(x) rownames(all.sigs)[which(x<0.1)])
  res<-sapply(names(sig.genes)[which(sapply(sig.genes,length)>5)],function(x) write.csv(sig.genes[[x]],file=paste(gsub("[( ) /]","",x),'diffExGenesqLessThan0.1.txt',sep=''),row.names=F,col.names=F,quote=F))


  }

#load eset
eset<-loadEset(rna.seq.data,pheno.file)
pset<-addResponseClass(eset)

adjfile='~/Code/fendR/dev/tf2gene.adj'
#load regulon
reg<-aracne2regulon(adjfile, eset, verbose = TRUE)

#now for each drug, identify signature, run viper.
drugs<-c("BRD-K11533227","dacarbazine")##two drugs for which we have diff ex genes
all.sigs<-sapply(drugs,function(drug) rowTtest(pset, pheno=drug,group1='High',group2='Low')$statistic)
rownames(all.sigs)<-rownames(exprs(pset))

drugs<-c("BRD-K11533227","dacarbazine")
null.sigs<-lapply(drugs,function(drug) ttestNull(pset, pheno=drug,group1='High',group2='Low',per=100))


#run viper
res1 <- msviper(all.sigs[,1], reg,null.sigs[[1]])
res2 <- msviper(all.sigs[,2], reg,null.sigs[[2]])

