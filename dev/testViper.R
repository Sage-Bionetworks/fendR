##
## Generate gene lists of differentially expressed genes between most and least
## responsive cell lines for a drug of interest
library(fendR)
library(viper)
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
eset<-loadEset(rna.seq.data,pheno.file,useEntrez=TRUE)
pset<-addResponseClass(eset)

#adjfile='~/Code/fendR/dev/tf2gene.adj'
#load regulon
#reg<-aracne2regulon(adjfile, eset, verbose = TRUE)
#load('dev/dnaseCCLEregulon.Rdata')
#now for each drug, identify signature, run viper.

drugs<-c("BRD-K11533227","dacarbazine")##two drugs for which we have diff ex genes
all.sigs<-sapply(drugs,function(drug) rowTtest(pset, pheno=drug,group1='High',group2='Low')$statistic)
rownames(all.sigs)<-rownames(exprs(pset))

drugs<-c("BRD-K11533227","dacarbazine")
null.sigs<-lapply(drugs,function(drug) ttestNull(pset, pheno=drug,group1='High',group2='Low',per=100))

library(aracne.networks)
net.names <- data(package="aracne.networks")$results[, "Item"]
all.networks <- lapply(net.names,function(x) get(x))
names(all.networks) <- net.names

#run viper
sig1<-all.sigs[,1]
names(sig1)<-rownames(all.sigs)
res1 <- viper(sig1, all.networks,method="none")#,null.sigs[[1]],method='ttest')
#
res2 <- msviper(all.sigs[,2], all.networks,method="none")# ,null.sigs[[2]])

#res1.all
#res2.all

