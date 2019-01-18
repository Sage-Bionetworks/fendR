##compare NF1 protein activity to drug sensitivity and gene expression

library(fendR)
library(plyr)
library(viper)

library(synapser)
synLogin()


CCLE=list(eset.file='syn12549491',viper.file='syn12549589')
Sanger=list(eset.file='syn12549635',viper.file='syn12549806')

file.list<-list(updatedCCLE=list(rna='syn11902828',pheno=list(byDrug='syn7466611',byGene='syn7466552')),
  updatedSanger=list(rna='syn9987858',pheno=list(byDrug='syn9987866',byGene='syn9988097')))

##CCLE files
eset.file=CCLE$eset.file
viper.file=CCLE$viper.file
drug.data=synGet(file.list$updatedCCLE$pheno$byDrug)$path
mut.data=file.list$updatedCCLE$pheno$byGene

#rna.data<-loadSampleData(eset)
eset.file=CCLE$eset.file

eset<-readRDS(synapser::synGet(eset.file)$path)

v.obj <- readRDS(synapser::synGet(viper.file)$path)

#rna.data<-loadSampleData(rf)
pheno.data<-loadPhenotypeData(drug.data)
mut.data<-loadPhenotypeData(synGet(mut.data)$path)

nf1.rna=exprs(eset)['4763',]
nf1.prot=v.obj['4763',]
nf1.mut<-subset(mut.data,Sample=='NF1')$Response%>%setNames(subset(mut.data,Sample=='NF1')$Phenotype)

all.samps=intersect(names(nf1.rna),intersect(names(nf1.prot),names(nf1.mut)))

cor(nf1.rna[all.samps],nf1.prot[all.samps])
df=data.frame(mRNA=nf1.rna[all.samps],protein=nf1.prot[all.samps],NF1Status=factor(nf1.mut,labels =c("WT","Mutated"))[all.samps])

#correlation of NF1 status, mRNA levels and viper levels
ggplot(df)+geom_point(aes(x=mRNA,y=protein,col=NF1Status))

#TODO: add in tests and title
ggplot(df%>%gather(key=Molecule,value=levels,-NF1Status))+geom_boxplot(aes(x=NF1Status,y=levels,fill=Molecule))


##now do the same for Sanger
eset.file=Sanger$eset.file
viper.file=Sanger$viper.file
drug.data=synGet(file.list$updatedSanger$pheno$byDrug)$path
mut.data=file.list$updatedSanger$pheno$byGene



eset<-readRDS(synapser::synGet(eset.file)$path)

v.obj <- readRDS(synapser::synGet(viper.file)$path)

#rna.data<-loadSampleData(rf)
pheno.data<-loadPhenotypeData(drug.data)
mut.data<-loadPhenotypeData(synGet(mut.data)$path)

nf1.rna=exprs(eset)['4763',]
nf1.prot=v.obj['4763',]
nf1.mut<-subset(mut.data,Sample=='NF1')$Response%>%setNames(subset(mut.data,Sample=='NF1')$Phenotype)

all.samps=intersect(names(nf1.rna),intersect(names(nf1.prot),names(nf1.mut)))

cor(nf1.rna[all.samps],nf1.prot[all.samps])
nf.df=data.frame(mRNA=nf1.rna[all.samps],protein=nf1.prot[all.samps],NF1Status=factor(nf1.mut,labels =c("WT","Mutated"))[all.samps])
nf.df$Sample=rownames(df)


ggplot(nf.df)+geom_point(aes(x=mRNA,y=protein,col=NF1Status))

#TODO: add in tests and title
ggplot(nf.df%>%gather(key=Molecule,value=levels,-NF1Status))+geom_boxplot(aes(x=NF1Status,y=levels,fill=Molecule))

#conclusion: sanger dataset is more accurate reflection of NF1 protein activity

##now create the joined table
jdf<-left_join(pheno.data,nf.df,by="Sample")

pvals<-jdf%>%group_by(Phenotype)%>%dplyr::summarize(pval=wilcox.test(Response~NF1Status)$p.value)

ggplot(subset(new.jdf,Phenotype=='Imatinib'))+geom_point(aes(x=protein,y=Response,col=NF1Status))

computeMolCorrelations<-function(mol.level){
  #first just look at protein correlation
  cor.val=cor(t(mol.level))['4763',]
  cor.t<-apply(mol.level,1,function(x) cor.test(x,mol.level['4763',])$p.value)

  cor.t<-apply(mol.level,1,function(x) cor.test(x,mol.level['4763',])$p.value)

  sigs<-which(p.adjust(cor.t)<0.05)
  ##now what?
}

#identify those drugs whose activity correlates with transcript/protein level
computeDrugCorrelation<-function(jdf){

    #now compute drug correlation
    nf1.cors=jdf%>%group_by(Phenotype)%>%dplyr::summarize(NF1cor=cor(Response,protein,use='pairwise.complete.obs'))%>%arrange(desc(NF1cor))
}
