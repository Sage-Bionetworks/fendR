##compare cell line datasets

require(tidyr)

#' create a single data frame taht contains drug response values matched for
#' cell line and drug
mapCelLLinesAndDrugs<-function(){
  tab<-read.table(system.file('SANGER_CCLE_map.tsv',package='fendR'),fill=T,header=T,sep='\t',na.strings=c("NA","N/A"),check.names=T)
  red.tab<-tab[which(!is.na(tab$CCLE.name)),]
  red.tab$CCLE.name=sapply(as.character(red.tab$CCLE.name),function(x) unlist(strsplit(x,'_'))[1])
  red.tab$COSMIC_ID<-as.character(red.tab$COSMIC_ID)

  pheno.file<-system.file('SANGER_dose_response_AUC.tsv',package='fendR')
  sanger.data<-loadPhenotypeData(pheno.file)
  sanger.data$Phenotype<-tolower(sapply(sanger.data$Phenotype,function(x) unlist(strsplit(x,split='_'))[1]))


  pheno.file<-system.file('CTRP_v20_AUC_vales_by_drug.tsv',package='fendR')
  ctrp.data<-loadPhenotypeData(pheno.file)
  ctrp.data$Phenotype<-tolower(ctrp.data$Phenotype)

  common.drugs<-intersect(ctrp.data$Phenotype,sanger.data$Phenotype)

  #now start joining
  full.sanger<-left_join(sanger.data,red.tab,by=c("Sample" = "COSMIC_ID"))
  matched<-full.sanger[which(!is.na(full.sanger$CCLE.name)),]
  matched<-matched[which(matched$Phenotype%in%common.drugs),]%>%mutate(SamPhen=paste(CCLE.name,Phenotype,sep='_'))%>%rename(SANGER.AUC=Response)

  res<-subset(ctrp.data,Phenotype%in%common.drugs)%>%mutate(SamPhen=paste(Sample,Phenotype,sep='_'))%>%rename(CTRP.AUC=Response)

  full.res<-inner_join(matched,res,by="SamPhen")

  ggplot(full.res)+geom_point(aes(x=SANGER.AUC,y=CTRP.AUC,col=Phenotype.x))
  full.res%>%group_by(Phenotype.x)%>%summarize(Drug.Cor=cor(SANGER.AUC,CTRP.AUC,use='pairwise.complete.obs',method='spearman'))%>%View

  full.res%>%group_by(CCLE.name)%>%summarize(Cell.Cor=cor(SANGER.AUC,CTRP.AUC,use='pairwise.complete.obs',method='spearman'))%>%View

  full.res%>%group_by(CCLE.name)%>%summarize(Cell.Cor=cor(SANGER.AUC,CTRP.AUC,use='pairwise.complete.obs',method='spearman'),DataPoints=length(intersect(which(!is.na(SANGER.AUC)),which(!is.na(CTRP.AUC)))))%>%View

  tab<-full.res%>%group_by(CCLE.name)%>%summarize(Cell.Cor=cor(SANGER.AUC,CTRP.AUC,use='pairwise.complete.obs',method='spearman'),DataPoints=length(intersect(which(!is.na(SANGER.AUC)),which(!is.na(CTRP.AUC)))))

  ggplot(tab)+geom_point(aes(x=DataPoints,y=Cell.Cor))
}

#' compareDrugValues compares
#'
compareDrugValues<-function(){


}
