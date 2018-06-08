##check to see the output of the CCLE/Sanger datasets

ccle.tab<-'syn12000477'
sanger.tab<-'syn12104377'
ntap.tab<-'syn12209124'

library(synapser)
synLogin()
library(tidyverse)


getDrugOverlap<-function(tab,title=''){
  q.res<-synTableQuery(paste('select distinct "Input Drug","Output Drugs",mu,beta,w,Quantiles from',tab))$asDataFrame()
  q.res$Recovered=apply(q.res,1,function(y) y[["Input Drug"]]%in%sapply(y[["Output Drugs"]],function(x) sapply(unlist(strsplit(x,split=',')),tolower)))
  library(ggplot2)
  counts<-q.res%>%group_by(mu,w,beta,Quantiles)%>%summarize(Fraction=length(which(Recovered))/length(Recovered))

  ggplot(counts)+geom_point(aes(x=beta,y=Fraction,color=w,size=mu),stat="identity")+facet_grid(Quantiles~.)+scale_fill_viridis_d()+ggtitle(paste('Number of drugs found in',title))
  ggsave(paste(tab,'drugsRecovered.png',sep='_'))

}

getDrugOverlap(sanger.tab,'Sanger Dataset')
getDrugOverlap(ccle.tab,'CCLE Dataset')
getDrugOverlap(ntap.tab,'NTAP Dataset')
