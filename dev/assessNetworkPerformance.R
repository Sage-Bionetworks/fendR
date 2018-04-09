##-------------------
## assessNetworkPerformance
## grabs the results of the PCSF algorithm and plots various performance measures
##
##-------------------


require(synapser)
require(ggplot2)
library(fendR)

#'
#'we want to compute the tanimoto distance between the proposed target
#'and the targets that were identified
computeTMDistance <-function(synTableId){


}


#'
#' we should also compute the overlap of predicted targets
#'
computeTargetOverlap <-function(synTableId){

}



#' as a more general summary - what drugs are showing up in which networks?
#'
drugDistributionByParameters <- function(synTableId="syn12000477"){
#what drugs are being selected?
  drug.tab <-getSelectedDrugByParameter(synTableId)
  new.res<- drug.tab%>%unite("Params", c(mu,beta,w,Quantiles),sep='_')
  dcounts<-new.res %>%group_by(Params,`Output Drugs`)%>%summarize(`NumDrugs`=n())

  ggplot(dcounts)+geom_bar(aes(x=`Output Drugs`,y=NumDrugs,fill=Params),stat='identity',position='dodge')+scale_fill_viridis_d()+theme(axis.text.x=element_text(angle=90,hjust=1))
}

#'
#'general function to get tidied table of drug names by parameters
getSelectedDrugByParameter <-function(synTableId="syn12000477"){
  require(tidyr)
  synapser::synLogin()
  tab.res <-synapser::synTableQuery(paste("select `Output Drugs`,w,beta,mu,Quantiles from",synTableId))$asDataFrame()

  new.res <- tab.res %>%mutate(`Output Drugs`=strsplit(as.character(`Output Drugs`),','))%>%unnest()
  return(new.res)

}




