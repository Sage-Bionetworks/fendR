##-------------------
## assessNetworkPerformance
## grabs the results of the PCSF algorithm and plots various performance measures
##
##-------------------


require(synapser)
require(ggplot2)
library(fendR)
require(tidyverse)
this.script='https://raw.githubusercontent.com/Sage-Bionetworks/fendR/master/dev/assessNetworkPerformance.R?token=ABwyOjNTIyrNV8HIeC66DIGnQWX90Uv1ks5a1M1-wA%3D%3D'

#'
#'we want to compute the tanimoto distance between the proposed target
#'and the targets that were identified
computeTMDistance <-function(synTableId="syn12000477",parId='syn12104372'){


}


#'
#' we should also compute the overlap of predicted targets
#'
computeTargetOverlap <-function(synTableId="syn12000477",parId='syn12104372'){

    drug.targs<-getDrugTargets()
}



#' as a more general summary - what drugs are showing up in which networks?
#'
drugDistributionByParameters <- function(synTableId="syn12000477",parId='syn12104372'){
#what drugs are being selected?
  drug.tab <-getSelectedDrugByParameter(synTableId)
  new.res<- drug.tab%>%unite("Params", c(mu,beta,w,Quantiles),sep='_')
  dcounts<-new.res %>%group_by(Params,`Output Drugs`)%>%summarize(`Times Selected`=n())
  icounts <- new.res %>% group_by(Params)%>% summarize(NumInputs=n_distinct(`Input Drug`))%>%inner_join(dcounts,by='Params') %>% mutate(FracSelected=`Times Selected`/NumInputs)

  fname=paste('drugsSelectedByparameter_',synTableId,'.png',sep='')
  p<-ggplot(icounts)+geom_bar(aes(x=`Output Drugs`,y=`FracSelected`,fill=Params),stat='identity',position='dodge')+scale_fill_viridis_d()+theme(axis.text.x=element_text(angle=90,hjust=1))
  ggsave(fname,p,limitsize=FALSE,width = par("din")[1],
    height = par("din")[2])
  synapser::synStore(File(fname,parentId=parId),used=synTableId,executed=this.script)
}

#'
#'general function to get tidied table of drug names by parameters
#'@param synTableID table id of pcsf output to collect
#'@return tidied data frame of drug information
#'@export
getSelectedDrugByParameter <-function(synTableId="syn12000477"){
  require(tidyr)
  synapser::synLogin()
  tab.res <-synapser::synTableQuery(paste("select `Input Drug`,`Output Drugs`,w,beta,mu,Quantiles from",synTableId))$asDataFrame()

  new.res <- tab.res %>%mutate(`Output Drugs`=strsplit(as.character(`Output Drugs`),','))%>%unnest()%>%dplyr::select(-ROW_ID,-ROW_VERSION)%>%unique()
  return(new.res)

}




