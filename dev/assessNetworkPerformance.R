##-------------------
## assessNetworkPerformance
## grabs the results of the PCSF algorithm and plots various performance measures
##
##-------------------


require(synapser)
require(ggplot2)
library(fendR)
require(tidyverse)
library(rJava)
library(rcdk)
library(fingerprint)
synLogin()
this.script='https://raw.githubusercontent.com/Sage-Bionetworks/fendR/master/dev/assessNetworkPerformance.R?token=ABwyOjNTIyrNV8HIeC66DIGnQWX90Uv1ks5a1M1-wA%3D%3D'

#'
#'we want to compute the tanimoto distance between the proposed target
#'and the targets that were identified
#'

is.smiles <- function(x, verbose = TRUE) { ##corrected version from webchem
  if (!requireNamespace("rcdk", quietly = TRUE)) {
    stop("rcdk needed for this function to work. Please install it.",
         call. = FALSE)
  }
  # x <- 'Clc(c(Cl)c(Cl)c1C(=O)O)c(Cl)c1Cl'
  if (length(x) > 1) {
    stop('Cannot handle multiple input strings.')
  }
  out <- try(rcdk::parse.smiles(x), silent = TRUE)
  if (inherits(out[[1]], "try-error") | is.null(out[[1]])) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

parseInputFingerprint <- function(input, fp.type) {
  if(is.smiles(input)==TRUE){
    input.mol <- parse.smiles(as.character(input))
    lapply(input.mol, do.typing)
    lapply(input.mol, do.aromaticity)
    lapply(input.mol, do.isotopes)
    fp.inp <- lapply(input.mol, get.fingerprint, type = fp.type)
  }else{
    print('Please input a valid SMILES string.')
  }
}

distance.minified <- function(fp1,fp.list){
  n <- length(fp1)
  f1 <- numeric(n)
  f2 <- numeric(n)
  f1[fp1@bits] <- 1

  sapply(fp.list, function(x){
    f2[x@bits] <- 1
    sim <- 0.0
    ret <- .C("fpdistance", as.double(f1), as.double(f2),
              as.integer(n), as.integer(1),
              as.double(sim),
              PACKAGE="fingerprint")
    return (ret[[5]])
  })
}

library(pbapply)

computeTMDistance <-function(synTableId="syn12000477",parId='syn12104372',drugMap ='syn11831632'){
  tab.res <-synapser::synTableQuery(paste("select `Input Drug`,`Output Drugs`,w,beta,mu,Quantiles from",synTableId))$asDataFrame()

  name.to.dte <- readRDS(synapser::synGet("syn11712154", version = 7)$path) %>%
    mutate(common_name = tolower(common_name)) %>%
    select(common_name, smiles, internal_id) %>%
    group_by(internal_id) %>%
    top_n(1) %>%
    mutate(count = n()) %>%
    slice(1) %>%
    select(-count)

  structure.map <- readRDS(synGet("syn11712148", version=9)$path)

  sapply(tab.res$ROW_ID, function(x){
    foo <- filter(tab.res, ROW_ID == x)
    input <- foo$`Input Drug`
    print(input)
    input.dte <- filter(name.to.dte, common_name == input) %>%
      group_by(common_name) %>%
      top_n(1, internal_id)

    input.dte <- input.dte$smiles
    input.fp <- parseInputFingerprint(input.dte, "circular")

    output <- strsplit(foo$`Output Drugs`, ",")[[1]] %>% tolower()

    output.dte <- filter(name.to.dte, common_name %in% output) %>%
      group_by(common_name) %>%
      top_n(1, internal_id)

    output.fp <- sapply(output.dte$smiles, parseInputFingerprint, "circular")

    sim <- distance.minified(input.fp[[1]], output.fp)
    mean.tanimoto <- mean(sim)
    sd.tanimoto <- sd(sim)
  })
}

res <- computeTMDistance()

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




