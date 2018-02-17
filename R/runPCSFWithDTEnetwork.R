#--------------------------------------------------------
# runPCSFwithDTEnetwork
#
# take set of viper-identified proteins and merge with DTE network to run pcsf
#
#---------------------------------------------------------


library(devtools)
load_all('/Users/sgosline/code/PCSF')
#library(PCSF)
loadDrugGraph <- function(){

  ##load drug-target networ
  library(synapseClient)
  synapseLogin()
  drug.graph<-readRDS(synGet('syn11802194')@filePath)
  return(drug.graph)
}

getDrugs <-function(drug.graph){
  drugs<-names(which(degree(drug.graph,mode="out")>0))
  drugs
}


getDrugIds <- function(drug_names){

  prefix="select * from syn11819696 where common_name="
  query=paste(prefix,paste(drug_names,collapse=" OR common_name="),sep='')
  res <- synTableQuery(query)@values$internal_id
  return(res)

}

getDrugNames <- function(drug_ids){

  prefix="select * from syn11819696 where internal_id="
  query=paste(prefix,paste(drug_names,collapse=" OR internal_name="),sep='')
  res <- synTableQuery(query)@values$common_name
  return(res)
}

buildNetwork <- function(drug.graph){

  ##add weights
  edge_attr(drug.graph,"weight")<-edge_attr(drug.graph,"mean_pchembl")


  ##load STRING
  data("STRING")

  ##build interactome
  ppi <- construct_interactome(STRING)

  ##now merge networks
  combined.graph<-as.undirected(drug.graph)+ppi

  return(combined.graph)
}


