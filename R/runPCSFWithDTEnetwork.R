#--------------------------------------------------------
# runPCSFwithDTEnetwork
#
# take set of viper-identified proteins and merge with DTE network to run pcsf
#
#---------------------------------------------------------


require(fendR)
library(devtools)
load_all('/Users/sgosline/code/PCSF')

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

getDrugNames <- function(drugs){

  mapping <- readRDS(synGet('syn11802195')@filePath)
  mapping
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


