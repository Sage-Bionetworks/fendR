#--------------------------------------------------------
# runPCSFwithDTEnetwork
#
# take set of viper-identified proteins and merge with DTE network to run pcsf
#
#---------------------------------------------------------


library(devtools)
load_all('/Users/sgosline/code/PCSF')
#library(PCSF)

#' \code{loadDrugGraph} Identifies drugs in a
#' @keywords
#' @export
#' @examples
#' @return
#'
loadDrugGraph <- function(){

  ##load drug-target networ
  require(synapseClient)
  synapseLogin()
  drug.graph<-readRDS(synGet('syn11802194')@filePath)
  return(drug.graph)
}

#' \code{getDrugs} Identifies drugs in a
#' @param drug.graph
#' @keywords
#' @export
#' @examples
#' @return
#'
getDrugs <-function(drug.graph){
  drugs<-names(which(degree(drug.graph,mode="out")>0))
  drugs
}


#' \code{getDrugIds} Identifies drugs in a
#' @param drug_names
#' @param pheno.file Tidied drug response
#' @keywords
#' @export
#' @examples
#' @return
#'
getDrugIds <- function(drug_names){
  require(synapseClient)
  synapseLogin()

  #remove problematic/combo screens
  parens=grep("(",drug_names,fixed=T)
  if(length(parens)>0)
    drug_names=drug_names[-parens]

  prefix="select * from syn11831632 where common_name='"
  query=paste(prefix,paste(drug_names,collapse="' OR common_name='"),sep='')
  res <- synTableQuery(paste(query,"'",sep=''))@values

  print(paste("Found",nrow(res),'drug internal ids for',length(drug_names),'common names'))
  colnames(res) <- c("ids","drugs")

  return(res)

}

#' \code{getDrugNames} Identifies drugs in a
#' @param drug_ids
#' @keywords
#' @export
#' @examples
#' @return
#'
getDrugNames <- function(drug_ids){
  require(synapseClient)
  synapseLogin()
  prefix="select * from syn11831632 where internal_id='"
  query=paste(prefix,paste(drug_names,collapse="' OR internal_id='"),sep='')
  res <- synTableQuery(paste(query,"'",sep=''))@values
  return(res)
}


#' \code{buildNetwork} Identifies drugs in a
#' @param drug.graph Tidied drug response
#' @keywords
#' @export
#' @examples
#' @return
#'
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


#' \code{getDrugsFromGraph} Identifies drugs in a
#' @param combined.graph Graph
#' @keywords
#' @export
#' @examples
#' @return
#'
getDrugsFromGraph <-function(drug.graph){
  #the goal of this is to extract the drugs from the graph,
  #which should be the only nodes with an indegree of zero
  names(which(degree(drug.graph,mode="in")==0))

}
