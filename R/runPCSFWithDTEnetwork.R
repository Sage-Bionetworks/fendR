#--------------------------------------------------------
# runPCSFwithDTEnetwork
#
# take set of viper-identified proteins and merge with DTE network to run pcsf
#
#---------------------------------------------------------

library(PCSF)
##make sure to install PCSF from her
##library(devtools)
##install_github("sgosline/PCSF",username='sgosline')

#' \code{loadDrugGraph} Identifies drugs in a
#' @keywords
#' @export
#' @examples
#' @return
#'
loadDrugGraph <- function(){

  ##load drug-target networ
  require(synapser)
  synLogin()
  drug.graph<-readRDS(synGet('syn11802194')$path)
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
  require(synapser)
  synLogin()
  require(dplyr)
  #remove problematic/combo screens
  parens=grep("(",drug_names,fixed=T)
  if(length(parens)>0)
    drug_names=drug_names[-parens]

  prefix="select * from syn11831632 where common_name='"
  query=paste(prefix,paste(drug_names,collapse="' OR common_name='"),sep='')
  res <- synapser::synTableQuery(paste(query,"'",sep=''))$asDataFrame()%>%dplyr::select(-ROW_ID,-ROW_VERSION)

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
  require(synapser)
  require(dplyr)
  synLogin()
  prefix="select * from syn11831632 where internal_id='"
  query=paste(prefix,paste(drug_ids,collapse="' OR internal_id='"),sep='')
  res <- synapser::synTableQuery(paste(query,"'",sep=''))$asDataFrame()%>%dplyr::select(-ROW_ID,-ROW_VERSION)
  colnames(res) <- c("ids","drugs")

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
#  edge_attr(drug.graph,"weight")<-edge_attr(drug.graph,"mean_pchembl")
  require(PCSF)

  rank.norm<-function(x){
    rank(x)/length(x)
  }

  red.graph<-delete_edges(drug.graph,E(drug.graph)[which(is.na(E(drug.graph)$mean_pchembl))])
  edge_attr(red.graph,"weight")<-rank.norm(edge_attr(red.graph,"mean_pchembl"))

  ##load STRING
  data("STRING")


  ##build interactome
  ppi <- construct_interactome(STRING)

  ##now merge networks
  u.drug.graph <- as.undirected(red.graph)
  combined.df<-rbind(igraph::as_data_frame(ppi),igraph::as_data_frame(u.drug.graph)[,c('from','to','weight')])
  combined.graph <- igraph::graph_from_data_frame(combined.df,directed=FALSE)


  return(combined.graph)
}


#' \code{getDrugsFromGraph} Identifies drugs in a
#' @param combined.graph Graph
#' @keywords
#' @export
#' @examples
#' @return drug names
#'
getDrugsFromGraph <-function(drug.graph){
  #the goal of this is to extract the drugs from the graph,
  #which should be the only nodes with an indegree of zero
  names(which(degree(drug.graph,mode="in")==0))

}

#' \code{runPcsfWithParams} Identifies drugs in a
#' @param combined.graph
#' @param terminals
#' @param dummies
#' @param w
#' @param b
#' @param mu
#' @param doRand
#' @keywords
#' @export
#' @examples
#' @return
runPcsfWithParams <- function(ppi,terminals, dummies, w=2, b=1, mu=5e-04,doRand=FALSE){

  if(doRand)
    res <- PCSF_rand(ppi,terminals,w=w,b=b,mu=mu,dummies=dummies,n=100,r=1)
  else
    res <- PCSF(ppi,terminals,w=w,b=b,mu=m,dummies=dummies)

  drug.inds<-which(V(res)$name%in%dummies)
  V(res)$type[drug.inds]<-'Compound'
  V(res)$name[drug.inds]<-getDrugNames(V(res)$name[drug.inds])[,2]
  return(res)

}


