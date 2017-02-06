##
## These commands will collect and load sample data to put it in a format that
## can be read into the fendR class.

## -->These might be able to be loaded into the fendR S3 object to be inherited by all the
## other classes


#' Load Network data
#'
#' \code{loadNetwork} takes a file path and formats it as network
#' @param Path to network file
#' @keywords network feather
#' @export
#' @import igraph
#' @return an iGraph object where edge weights represent distance between nodes (smaller means *more* association)
#' @examples
loadNetwork <- function(fname){
  library(igraph)
  tab<-read.table(fname)
  net<-graph_from_data_frame(tab,directed=F)
  E(net)$weight<-1-min(tab[,3],1)
  return(net)
}

#' Load sample data
#'
#' \code{loadSampleData} takes a file path and loads it into a data frame
#' @param Path to sample file
#' @keywords
#' @export
#' @import tidyr
#' @examples
#' @return tidied data frame with columns 'Gene','Sample' and 'Value'
loadSampleData <- function(fname){
  library(tidyr)
  tab<-read.table(fname)
  tab$Gene<-rownames(tab)
  res<-gather(tab,"Sample","Value",1:(ncol(tab)-1))
  return(res)
}

#' Load phenotype data
#'
#' \code{loadPhenotypeData} loads in the phenotype data with the sample identifiers as rows and any phenotype of interest as the column
#' @param Path to phenotype file
#' @keywords phenotype
#' @export
#' @import tidyr
#' @examples
#' @return tidied data frame with columns 'Sample','Phenotype','Response'
loadPhenotypeData <- function(fname){
  library(tidyr)
  tab<-read.table(fname)
  tab$Sample<-rownames(tab)
  res<-gather(tab,"Phenotype","Response",1:(ncol(tab)-1))
  return(res)
}

#' Load target/snp data
#'
#' \code{loadTargetData} loads in the target/phenotype data with the first column representing the drug of choice and the second represneting the gene
#' @param Path to  file
#' @keywords phenotype
#' @export
#' @examples
#' @return Data frame with at least 2 column
loadTargetData <- function(fname){
  tab<-read.table(fname,header=T)
  colnames(tab)<-c("Phenotype","Gene")
  tab
}
