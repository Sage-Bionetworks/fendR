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








#' edgeList2matrix convertion and write out
#'
#' \code{edgeList2matrix} loads in the target/phenotype data with the first column representing the drug of choice and the second represneting the gene
#' @param elPath file path to edge list file with three colums, geneA, geneB and edge
#' @param outType bigMatrix or feather
#' @param outPath directory for output if NULL it will be same as input file
#' @import data.table tibble tidyr plyr bigmemory feather
#' @return Data frame with at least 2 column
edgeList2matrix =function(elPath, outType = "bigMatrix", outPath=NULL) # el (edge list) should be a data.tableish represenation of an edgelist with 1st 2 colums being gene names
{
  suppressPackageStartupMessages(library("data.table"))	
  suppressPackageStartupMessages(library("tibble"))	
  suppressPackageStartupMessages(library("tidyr"))	
  suppressPackageStartupMessages(library("plyr"))	
  el    <- fread(elPath, data.table = T); colnames(el) <- c("geneA","geneB","edge"); gc()              # names are needed for the plyr and tidyr calls
  el    <- as_tibble(el)
  el    <- arrange(el, geneA, geneB); gc()                                                             # need it sorted before make into symetric matrix
  el    <- spread(el, key= "geneB", value = "edge"); gc()                                              # tidyr call to get matrix
  temp  <- as.matrix(el[, 2:ncol(el)]); gc()                                                           # necessary for filling in lower tri 
  gName <- c(el$geneA[1],colnames(el)[-1]); rm(el);gc()
  ret   <- matrix(1,ncol(temp)+1, ncol(temp)+1); colnames(ret) <- gName
  
  ret[upper.tri(ret)] <- temp[upper.tri(temp,diag = T)]
  temp <- t(temp)
  ret[lower.tri(ret)] <- temp[lower.tri(temp,diag = T)]; rm(temp); gc()
  
  if(is.null(outPath)){outPath = dirname(elPath)}
  if(outType =="bigMatrix")
  {
    suppressPackageStartupMessages(library("bigmemory"))	
    write.big.matrix(as.big.matrix(ret), filename = file.path(outPath,gsub("\\.txt$","BM\\.txt",basename(elPath))))
  }
  if(outType =="feather")
  {  
    suppressPackageStartupMessages(library("feather"))	
    write_feather(as.data.frame(ret), path = file.path(outPath,gsub("txt$","feather",basename(elPath))))
  }
}







