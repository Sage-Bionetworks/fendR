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
#' @param outPath directory for output if NULL it will be same as input file
#' @import data.table tibble tidyr plyr bigmemory feather
#' @export
#' @return NULL if successfull,error otherwise
edgeList2matrix =function(elPath, outPath=NULL) # el (edge list) should be a data.tableish represenation of an edgelist with 1st 2 colums being gene names
{
  suppressPackageStartupMessages(library("data.table"))
  suppressPackageStartupMessages(library("tibble"))
  suppressPackageStartupMessages(library("tidyr"))
  suppressPackageStartupMessages(library("plyr"))
  suppressPackageStartupMessages(library("feather"))
  
  writeLines("This function reads in a fully connected matrix as 4 column edge list:\n\"geneA geneB goldStandardFlag edge\"\n\nlike those found at http://fntm.princeton.edu/download/\nand writes out a NxN feather file\nit can take a while (like 30 min on a r3.4xlarge ec2 instance)...\nEnjoy!\n\nP.S. you may want to clean up your environment before running this")
  
  el    <- data.table::fread(elPath, data.table = T); colnames(el) <- c("geneA","geneB","goldStandardFlag", "edge"); gc() # names are needed for the plyr and tidyr calls
  el    <- el[,-3]; gc()
  el    <- tibble::as_tibble(el); gc()
  elRev <- el; elRev$geneA <- el$geneB; elRev$geneB <- el$geneA
  el    <- base::rbind(el, elRev); rm(elRev); gc();
  
  print("spreading df to matrix")
  print(date())
  el       <- tidyr::spread(el, key= "geneB", value = "edge"); gc() # tidyr call to get tible matrix
  el       <- as.matrix(el[,-1])
  diag(el) <- 1
  print(date())
  
  if(is.null(outPath)){outPath = dirname(elPath)}
  print(paste("writing out to", outPath))
  feather::write_feather(as.data.frame(el), path = file.path(outPath,gsub("txt$","feather",basename(elPath))))
}







