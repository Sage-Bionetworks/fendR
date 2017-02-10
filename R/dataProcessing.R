##
## These commands will collect and load sample data to put it in a format that
## can be read into the fendR class.

## -->These might be able to be loaded into the fendR S3 object to be inherited by all the
## other classes


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
  tab<-read.table(fname,stringsAsFactors =FALSE)
  tab$Gene<-rownames(tab)
  res<-tidyr::gather(tab,"Sample","Value",1:(ncol(tab)-1))
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
  tab<-read.table(fname,header=T,stringsAsFactors =FALSE)
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

  writeLines("This function reads in a fully connected matrix as 3 column edge list:\n\"geneA geneB edge\"\nand writes out a NxN feather file...  Enjoy!")

  el    <- fread(elPath, data.table = T); colnames(el) <- c("geneA","geneB","edge"); gc()              # names are needed for the plyr and tidyr calls
  el    <- as_tibble(el)
  print("spreading df to matrix")
  el    <- spread(el, key= "geneB", value = "edge"); gc()                                              # tidyr call to get matrix
  print("filling in lower tri")
  temp  <- as.matrix(el[, 2:ncol(el)]); gc()                                                           # necessary for filling in lower tri
  gName <- c(el$geneA[1],colnames(el)[-1]); rm(el);gc()
  ret   <- matrix(1,ncol(temp)+1, ncol(temp)+1); colnames(ret) <- gName

  ret[upper.tri(ret)] <- temp[upper.tri(temp,diag = T)]
  temp <- t(temp)
  ret[lower.tri(ret)] <- temp[lower.tri(temp,diag = T)]; rm(temp); gc()

  if(is.null(outPath)){outPath = dirname(elPath)}
  print(paste("writing out to", outPath))
  write_feather(as.data.frame(ret), path = file.path(outPath,gsub("txt$","feather",basename(elPath))))
}







