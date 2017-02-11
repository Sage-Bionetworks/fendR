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
  tab<-read.table(fname,stringsAsFactors =FALSE)
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
  tab<-read.table(fname,stringsAsFactors =FALSE)
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

#' translateMatrixIdentifiers
#'
#' \code{translateMatrixIdentifiers} translates the column and row identifers based on a translation table.  This may drop and/or duplicate columns and rows.
#' @param matrix  A symmetric matrix
#' @param translationTable  A data.frame that maps from names in column 'from' to names in column 'to'
#' @param agg.fun A string name of a function that will be used to aggregrate multiple columns (from) into a single column (to).  Currently: mean, min, or max.
#' @export
#' @return A symmetric matrix with columns and rows translated (including dropped or duplicated) according to the translationTable.
translateMatrixIdentifiers <- function(matrix, translationTable, agg.fun = "mean", debug = FALSE) 
{
  if(!("from" %in% colnames(translationTable))) {
    stop("translationTable must have 'from' column\n")
  }

  if(!("to" %in% colnames(translationTable))) {
    stop("translationTable must have 'to' column\n")
  }

  supported.agg.funcs <- c("mean", "max", "min")
  if(!(agg.fun %in% supported.agg.funcs)) {
    stop(paste0("Supported aggregate functions are: ", paste0(supported.agg.funcs, collapse=","), "\n"))
  }

  ## Note that the matrix is symmetric (i.e., rows = cols)
  cols <- colnames(matrix)

  ## Drop cols and rows that do not have a translation
  cols.with.translation.flag <- (cols %in% mgi2HsMap$from)
  cols.with.translation <- cols[cols.with.translation.flag]
  cols.without.translation <- cols[!cols.with.translation.flag]
  cat(paste0(length(cols.with.translation), " columns (of ", length(cols), ") have at least one translation\n"))

  ## Since rows and columns are symmetric
  rows.with.translation.flag <- cols.with.translation.flag

  matrix <- matrix[rows.with.translation.flag, cols.with.translation.flag]
  gc()

  matrix <- as.matrix(matrix)
  gc()

  rownames(matrix) <- colnames(matrix)

  ## Subset the translation table to those entries used in the matrix
  translationTable <- translationTable[translationTable$from %in% colnames(matrix),]
  translationTable <- unique(translationTable)

  ## Create a map from col/row name to its index
  df <- data.frame(name = colnames(matrix), index = 1:ncol(matrix))

  ## Annotate the translation table with the relevant index of the entry in the 'from' column
  tbl <- merge(translationTable, df, by.x = "from", by.y = "name")

  indices <- as.numeric(tbl$index)

  agg <- NULL
  ret <- NULL
  use.matrix.utils <- FALSE
  if(use.matrix.utils) {
    suppressPackageStartupMessages(library("Matrix.utils"))
    if(agg.fun != "mean") {
      stop("aggregate.Matrix only handles fun=mean\n")
    }
    agg <- t(aggregate.Matrix(t(matrix[indices,indices]), groupings=list(tbl$to), fun=agg.fun))
    if(debug == TRUE) {
      ## Output a column in the original matrix that is duplicated in the output matrix
      tmp <- tbl[!duplicated(tbl$to, fromLast=TRUE) & !duplicated(tbl$to, fromLast=FALSE),]
      tmp <- tmp[duplicated(tmp$from, fromLast=TRUE) | duplicated(tmp$from, fromLast=FALSE),]
      tmp <- tmp[order(tmp$from),]
      tmp <- head(tmp, 2)
      cat("This single column:\n")
      print(head(matrix[,tmp$from[1],drop=FALSE],5))
      cat("Should be duplicated in each of these two columns:\n")
      print(head(agg[,(colnames(agg) %in% tmp$to),drop=FALSE],6))

      tmp <- tbl[duplicated(tbl$to, fromLast=TRUE) | duplicated(tbl$to, fromLast=FALSE),]
      tmp <- tmp[!duplicated(tmp$from, fromLast=TRUE) & !duplicated(tmp$from, fromLast=FALSE),]
      tmp <- tmp[order(tmp$to),]
      tmp <- head(tmp, 2)
      cat("\nThese two columns:\n")
      print(head(matrix[,tmp$from,drop=FALSE],5))
      cat(paste0("Should be summarized by ", agg.fun, " into this single column:\n"))
      print(head(agg[,(colnames(agg) %in% tmp$to),drop=FALSE],6))
    }
    ret <- aggregate.Matrix(agg, groupings=list(tbl$to), fun=agg.fun)
  } else {
    suppressPackageStartupMessages(library("dplyr"))
    df <- as.data.frame(matrix[indices, indices])
    rownames(df) <- 1:nrow(df)
    colnames(df) <- rownames(df)
    df$group <- tbl$to
    agg <- t((df %>% group_by(group) %>% summarise_each(funs_(agg.fun))))
    if(debug == TRUE) {
      ## Output a column in the original matrix that is duplicated in the output matrix
      tmp <- tbl[!duplicated(tbl$to, fromLast=TRUE) & !duplicated(tbl$to, fromLast=FALSE),]
      tmp <- tmp[duplicated(tmp$from, fromLast=TRUE) | duplicated(tmp$from, fromLast=FALSE),]
      tmp <- tmp[order(tmp$from),]
      tmp <- head(tmp, 2)
      cat("This single column:\n")
      print(head(matrix[,tmp$from[1],drop=FALSE],5))
      cat("Should be duplicated in each of these two columns:\n")
      print(head(agg[,(agg[1,] %in% tmp$to),drop=FALSE],6))

      tmp <- tbl[duplicated(tbl$to, fromLast=TRUE) | duplicated(tbl$to, fromLast=FALSE),]
      tmp <- tmp[!duplicated(tmp$from, fromLast=TRUE) & !duplicated(tmp$from, fromLast=FALSE),]
      tmp <- tmp[order(tmp$to),]
      tmp <- head(tmp, 2)
      cat("\nThese two columns:\n")
      print(head(matrix[,tmp$from,drop=FALSE],5))
      cat(paste0("Should be summarized by ", agg.fun, " into this single column:\n"))
      print(head(agg[,(agg[1,] %in% tmp$to),drop=FALSE],6))
    }
    colnames(agg) <- agg[1,]
    agg <- (agg[-1,])
    class(agg) <- "numeric"
    df <- as.data.frame(agg)
    df$group <- tbl$to
    ret <- ((df %>% group_by(group) %>% summarise_each(funs_(agg.fun))))
    labels <- ret$group
    ret <- ret[,-1]
    ret <- as.matrix(ret)
    class(ret) <- "numeric"
    rownames(ret) <- labels
    colnames(ret) <- rownames(ret)
  }

  ## Confirm that the resulting matrix is symmetric
  if((nrow(ret) != ncol(ret)) || (!(all(rownames(ret) == colnames(ret))))) {
    stop("Translated matrix is not symmetric!\n")
  }

  ret
}






