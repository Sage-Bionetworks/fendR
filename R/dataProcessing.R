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
#' @examples
#' @return tidied data frame with columns 'Gene','Sample' and 'Value'
loadSampleData <- function(fname){
  if(length(grep("gz$",fname))>0)
    tab<-read.table(gzfile(fname),stringsAsFactors =FALSE,check.names=F,header=T)
  else
    tab<-read.table(fname,stringsAsFactors =FALSE,check.names=F,header=T)
  tab$Gene<-rownames(tab)
  res<-NULL
  try(
    res<-tidyr::gather(data.frame(tab,check.names=F),"Sample","Value",1:(ncol(tab)-1))
  )
  if(is.null(res))
    res<-tidyr::gather(data.frame(tab,check.names=T),"Sample","Value",1:(ncol(tab)-1))

  return(res)
}

#' Load phenotype data
#'
#' \code{loadPhenotypeData} loads in the phenotype data with the sample identifiers as rows and any phenotype of interest as the column
#' @param Path to phenotype file
#' @keywords phenotype
#' @export
#' @examples
#' @return tidied data frame with columns 'Sample','Phenotype','Response'
loadPhenotypeData <- function(fname){
  tab<-read.table(fname,sep ='\t',check.names=F)
  tab$Sample<-rownames(tab)
  res<-tidyr::gather(tab,"Phenotype","Response",1:(ncol(tab)-1))
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
  tab<-read.table(fname,header=T,sep='\t',stringsAsFactors =FALSE,check.names=F)
  if(ncol(tab)==2)
    colnames(tab)<-c("Phenotype","Gene")
  else
    tab<-tab%>%select(Phenotype,Gene)
  tab
}



#' edgeList2matrix convertion and write out
#'
#' \code{edgeList2matrix} loads in the target/phenotype data with the first column representing the drug of choice and the second represneting the gene
#' @param elPath file path to edge list file with three colums, geneA, geneB and edge
#' @param outPath directory for output if NULL it will be same as input file
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

#' drugNameMapping
#'
#' \code{drugNameMapping} collects information from the CCLE, SANGER and NTAP datasets to create a table that maps one drug to the others
#' @param
#' @param
#' @export
#' @return a data frame with each drug and its aliass on each row
drugNameMapping <- function(){}

#' cellLineMapping
#'
#' \code{cellLineMapping} collects information from CCLE and SANGER datasets to ensure to create a table that maps cell lines from one to another
#'
cellLineMapping <-function(){

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
  cols.with.translation.flag <- (cols %in% translationTable$from)
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






