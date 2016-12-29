##
## These commands will collect and load data to put it in a format that
## can be read into the fendR class.

## -->These might be able to be loaded into the fendR S4 object to be inherited by all the
## other classes

#' Load Network data
#'
#' \code{loadNetwork} takes a file path and formats it as network
#' @param Path to network file
#' @keywords network feather
#' @export
#' @return ? matrix? data frame?
#' @examples
loadNetwork <- function(fname){

}

#' Load sample data
#'
#' \code{loadSampleData} takes a file path and loads it into a data frame
#' @param Path to sample file
#' @keywords
#' @export
#' @examples
#' @return Data Frame
loadSampleData <- function(fname){
  read.table(fname)
}

#' Load phenotype data
#'
#' \code{loadPhenotypeData} loads in the phenotype data with the sample identifiers as rows and any phenotype of interest as the column
#' @param Path to phenotype file
#' @keywords phenotype
#' @export
#' @examples
#' @return Data frame with at least 1 column
loadPhenotypeData <- function(fname){
  read.table(fname)
}
