#------------------
#
#analyze results, create browseable visualization
#------------------

#' \code{getGraphSaveVis} Grabs a graph object stored on synapse, runs enrichment and saves visualiation
#' @param synId synapse id of file
#' @keywords
#' @export
#' @examples
#' @return not sure
getGraphSaveVis <- function(synId){
  require(synapser)
  synLogin()
  file<-readRDS(synapser::synGet(synId)$path)

  require(PCSF)
  res<-enrichment_analysis(file)
}
