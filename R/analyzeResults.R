#------------------
#
#analyze results, create browseable visualization
#do functional analysis restul
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
  res
}

getGeneExpressionEnrich<-function(tableId,synId){

}

#' \code{getViperProtExpressionEnrich} Gets viper proteins from results and computes go BP enrichment
#' @param synId synapse id of stored network result
#' @param tableId of stored fendR results
#' @keywords
#' @export
#' @examples
#' @return data frame of all go enrichment with adjusted p-value less than 0.1
getViperProtExpressionEnrich <-function(synId,tableId){
  require(synapser)
  synLogin()
  vp=synTableQuery(paste("select distinct 'Viper Proteins' from ",tableId,' where "PCSF Result"=',"'",synId,"'",sep=""))$asDataFrame()
  require(enrichR)
  res<-enrichr(unlist(strsplit(vp[1,1],split=',')),database='GO_Biological_Process_2013')
  subset(res$GO_Biological_Process_2013,Adjusted.P.value<0.1)
}

#' \code{plotDrugsAcrossData} Gets drugs from network and plots across input data
#' @param synId synapse id of stored network result
#' @param tableId of stored fendR results
#' @keywords
#' @export
#' @examples
#' @return data frame of all go enrichment with adjusted p-value less than 0.1
plotDrugsAcrossData<-function(synId,tableId){
  require(synapser)
  synLogin()
  vp=synTableQuery(paste("select distinct 'Original eSet' from ",tableId,' where "PCSF Result"=',"'",synId,"'",sep=""))$asDataFrame()

  #first get drugs responses
  require(Biobase)
  eset <-readRDS(synapser::synGet(vp[1,1])$path)
  p.data<-pData(eset)

  #now get selected drugs
  network<-readRDS(synapser::synGet(synId)$path)


}
