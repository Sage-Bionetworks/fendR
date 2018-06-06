#------------------
#
#analyze results, create browseable visualization
#do functional analysis restul
#------------------

this.script <<-'https://raw.githubusercontent.com/Sage-Bionetworks/fendR/master/R/analyzeResults.R?token=ABwyOsi9wIqVX7UHkp8AOLwoQBZ95HmGks5bIJM3wA%3D%3D'
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

  enrich<-do.call('rbind',lapply(1:length(res$enrichment),function(x){
      tres=res$enrichment[[x]]
      tres$Cluster=rep(x,nrow(tres))
      tres
  }))

  fname=paste(synId,'graphEnrichment.tsv',sep='_')
  write.table(enrich,file=fname,sep='\t')

  foldname=paste(synId,"results",sep='_')
  fold.id<-synStore(Folder(name=foldname,parentId=synGet(synId)$properties$parentId))


  synStore(File(fname,parentId=fold.id$properties$id),used=c(synId),executed=this.script)

  gname=paste(synId,'graphPlot.png',sep='_')
  png(gname)
  plot(res$subnet)
  dev.off()
  synStore(File(gname,parentId=fold.id$properties$id),used=c(synId),executed=this.script)
}

getGeneExpressionEnrich<-function(tableId,synId){
  foldname=paste(synId,"results",sep='_')
  fold.id<-synStore(Folder(name=foldname,parentId=synGet(synId)$properties$parentId))

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
  tres<-subset(res$GO_Biological_Process_2013,Adjusted.P.value<0.1)
  fname=paste('viperProteExpresionEnrichment',synId,'from',tableId,'p0.1.tsv',sep='_')
  write.table(tres,fname,sep='\t')

  foldname=paste(synId,"results",sep='_')
  fold.id<-synStore(Folder(name=foldname,parentId=synGet(synId)$properties$parentId))

  synStore(File(fname,parentId=fold.id$properties$id),used=c(synId,tableId),executed=this.script)
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
  #network<-readRDS(synapser::synGet(synId)$path)
  vd=synTableQuery(paste("select distinct 'Output Drugs' from ",tableId,' where "PCSF Result"=',"'",synId,"'",sep=""))$asDataFrame()
  od<-sapply(unlist(strsplit(vd[1,1],split=',')),tolower)
  overlap<-intersect(colnames(p.data),od)
  print(paste('found',length(overlap),'drugs that were tested alreadys out of',length(od),'in network'))

  require(ggplot2)
  #retidy up the phenotypic data
  red.p<-p.data[,c('nf1 genotype',overlap)]
  red.p$Sample <- rownames(red.p)
  red.p<-red.p%>%gather(Drug,Response,overlap)
  ggplot(red.p)+geom_boxplot(aes(x=Drug,y=Response,col=`nf1 genotype`))

  fname<-paste("drugResultsFrom",synId,'networkStoredIn',tableId,'table.png',sep='_')
  ggsave(fname)
  foldname=paste(synId,"results",sep='_')
  fold.id<-synStore(Folder(name=foldname,parentId=synGet(synId)$properties$parentId))

  res<-synStore(File(fname,parentId=fold.id$properties$id),used=c(vp[1,1],synId,tableId),executed=this.script)

  red.p
}


#' \code{doAllNetworkAssess} does all the network stuff
#' @param synId synapse id of stored network result
#' @param tableId of stored fendR results
#' @keywords
#' @export
#' @examples
#' @return
doAllNetworkAssess <-function(synId,tableId){
  getGraphSaveVis(synId)
  plotDrugsAcrossData(synId,tableId)
  getViperProtExpressionEnrich(synId,tableId)
}
