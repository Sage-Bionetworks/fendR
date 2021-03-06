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
  this.script <<-'https://raw.githubusercontent.com/Sage-Bionetworks/fendR/master/R/analyzeResults.R?token=ABwyOsi9wIqVX7UHkp8AOLwoQBZ95HmGks5bIJM3wA%3D%3D'

  require(synapser)
  synapser::synLogin()
  sf<-readRDS(synapser::synGet(synId)$path)

  require(PCSF)
  res<-PCSF::enrichment_analysis(sf)

  enrich<-do.call('rbind',lapply(1:length(res$enrichment),function(x){
      tres=res$enrichment[[x]]
      tres$Cluster=rep(x,nrow(tres))
      tres
  }))

  fname=paste(synId,'graphEnrichment.tsv',sep='_')
  write.table(enrich,file=fname,sep='\t')

  foldname=paste(synId,"results",sep='_')
  fold.id<-synapser::synStore(Folder(name=foldname,parentId=synapser::synGet(synId)$properties$parentId))


  synStore(File(fname,parentId=fold.id$properties$id),used=c(synId),executed=this.script)

  gname=paste(synId,'graphPlot.html',sep='_')
  #plot(res$subnet)
  #png(filename=gname)
  net<-plot(res$subnet)
  visNetwork::visSave(net,gname)
  #dev.off()
  synapser::synStore(File(gname,parentId=fold.id$properties$id),used=c(synId),executed=this.script)
}

getGeneExpressionEnrich<-function(tableId,synId){
  foldname=paste(synId,"results",sep='_')
  fold.id<-synapser::synStore(Folder(name=foldname,parentId=synapser::synGet(synId)$properties$parentId))

}

#' \code{getViperProtExpressionEnrich} Gets viper proteins from results and computes go BP enrichment
#' @param synId synapse id of stored network result
#' @param tableId of stored fendR results
#' @keywords
#' @export
#' @examples
#' @return data frame of all go enrichment with adjusted p-value less than 0.1
getViperProtExpressionEnrich <-function(synId,tableId){
  this.script <<-'https://raw.githubusercontent.com/Sage-Bionetworks/fendR/master/R/analyzeResults.R?token=ABwyOsi9wIqVX7UHkp8AOLwoQBZ95HmGks5bIJM3wA%3D%3D'

  require(synapser)
  synLogin()
  vp=synapser::synTableQuery(paste("select distinct 'Viper Proteins' from ",tableId,' where "PCSF Result"=',"'",synId,"'",sep=""))$asDataFrame()
  require(enrichR)
  res<-enrichr(unlist(strsplit(vp[1,1],split=',')),database='GO_Biological_Process_2013')
  tres<-subset(res$GO_Biological_Process_2013,Adjusted.P.value<0.1)
  fname=paste('viperProteExpresionEnrichment',synId,'from',tableId,'p0.1.tsv',sep='_')
  write.table(tres,fname,sep='\t')

  foldname=paste(synId,"results",sep='_')
  fold.id<-synStore(Folder(name=foldname,parentId=synapser::synGet(synId)$properties$parentId))

  synapser::synStore(File(fname,parentId=fold.id$properties$id),used=c(synId,tableId),executed=this.script)
}

#' \code{plotDrugsAcrossData} Gets drugs from network and plots across input data
#' @param synId synapse id of stored network result
#' @param tableId of stored fendR results
#' @keywords
#' @export
#' @examples
#' @return data frame of all go enrichment with adjusted p-value less than 0.1
plotDrugsAcrossData<-function(synId,tableId,genotype='nf1 genotype'){
  this.script <<-'https://raw.githubusercontent.com/Sage-Bionetworks/fendR/master/R/analyzeResults.R?token=ABwyOsi9wIqVX7UHkp8AOLwoQBZ95HmGks5bIJM3wA%3D%3D'

  require(synapser)
  synLogin()
  vp=synTableQuery(paste("select distinct 'Original eSet' from ",tableId,' where "PCSF Result"=',"'",synId,"'",sep=""))$asDataFrame()

  #first get drugs responses
  require(Biobase)
  eset <-readRDS(synapser::synGet(vp[1,1])$path)
   #now get selected drugs
  #network<-readRDS(synapser::synGet(synId)$path)
  vd=synTableQuery(paste("select distinct 'Output Drugs' from ",tableId,' where "PCSF Result"=',"'",synId,"'",sep=""))$asDataFrame()
  od<-sapply(unlist(strsplit(vd[1,1],split=',')),tolower)

  res=plotDrugs(eset,od,genotype)
  fname<-paste("drugResultsFrom",synId,'networkStoredIn',tableId,'table.png',sep='_')
  ggsave(fname)
  foldname=paste(synId,"results",sep='_')
  fold.id<-synapser::synStore(Folder(name=foldname,parentId=synapser::synGet(synId)$properties$parentId))

  res<-synapser::synStore(File(fname,parentId=fold.id$properties$id),used=c(vp[1,1],synId,tableId),executed=this.script)

}

#' \code{plotDrugs} carries out the plotting mechanism for  a single network
#' @param eset expressionset
#' @param drugList list of drugs
#' @export
plotDrugs <-function(eset.file,drugList,genotype){
  eset<-readRDS(synapser::synGet(eset.file)$path)
  p.data<-pData(eset)
  dnames=toupper(sapply(colnames(p.data),function(x) unlist(strsplit(x,split='_'))[1]))
  overlap<-intersect(dnames,toupper(drugList))

  print(paste('found',length(overlap),'drugs that were tested alreadys out of',length(drugList),'in network'))


  if(length(overlap)==0)
    return(NULL)
  require(ggplot2)
  #retidy up the phenotypic data
  dd.names<-colnames(p.data)[which(dnames%in%overlap)]
  red.p<-p.data[,c(genotype,dd.names)]
  red.p$Sample <- rownames(red.p)
  red.p<-red.p%>%gather(Drug,Response,dd.names)
  red.p$mutation<-as.factor(red.p[[genotype]])
  red.p=subset(red.p,!is.na(mutation))
  ggplot(red.p)+geom_boxplot(aes(x=mutation,y=Response,col=mutation))+facet_grid(.~Drug)+ggtitle(paste(genotype,'mutation by drug response'))
  fname=paste(eset.file,paste(drugList,collapse='_'),genotype,'genotype.png',sep='_')
  ggsave(filename=fname)
  tab<-red.p%>%group_by(Drug)%>%do(broom::tidy(wilcox.test(Response~mutation,.,na.rm=T)))%>%dplyr::select(Drug,p.value)
  tab$figFile=rep(fname,nrow(tab))
  tab
  }


#' \code{doAllNetworkAssess} does all the network stuff
#' @param synId synapse id of stored network result
#' @param tableId of stored fendR results
#' @keywords
#' @export
#' @examples
#' @return
doAllNetworkAssess <-function(synId,tableId){
  print(paste('assessing',synId,'and',tableId))
  fendR::getGraphSaveVis(synId)
  fendR::plotDrugsAcrossData(synId,tableId)
  fendR::getViperProtExpressionEnrich(synId,tableId)
}
