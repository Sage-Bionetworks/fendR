
##this is the forest implemetnation of the fendR class

######################################################################
# Create the forestFendR class
#
# This is used to represent a forest-based implementation of the fendR framework

#' An S3 class to represent a Forest-based implementation of the fendR predictive
#' network algorithm
#' @inheritParams fendR
#' @param forestPath Path to python forest scripts
#' @export
#' @return forestFendR object
forestFendR<-function(network, featureData, phenoData, forestPath){
 me <-fendR(network,featureData,phenoData)
 me$forestPath <- forestPath
 class(me) <- append(class(me),'forestFendR')
 return(me)
}

######################################################################
# Set methods implemented by forestFendR

#' Engineer Features from Network
#' \code{createNewFeaturesFromNetwork} takes the gene-based measurements and alters their #' score using a network
#' @param object That contains a data frame and network
#' @keywords
#' @export
#' @return data.frame representing new gene by sample matrix that is augmented
createNewFeaturesFromNetwork.forestFendR<-function(object){
    ##summarize phenotypic features across samples
    geneSums<-rowSums(object$featureData)
    geneSums<-subset(geneSums,geneSums>0)
    #these become our nodeWeights for forest

    #write network, weight, config file to working directory
    prizeFileName='prizes.txt'
    write.table(geneSums,prizeFileName,quote=F,col.names=F)
    netFileName=object$network

    #then call forest.py with arguments


    #optimize arguments for number of trees and size.
    mu.range<-seq(0.001,0.004,by=0.001)
    beta.range<-seq(100,200,by=25)
    w.range<-seq(1,5,by=1)

    #iterate over all possible parameters
    all.graphs<-lapply(as.character(mu.range),function(mu){
      res<-lapply(as.character(beta.range),function(beta,mu){
        res<-lapply(as.character(w.range),function(w,beta,mu){
            runForestWithParams(object$forestPath,mu,beta,w,prizeFileName,netFileName,depth=10)

        },beta,mu)
        names(res)<-as.character(w.range)
        return(res)
      },mu)
      names(res)<-as.character(beta.range)
      return(res)
      })
    names(all.graphs)<-as.character(mu.range)


          #assign new features based on shortest path, using igraph
  return(all.graphs)

}
####################
#below are unexported helper functions to make processing forest easier
#in the future we can replace these with C calls to msgsteiner.
#' Run Forest using python code
#' \code{runForestWithParams} runs python code and returns graph objects
#' @param object That contains a data frame and network
#' @keywords
#' @return list of graph objects
runForestWithParams <- function(forestPath,mu,beta,w,prizeFileName,netFileName,depth=10){

  #create a tmp directory
  paramstr=paste('mu',mu,'beta',beta,'w',w,sep='_')
  dirname=paste(forestPath,'fendROutput',sep='/')

  if(!dir.exists(dirname))
    dir.create(dirname)

  #write conf file
  cf=file(paste(dirname,'conf.txt',sep='/'))
  writeLines(c(paste('w =',w),paste('b =',beta),paste('D =',depth),paste('mu =',mu)),cf)
  close(cf)

  cmd=paste('/usr/local/bin/python ',forestPath,'/scripts/forest.py --prize ',prizeFileName,
    ' --edge ',netFileName,
    ' --conf ',paste(dirname,'conf.txt',sep='/'),
    ' --outpath ',dirname,
    ' --outlabel ',paramstr,
    ' --msgpath ',forestPath,'/msgsteiner',
    sep='')

  #run code
  print(paste('Running forest:',cmd))
  res <- 1
  try(res<-system(cmd))

  if(res==0){
  #collect files and load into igraph -- which files do we want?
  net <- read.table(paste(dirname,'/',paramstr,'_optimalForest.sif',sep=''),sep='\t')
  nodes <-unique(c(nodes[,1],nodes[,3]))

  #return igraphs
  opt_graph <- graph_from_edgelist(net[,c(1,3)],directed = FALSE)
  print(paste("Returning graph with",length(E(opt_graph)),'edges and',length(V(opt_graph)),'nodes'))
  return(opt_graph)
  }else{
    print("Unable to run Forest, returning empty graph")
    return(empty_graph(0,directed=FALSE))
  }
}



augmentFeaturesWithGraph<-function(object,graph){

}

#' \code{scoreDataFromModel} takes the new model and predicts a phenotype from an input set
#' @param model
#' @param unseenFeatures
#' @keywords
#' @export
#' @return list of scores for each of the columns of the unseen feature data frame
scoreDataFromModel.forestFendR<-function(object){

}


