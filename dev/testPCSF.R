##test pcsf code with drug/ppi/viper combination

library(fendR)
library(devtools)
load_all('/Users/sgosline/code/PCSF')
##merge two networks into one
mergePpiDrugNets<-function(){
  data("STRING")

  ppi <- construct_interactome(STRING)

}

library(synapseClient)
synapseLogin()
drug.graph<-readRDS(synGet('syn11802194')@filePath)





