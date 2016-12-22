#this file creates examples that run fendR

library(fendR)

#load in network, files from system.file included in package
gene.data<-loadSampleData(system.file('',package='fendR'))
pheno.data<-loadPhenoData(system.file('',package='fendR'))
network.data<-loadNetwork(system.file('',package='fendR'))

#create new forest class with data
forestObj<-forestFendr(gene.data,pheno.data,network.data)

#do some LOO magic
