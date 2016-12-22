#this file creates examples that run fendR

library(fendR)

#load in network, files from system.file included in package
gene.data<-loadSampleData(system.file('CCLE_binary_mutation_matrix_ucscGenesFromCBioPortal.tsv',package='fendR'))
pheno.data<-loadPhenotypeData(system.file('CTRP_v20_AUC_vales_by_drug.tsv',package='fendR'))
network.data<-loadNetwork(system.file('',package='fendR'))

#create new forest class with data
forestObj<-forestFendr(gene.data,pheno.data,network.data)

#do some LOO magic
