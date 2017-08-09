###build predictive model of dataset using just regression/rf techniques
library(fendR)


##goal is to try
ntapPerf<-function(use.rna=FALSE){
  synapseLogin()
  ##should we load the data or not? seems like a waste of time at this point
  network.file<-'https://github.com/fraenkel-lab/OmicsIntegrator/raw/master/data/iref_mitab_miscore_2013_08_12_interactome.txt'

  gene.data<-read.table(synGet('syn8304629')@filePath,header=T)
  #gene.data<-loadSampleData(gene.file)

 # rna.seq.data<-
  rna.data<-read.table(synGet('syn8304627')@filePath,header=T)

  pheno.data<-read.table(synGet('syn8304620')@filePath,header=T)
  pheno.data<-pheno.data[-which(is.na(pheno.data$Phenotype)),]

#  pheno.data<-loadPhenotypeData(pheno.file)

  target.data<-read.table(synGet('syn8304095')@filePath,header=T)
  target.data<-target.data[-which(is.na(target.data$Gene)),]
#  target.data<-loadTargetData(target.file)

  target.genes <- unique(target.data$Gene)
  if(use.rna)
    gene.data <- rna.data

  target.genes <- unique(target.data$Gene)
  fObj <- basicFendR(networkFile=network.file,
    featureData=gene.data,
    sampleOutcomeData=pheno.data,
    phenoFeatureData = target.data,
    targetGenes=target.genes
  )

  ##bracket phenotypes into groups of 10-20 to make it easier to compute?
  all.phenos<-unique(pheno.data$Phenotype)
  require(caret)
  testPhenos<-split(all.phenos,ceiling(seq_along(all.phenos)/(length(all.phenos)/30)))
  # kfold<-sapply(caret::createFolds(all.samps,k=k),function(x) all.samps[x])#sapply(caret::createFolds(all.phenos,k=20),function(x) all.phenos[x])
  res.loop<-sapply(names(testPhenos),function(tp) {
    print(paste('Running cross validation comparison for',tp))
    drugs<-testPhenos[[tp]]
    res<-crossValidationCompare(fObj,
      modelCall='lm',
      modelArgs=list(),
      testPheno=drugs,
      k=3,
      sampleIndependent=TRUE)
    plotModelResults(res,prefix=paste0('NTAP_lmFor',tp,ifelse(use.rna,'rnaData','mutData'),'k3'))
    res
  })
}

sangerPerf<-function(use.rna=FALSE){
  network.file<-'https://github.com/fraenkel-lab/OmicsIntegrator/raw/master/data/iref_mitab_miscore_2013_08_12_interactome.txt'

  ##should we load the data or not? seems like a waste of time at this point
  gene.file<-system.file('SANGER_binaryEventMatrix.tsv',package='fendR')
  gene.data<-loadSampleData(gene.file)

  rna.seq.data<-system.file('SANGER_brainarray_rma_expr.tsv.gz', package='fendR')
  rna.data<-loadSampleData(rna.seq.data)

  pheno.file<-system.file('SANGER_dose_response_AUC.tsv',package='fendR')
  pheno.data<-loadPhenotypeData(pheno.file)

  target.file<-system.file('SANGER_drugGeneTargets.tsv',package='fendR')
  target.data<-loadTargetData(target.file)

  if(use.rna)
    gene.data <- rna.data

  target.genes <- unique(target.data$Gene)
  fObj <- basicFendR(networkFile=network.file,
    featureData=gene.data,
    sampleOutcomeData=pheno.data,
    phenoFeatureData = target.data,
    targetGenes=target.genes
  )

  ##bracket phenotypes into groups of 10-20 to make it easier to compute?
  all.phenos<-unique(pheno.data$Phenotype)
  require(caret)
  testPhenos<-split(all.phenos,ceiling(seq_along(all.phenos)/(length(all.phenos)/15)))#
  #sapply(caret::createFolds(all.phenos,k=20),function(x) all.phenos[x])
  res.loop<-sapply(names(testPhenos),function(tp) {
    print(paste('Running cross validation comparison for',tp))
    drugs<-testPhenos[[tp]]
    res<-crossValidationCompare(fObj,
      modelCall='lm',
      modelArgs=list(),
      testPheno=drugs,
      k=3,
      sampleIndependent=TRUE)
    plotModelResults(res,prefix=paste0('Sanger_lmFor',tp,ifelse(use.rna,'rnaData','mutData'),'k10'))
    res
  })
}

cclePerf<-function(use.rna=FALSE){


  network.file<-'https://github.com/fraenkel-lab/OmicsIntegrator/raw/master/data/iref_mitab_miscore_2013_08_12_interactome.txt'

  ##should we load the data or not? seems like a waste of time at this point
  gene.file<-system.file('CCLE_binary_mutation_matrix_ucscGenesFromCBioPortal.tsv',package='fendR')
  gene.data<-loadSampleData(gene.file)

  rna.seq.data<-system.file('CCLE_medianZscore_rnaSeq_ucscGenesFromCbioPortal.tsv.gz', package='fendR')
  rna.data<-loadSampleData(rna.seq.data)

  pheno.file<-system.file('CTRP_v20_AUC_vales_by_drug.tsv',package='fendR')
  pheno.data<-loadPhenotypeData(pheno.file)

  target.file<-system.file('CTRP_v20_drug_target_vals.tsv',package='fendR')
  target.data<-loadTargetData(target.file)

  target.genes <- unique(target.data$Gene)

  if(use.rna)
    gene.data <- rna.data

  fObj <- basicFendR(networkFile=network.file,
    featureData=gene.data,
    sampleOutcomeData=pheno.data,
    phenoFeatureData = target.data,
    targetGenes=target.genes
  )
  fObj <- loadNetwork(fObj) ###only need to load graph once


  ##bracket phenotypes into groups of 10-20 to make it easier to compute?
  all.phenos<-unique(pheno.data$Phenotype)
  require(caret)
  testPhenos<-split(all.phenos,ceiling(seq_along(all.phenos)/(length(all.phenos)/15)))#

  #testPhenos<-sapply(caret::createFolds(all.phenos,k=20),function(x) all.phenos[x])
  res.loop<-sapply(names(testPhenos),function(tp) {
    print(paste('Running cross validation comparison for',tp))
    drugs<-testPhenos[[tp]]
    res<-crossValidationCompare(fObj,
      modelCall='lm',
      modelArgs=list(),
      testPheno=drugs,
      k=10,
      sampleIndependent=TRUE)
   plotModelResults(res,prefix=paste0('CCLE_lmFor',tp,ifelse(use.rna,'rnaData','mutData'),'k10'))
   res
  })

}

sangerPerf(use.rna=TRUE)
#cclePerf(use.rna=TRUE)
#ntapPerf(use.rna=TRUE)
