###build predictive model of dataset using just regression/rf techniques
library(fendR)


##goal is to try
ntapPerf<-function(){

}

sangerPerf<-function(){
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
  testPhenos<-sapply(caret::createFolds(all.phenos,k=20),function(x) all.phenos[x])
  res.loop<-sapply(names(testPhenos),function(tp) {
    print(paste('Running cross validation comparison for',tp))
    drugs<-testPhenos[[tp]]
    res<-crossValidationCompare(fObj,
      modelCall='lm',
      modelArgs=list(),
      testPheno=drugs,
      k=10,
      sampleIndependent=TRUE)
    plotModelResults(res,prefix=paste0('Sanger_lmFor',tp,'k10'))
    res
  })
}

cclePerf<-function(){


  network.file<-'https://github.com/fraenkel-lab/OmicsIntegrator/raw/master/data/iref_mitab_miscore_2013_08_12_interactome.txt'

  ##should we load the data or not? seems like a waste of time at this point
  gene.file<-system.file('CCLE_binary_mutation_matrix_ucscGenesFromCBioPortal.tsv',package='fendR')
  gene.data<-loadSampleData(gene.file)

  rna.seq.data<-system.file('CCLE_medianZscore_rnaSeq_ucscGenesFromCbioPortal.tsv', package='fendR')
  rna.data<-loadSampleData(rna.seq.data)

  pheno.file<-system.file('CTRP_v20_AUC_vales_by_drug.tsv',package='fendR')
  pheno.data<-loadPhenotypeData(pheno.file)

  target.file<-system.file('CTRP_v20_drug_target_vals.tsv',package='fendR')
  target.data<-loadTargetData(target.file)

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
  testPhenos<-sapply(caret::createFolds(all.phenos,k=20),function(x) all.phenos[x])
  res.loop<-sapply(names(testPhenos),function(tp) {
    print(paste('Running cross validation comparison for',tp))
    drugs<-testPhenos[[tp]]
    res<-crossValidationCompare(fObj,
      modelCall='lm',
      modelArgs=list(),
      testPheno=drugs,
      k=10,
      sampleIndependent=TRUE)
   plotModelResults(res,prefix=paste0('CCLE_lmFor',tp,'k10'))
   res
  })

}

sangerPerf()
