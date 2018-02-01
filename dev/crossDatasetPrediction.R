##
## predict from one dataset to another
##
library(fendR)
library(synapseClient)
library(dplyr)

all.datasets=c('CCLE','Sanger','NTAP')

#build fendR from dataset
buildFendRFromDataset<-function(ds.name,use.rna=TRUE,responseNorm=NA,geneNorm=NA){
  if(!ds.name%in%all.datasets){
    print(paste(ds.name,'not available, select from:',paste(all.datasets,collapse=',')))
    return(NULL)
  }

  network.file<-'https://github.com/fraenkel-lab/OmicsIntegrator/raw/master/data/iref_mitab_miscore_2013_08_12_interactome.txt'

  if(ds.name=='Sanger'){
    gene.file<-system.file('SANGER_binaryEventMatrix.tsv',package='fendR')
    gene.data<-loadSampleData(gene.file)

    rna.seq.data<-system.file('SANGER_brainarray_rma_expr.tsv.gz', package='fendR')
    rna.data<-loadSampleData(rna.seq.data)

    pheno.file<-system.file('SANGER_dose_response_AUC.tsv',package='fendR')
    pheno.data<-loadPhenotypeData(pheno.file)

    target.file<-system.file('SANGER_drugGeneTargets.tsv',package='fendR')
    target.data<-loadTargetData(target.file)
  }else if(ds.name=='CCLE'){
    gene.file<-system.file('CCLE_binary_mutation_matrix_ucscGenesFromCBioPortal.tsv',package='fendR')
    gene.data<-loadSampleData(gene.file)

    rna.seq.data<-system.file('CCLE_medianZscore_rnaSeq_ucscGenesFromCbioPortal.tsv.gz', package='fendR')
    rna.data<-loadSampleData(rna.seq.data)

    pheno.file<-system.file('CTRP_v20_AUC_vales_by_drug.tsv',package='fendR')
    pheno.data<-loadPhenotypeData(pheno.file)

    target.file<-system.file('CTRP_v20_drug_target_vals.tsv',package='fendR')
    target.data<-loadTargetData(target.file)
  }else if(ds.name=='NTAP'){
    require(synapseClient)
    synapseLogin()
    gene.data<-read.table(synGet('syn8304629')@filePath,header=T)
    #gene.data<-loadSampleData(gene.file)

    # rna.seq.data<-
    rna.data<-read.table(synGet('syn8304627')@filePath,header=T)

    pheno.data<-read.table(synGet('syn8304620')@filePath,header=T)
    pheno.data<-pheno.data[-which(is.na(pheno.data$Phenotype)),]

    #  pheno.data<-loadPhenotypeData(pheno.file)

    target.data<-read.table(synGet('syn8304095')@filePath,header=T)
    target.data<-target.data[-which(is.na(target.data$Gene)),]
  }
  if(use.rna)
    gene.data <- rna.data

  target.genes<-unique(target.data$Gene)
  fObj <- basicFendR(networkFile=network.file,
    featureData=gene.data,
    sampleOutcomeData=pheno.data,
    phenoFeatureData = target.data,
    targetGenes=target.genes,
    responseNorm=responseNorm,
    geneNorm=geneNorm
  )


  return(fObj)
}

#add in quick zscore function
zscore <-function(x){
  (x-mean(x,na.rm=T))/sd(x,na.rm=T)
}



