
## This is an implementation of the nearest-neighbor network (n3) fendR class.
## It is based on the methods used by Yuanfang Guan in DREAM challenges.

######################################################################
# Create the n3FendR class
#
# This is used to represent the nearest-neighbor network implementation of the fendR framework.

#' An S3 class to represent the nearest-neighbor network implementation of the fendR predictive
#' network algorithm.
#' @param network the file name of a feather object representing a dense matrix
#' @param featureData a data.frame that contains rows representing genes and columns representing samples
#' @param sampleOutcomeData a data.frame representing at least one column of phenotype and rows representing samples
#' @param phenoFeatureData a data.frame where rows represent genes and columns represent a relationship between phenotype and gene
#' @param target.genes a vector of target gene names (a subset of those in featureData) that will be used to sparsify the network
#' @param testDrugs a vector of drug names to be used to subset the data (e.g., for debugging/efficiency purposes)
#' @inheritParams fendR
#' @export
#' @return n3FendR object
n3FendR <- function(network, featureData, phenoFeatureData, sampleOutcomeData, target.genes, testDrugs = NA) {

  phenos <- intersect(sampleOutcomeData$Phenotype, phenoFeatureData$Phenotype)
  if(!is.na(testDrugs) && any(testDrugs %in% phenos)) {
    cat(paste0("Reducing scope to only focus on ", paste(testDrugs, collapse=', '), " drugs\n"))
    phenos <- phenos[phenos %in% testDrugs]
  }
  sampleOutcomeData <- subset(sampleOutcomeData, Phenotype %in% phenos)
  phenoFeatureData <- subset(phenoFeatureData, Phenotype %in% phenos)

  me <-fendR(network, featureData, phenoFeatureData, sampleOutcomeData)
  cls.me <- class(me)

  ## Get the genes in the network
  all.network.genes <- names(feather_metadata(network)$types)

##  full.gene.set <- unique(intersect(featureData$Gene, intersect(phenoFeatureData$Gene, all.network.genes)))
  full.gene.set <- unique(intersect(featureData$Gene, all.network.genes))
  num.genes.in.full.feature.space <- length(full.gene.set)
  reduced.gene.set <- full.gene.set
  if(!is.na(target.genes) && !is.null(target.genes)) {
    if(!any(target.genes %in% full.gene.set)) {
      stop("None of target genes are in the featureData")
    }
    reduced.gene.set <- full.gene.set[full.gene.set %in% target.genes]
    cat(paste0("Reducing feature space from ", num.genes.in.full.feature.space, " to ", length(reduced.gene.set), "\n"))
  } else {
    cat(paste0("Using all ", num.genes.in.full.feature.space, " genes in featureData.\n"))
    cat("Not sparsifying data or network\n")
  }

  me <- append(me, list(target.genes = target.genes, reduced.gene.set = reduced.gene.set, all.network.genes = all.network.genes))
  class(me) <- c("n3FendR", cls.me)
  return(me)
}

######################################################################
# Set methods implemented by n3FendR



#' Engineer Features from Network
#' \code{createNewFeaturesFromNetwork} takes the phenotype scores and updates the gene scores based on phenotype
#' @param object That contains a data frame and network
#' @param num.degrees Consider nearest neighbors separated by num.degrees (only 0 or 1 are implemeneted)
#' @keywords
#' @export
#' @return list of gene features for each phenotype/drug response
createNewFeaturesFromNetwork.n3FendR <- function(object, testDrugs = NA, num.degrees = 1) {
   # suppressPackageStartupMessages(library(plyr))
  #  suppressPackageStartupMessages(library(feather))

    sampleOutcomeData <- object$sampleOutcomeData
    phenoFeatureData <- object$phenoFeatureData
    featureData <- object$featureData

    if(!(num.degrees %in% c(0,1))) {
      stop(paste0("Must specify num.degrees = 0 or 1\n"))
    }

    ## Sparsify the data (featureData, phenoFeatureData, and network) by only considering a subset of
    ## curated target.genes -- these are likely drug targets and/or "cancer genes," etc.
    ## NB: by doing this up front/here we do not propagate mutations in non target genes
    ## to the target genes.

    featureData <- subset(featureData, Gene %in% object$reduced.gene.set)
    ## phenoFeatureData <- subset(phenoFeatureData, Gene %in% object$reduced.gene.set)
    ## Subset to phenotypes having both feature data and outcome data
    phenos <- intersect(sampleOutcomeData$Phenotype, phenoFeatureData$Phenotype)
    all.phenos <- union(sampleOutcomeData$Phenotype, phenoFeatureData$Phenotype)
    cat(paste0("Found ", length(phenos), " phenotypes that have feature data and outcome data out of ", length(all.phenos), "\n"))

    if(!is.na(testDrugs) && any(testDrugs %in% phenos)) {
      stop("This should be done in the constructor.")
      cat(paste0("Reducing scope to only focus on ", paste(testDrugs, collapse=', '), " drugs\n"))
      phenos <- phenos[phenos %in% testDrugs]
    }

    ## Convert featureData in tidy format to a matrix m (with rows = samples, columns = genes)
    ## having entries m_sg for sample s and gene g
    m <- tidyr::spread(featureData, key="Sample", value="Value")
    rownames(m) <- m$Gene
    m <- m[, colnames(m) != "Gene"]
    m <- t(m)
    m <- as.matrix(m[, object$reduced.gene.set])

    ## If any of the anchors do not show up in the mutation matrix (e.g., this
    ## might happen if mutations are discovered using targeted sequencing)
    ## then add them
    all.anchors <- unique(unlist(ldply(phenos, .parallel = TRUE, .fun = function(pheno) {
      
      ## For each gene target of that phenotype/drug
      anchors <- as.character(subset(phenoFeatureData, Phenotype == pheno)$Gene)
      anchors
    })))
    all.anchors <- all.anchors[all.anchors %in% object$all.network.genes]
    
    if(num.degrees == 1 ) {
      for(anchor in all.anchors) {
        if(!(anchor %in% colnames(m))) {
          m <- cbind(m, rep(0, nrow(m)))
          colnames(m)[ncol(m)] <- anchor
        }
      }
    }
    
    features <- colnames(m)
    
    ## For each anchor (e.g., drug target) a of phenotype (e.g., drug) p, define a feature matrix m^a 
    ## by propagating mutations in the original matrix m to nearest neighbors.  The network induced by anchor a is
    ## represented by the edge list e^a_gg', where e^a_gg' > 0 only if g and/or g' = a and/or g = g'
    ## The feature matrix m^d _for the phenotype_ is then defined over all anchors a for that phenotype
    ## e.g., by assigning
    ## m^d_sg = max_a m^a_sg, or
    ## m^d_sg = min_a m^a_sg, or
    ## m^d_sg = mean_a m^a_sg
    ##
    ## If g = anchor a, then
    ## m^a_sg = \sum_g' e^a_gg' m_sg'  (i.e., matrix x vector: m^a_,g = m e^a)
    ## If g != anchor a, then
    ## m^a_sg = delta_a,g' e_gg' m_sg' + delta_g,g' e_gg' m_sg'
    ## i.e., for a non-anchor g, the only contributions are from g and from the anchor a

    ## For each phenotype/drug
    pheno.features <- plyr::ldply(phenos, .parallel = TRUE, .fun = function(pheno) {

      ## For each gene target of that phenotype/drug
      anchors <- as.character(subset(phenoFeatureData, Phenotype == pheno)$Gene)
      anchors <- anchors[anchors %in% object$all.network.genes]
      
      cat(paste0("Performing nearest neighbor proportion for phenotype ", pheno, " with target(s): ", paste(anchors, collapse=","), "\n"))

      ## For each gene target/anchor of that phenotype/drug, propagate mutation to nearest neighbors
      anchor.responses <- lapply(anchors, function(anchor) {
        ## Read in column of feather network
        ## NB: we have ensured above this gene is in the network, so this read will not fail.
        anchor.weights <- as.data.frame(read_feather(object$network, columns=c(anchor)))[,1]
        names(anchor.weights) <- object$all.network.genes
        anchor.weights <- as.vector(anchor.weights[features])
        names(anchor.weights) <- features
        m.a <- matrix(data = 0, nrow = nrow(m), ncol = ncol(m))
        rownames(m.a) <- rownames(m)
        colnames(m.a) <- colnames(m)

        if(num.degrees == 0 ) {
          ## For zeroth order, simply weight each mutation by its distance to the anchor
          ## i.e., m^a_sg = e_ga m_sg  (no sum)
          for(g in names(anchor.weights)) {
            m.a[, g] <- m[, g] * anchor.weights[g]
          }
        } else if(num.degrees == 1) {
          ## m^a_sg = \sum_g' e^a_gg' m_sg'  (i.e., matrix x vector: m^a_,g = m e^a)
          ## Handle the anchor
          m.a[, anchor] <- m %*% anchor.weights
          ## For each non-anchor neighbor
          ## If g != anchor a, then
          ## m^a_sg = delta_a,g' e_gg' m_sg' + delta_g,g' e_gg' m_sg'
          ## i.e., for a non-anchor g, the only contributions are from g and from the anchor a
          non.anchors <- names(anchor.weights)[names(anchor.weights) != anchor]
          for(nn in non.anchors) {
            ## Below assumes that e_gg = 1
            ## m.a[, nn] <- ( m[, anchor] * anchor.weights[nn] ) + ( m[, nn] * 1 )
            ## Don't assume that e_gg = 1.  Instead, we need to read the column for this
            ## nearest neighbor nn from the feather object
            nn.weights <- as.data.frame(read_feather(object$network, columns=c(nn)))[,1]
            names(nn.weights) <- object$all.network.genes
            nn.weights <- as.vector(nn.weights[features])
            names(nn.weights) <- features
            m.a[, nn] <- ( m[, nn] * nn.weights[nn] ) + ( m[, anchor] * nn.weights[anchor] )
          }
        }
        m.a
      })

      ## anchor.response is a list of matrices m.a over anchors a.
      ## Define a single aggregate matrix m.p over the phenotype p.
      ## Obvious ways to do this include:
      ## 1. Taking elementwise means across the m.a matrices
      ## m.d <- apply(simplify2array(anchor.responses), c(1,2), mean)
      ## 2. Taking the elementwise min across the m.a matrices
      ## m.d <- apply(simplify2array(anchor.responses), c(1,2), min)
      ## 3. Taking the elementwise max across the m.a matrices
      m.d <- apply(simplify2array(anchor.responses), c(1,2), max)

      ## Return the results in matrix, rather than tidy, format
      ## Append the phenotype name to the sample
      m.d <- as.data.frame(m.d)
      m.d$Sample <- rownames(m.d)
      m.d$Phenotype <- pheno
      rownames(m.d) <- NULL
      gc()
      return(m.d)

      ## Convert to matrix to tidy format
      df <- as.data.frame(t(m.d))
      df$Gene <- rownames(df)
      df <- tidyr::gather(df, key="Sample", value="Value", -Gene)
      df$Phenotype <- pheno
      gc()
      df
    })

    ## Store the remapped features in the object
    object$remappedFeatures <- pheno.features

    return(object)
}

#' Return feature matrix, consisting of the original features, and the associated responses in a separate data.frame.
#' @description Creates a matrix suitable for modeling with the formula Response ~ from the original features.
#' @keywords model fendR
#' @keywords limit.to.engineered.genes Boolen indicating whether the returned features should be limited to those used in the engineered features (which may be restricted beyond the target.genes passed to the constructor by those in the network)
#' @export
#' @import plyr Matrix
#' @return A list with slots "feature.mat" (the feature matrix) and "response.mat" (a data.frame with columns "Sample", "Phenotype", and "Response").  The rows of feature.mat and response.mat are in 1-to-1 correspondence.
originalSparseResponseMatrix.n3FendR <- function(object, phenotype=c(), limit.to.engineered.genes = FALSE){

  suppressPackageStartupMessages(library(plyr))
  suppressPackageStartupMessages(library(Matrix))

  if(length(phenotype)>0) {
    stop("Restrict drugs/phenotypes in the constructor; not here")
  }

  ## NB: for now, we will not subset the data to the target genes.
  ## We might want to do that for sake of comparison.

  ## Convert featureData in tidy format to a matrix m (with rows = samples, columns = genes)
  ## having entries m_sg for sample s and gene g
  featureData <- object$featureData
  m <- tidyr::spread(featureData, key="Sample", value="Value")
  rownames(m) <- m$Gene
  m <- m[, colnames(m) != "Gene"]
  m <- t(m)

  if(limit.to.engineered.genes) {
    m <- (m[, object$reduced.gene.set])
  }
  m <- as.matrix(m)

  phenos <- intersect(object$sampleOutcomeData$Phenotype, object$phenoFeatureData$Phenotype)

  num.phenos <- length(phenos)

  kr <- kronecker(Diagonal(num.phenos), Matrix(m, sparse=TRUE))

  response <- expand.grid(rownames(m), phenos)
  colnames(response) <- c("Sample", "Phenotype")
  response <- merge(response, object$sampleOutcomeData, by = c("Sample", "Phenotype"), all.x = TRUE)

  return(list("feature.mat" = kr, "response.mat" = response))
}

#' Return feature matrix, consisting of the original features, and the associated responses in column Response.
#' @description Creates a matrix suitable for modeling with the formula Response ~ from the original features.
#' @keywords model fendR
#' @keywords limit.to.engineered.genes Boolen indicating whether the returned features should be limited to those used in the engineered features (which may be restricted beyond the target.genes passed to the constructor by those in the network)
#' @export
#' @return A matrix suitable for modeling with the formula Response ~ from the original features, whose columns are genes and the Response and whose rows are "phenotype" . "sample".
originalResponseMatrix.n3FendR <- function(object, phenotype=c(), limit.to.engineered.genes = FALSE){

  suppressPackageStartupMessages(library(plyr))
  suppressPackageStartupMessages(library(Matrix))

  if(length(phenotype)>0) {
    stop("Restrict drugs/phenotypes in the constructor; not here")
  }

  ## NB: for now, we will not subset the data to the target genes.
  ## We might want to do that for sake of comparison.

  ## Convert featureData in tidy format to a matrix m (with rows = samples, columns = genes)
  ## having entries m_sg for sample s and gene g
  featureData <- object$featureData
  m <- tidyr::spread(featureData, key="Sample", value="Value")
  rownames(m) <- m$Gene
  m <- m[, colnames(m) != "Gene"]
  m <- t(m)

  if(limit.to.engineered.genes) {
    m <- as.matrix(m[, object$reduced.gene.set])
  }

  phenos <- intersect(object$sampleOutcomeData$Phenotype, object$phenoFeatureData$Phenotype)

  ## Create a feature matrix that simply appends the original feature matrices vertically.

  ## But also append the sample and phenotype.
  pheno.features <- ldply(phenos, .parallel = TRUE, .fun = function(pheno) {
    m.d <- as.data.frame(m)
    m.d$Sample <- rownames(m.d)
    m.d$Phenotype <- pheno

    return(m.d)
  })

  newdf <- merge(pheno.features, object$sampleOutcomeData, by = c("Sample", "Phenotype"), all = FALSE)

  return(newdf)
}

#' Return re-engineered feature matrix, consisting of re-engineered features and the associated responses in column Response.
#' @description Creates a matrix suitable for modeling with the formula Response ~ from the re-engineered features.
#' @keywords model fendR
#' @export
#' @return A matrix suitable for modeling with the formula Response ~ from the re-engineered features, whose columns are genes and the Response and whose rows are "phenotype" . "sample".
engineeredResponseMatrix.n3FendR <- function(object,phenotype=c()){
  if(length(phenotype)>0) {
    stop("Restrict drugs/phenotypes in the constructor; not here")
  }

  newdf <- merge(object$remappedFeatures, object$sampleOutcomeData, by = c("Sample", "Phenotype"), all = FALSE)

  return(newdf)
}

engineeredSparseResponseMatrix.n3FendR <- function(object,phenotype=c()){
  if(length(phenotype)>0) {
    stop("Restrict drugs/phenotypes in the constructor; not here")
  }
  newdf <- merge(object$remappedFeatures, object$sampleOutcomeData, by = c("Sample", "Phenotype"), all = FALSE)

  cols <- c("Sample", "Phenotype", "Response")
  response <- newdf[,cols]
  mat <- newdf[, !(colnames(newdf) %in% cols)]
  mat <- Matrix(as.matrix(mat), sparse=TRUE)

  return(list("feature.mat" = mat, "response.mat" = response))
}
