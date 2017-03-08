## Run the n3 fendR class against CCLE data

suppressPackageStartupMessages(library(fendR))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(synapseClient))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(gplots))


suppressPackageStartupMessages(library("parallel"))

num.cores <- detectCores()
num.processes <- 1
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  num.processes <- num.cores - 1
  cat(paste("Registering ", num.processes, " cores.\n", sep=""))
  registerDoMC(cores=num.processes)
}


synapseLogin()

##  fit elastic net for a given alpha.
##  NB: foldid argument guarantees the same folds are used across invocations of this function.
##      this allows us to use the same folds to evaluate different alphas.
##  Further: cv.glmnet randomly selects folds, so setting the seed is necessary to ensure comparability
##           across alpha values.
##  And, like, one more:  glmnet's standardize probably doesn't mean what you think it does.  You should standardize x
##                        explicitly if you need to (and be certain to do so before predicting as well).
##                        In particular, ?glmnet reports that 'The coefficients are always returned on the original scale.'
fit.elastic.net.alpha <- function(alpha, x, y, foldid, nfolds, family = "gaussian", type.measure = "mse", seed = 1234) {
    set.seed(seed)
    N <- nrow(x)
    cv <- cv.glmnet(x, y, family = family, type.measure = type.measure, foldid = foldid, nfolds = nfolds, alpha = alpha, standardize = FALSE)    
    cv
}

## Read in CCLE mutation and drug response data.
gene.file <- system.file('CCLE_binary_mutation_matrix_ucscGenesFromCBioPortal.tsv', package='fendR')
gene.data <- loadSampleData(gene.file)
pheno.file <- system.file('CTRP_v20_AUC_vales_by_drug.tsv',package='fendR')
pheno.data <- loadPhenotypeData(pheno.file)

## Read in drug targets.
target.file <- system.file('CTRP_v20_drug_target_vals.tsv', package='fendR')
target.data <- loadTargetData(target.file)

## Download Yuanfang Guan's network (translated from mouse MGI ids to human Hugo symbols by
## aggregrating multiple columns to one using their mean)
network.file <- getFileLocation(synGet("syn8265096"))

target.genes <- unique(target.data$Gene)

testDrugs=c("selumetinib","sorafenib","vorinostat")
testDrugs <- NA

## Create the n3 fendR object
fObj <- n3FendR(network = network.file,
  featureData = gene.data,
  sampleOutcomeData = pheno.data,
  phenoFeatureData = target.data,
  target.genes = target.genes,
  testDrugs = testDrugs
)

## Create the engineered features (nearest neighbor degree = 0)
cat("Enginerring NN = 0 features\n")
fObj.nn0 <- createNewFeaturesFromNetwork(fObj, num.degrees = 0)

## Create the engineered features (nearest neighbor degree = 1)
cat("Enginerring NN = 1 features\n")
fObj.nn1 <- createNewFeaturesFromNetwork(fObj, num.degrees = 1)

## Get the engineered feature matrix, which includes the response in column Response.
cat("Get engineered NN = 0 sparse matrix\n")
m.eng.nn0 <- engineeredSparseResponseMatrix(fObj.nn0)

cat("Get engineered NN = 1 sparse matrix\n")
m.eng.nn1 <- engineeredSparseResponseMatrix(fObj.nn1)

## Get the original feature matrix, but limited to genes included in the engineered features. 
cat("Extracting original features from reduced gene space\n")
m.orig.reduced <- originalResponseMatrix(fObj.nn0, limit.to.engineered.genes = TRUE)

m.orig.sparse.reduced <- originalSparseResponseMatrix(fObj.nn0, limit.to.engineered.genes = TRUE)

## Extract the sample names from the feature matrix
m.eng.samples <- m.eng.nn0$response.mat$Sample
unique.samples <- unique(m.eng.samples)

## Get the original feature matrix.
fit.original.features <- FALSE
if(fit.original.features) {
  cat("Extracting original features\n")
  m.orig <- originalResponseMatrix(fObj.nn0)
}

## Tune elastic net alpha parameter and return the best model
tune.elastic.net <- function(feature.mat, response.vec, pdf.prefix = "elastic-tune") {
  flag <- !is.na(response.vec)
  mat <- feature.mat[flag,]
  resp <- response.vec[flag]
  alphas <- seq(from = 0, to = 1, by = 0.05)
  nfolds <- 5
  N <- length(resp)
  foldid <- sample(rep(seq(nfolds), length = N))
  models <- llply(alphas, .parallel = TRUE, 
                  .fun = function(alpha) {
                           cat(paste0("Fitting elastic net for alpha = ", alpha, "\n"))
                           fit.elastic.net.alpha(alpha = alpha, x = mat, y = resp, foldid = foldid, nfolds = nfolds)
                         })

  ## Generate a heatmap of mean cross-validation error (cvm) as a function of alpha and lambda _index_ (not lambda)
  pdf(paste0(pdf.prefix, ".pdf"), onefile=TRUE)
  ## Get the cross-validated errors
  ncols <- 2
  nrows <- max(2, ceiling(length(alphas)/ncols))
  nrows <- 2
  nrows <- 4
  par(mfrow=c(nrows, ncols))
  for(i in 1:length(alphas)) {
    plot(models[[i]], main = paste0("alpha = ", alphas[i]))
  }
  ##cv.errs <- as.data.frame(lapply(models, function(x) { return(x$cvm) }))
  ##colnames(cv.errs) <- 1:ncol(cv.errs)
  ##heatmap.2(as.matrix(cv.errs), Rowv = F, Colv = F, scale = "none", trace = "none", dendrogram = "none", labCol = alphas, xlab = "alpha", ylab = "lambda index")
  
  ## Find the max AUC based on lambda.1se (NB: this is returned in cvm)
  best.cv.errs <- unlist(lapply(models, function(x) x$cvm[which(x$lambda == x$lambda.1se)]))
  plot(alphas, best.cv.errs, main="lambda.1se")

  best.cv.errs <- unlist(lapply(models, function(x) x$cvm[which(x$lambda == x$lambda.min)]))
  plot(alphas, best.cv.errs, main="lambda.min")
  d <- dev.off()
    
  ## Return the model corresponding to the best.alpha (i.e., that which gives the best error)
  best.alpha <- alphas[which(best.cv.errs == min(best.cv.errs))]
  cat(paste0("Selecting alpha = ", best.alpha, "\n"))
  best.model <- models[[which(alphas == best.alpha)]]
  return(list(model = best.model, alpha = best.alpha))
}

cat("Saving image\n")
save.image(".Rdata")

## Randomly partition samples into training and test sets
set.seed(1234)
training.fraction <- 0.7
training.samples <- sample(unique.samples, replace = FALSE, size = ceiling(training.fraction * length(unique.samples)))
test.samples <- unique.samples[!(unique.samples %in% training.samples)]

## Train glmnet on the original feature matrix
cat("Training glmnet on original features\n")
orig.training.flag <- (m.orig.sparse.reduced$response.mat$Sample %in% training.samples) & (!is.na(m.orig.sparse.reduced$response.mat$Response))
orig.training.mat <- m.orig.sparse.reduced$feature.mat[orig.training.flag, ]
orig.training.response <- m.orig.sparse.reduced$response.mat$Response[orig.training.flag]
orig.best.model <- tune.elastic.net(orig.training.mat, orig.training.response, pdf.prefix = "orig-elastic-tune")

save.image(".Rdata")

## Train glmnet on the engineered feature matrix
cat("Training glmnet on engineered NN = 0 features\n")
eng.nn0.training.flag <- (m.eng.nn0$response.mat$Sample %in% training.samples) & (!is.na(m.eng.nn0$response.mat$Response))
eng.nn0.training.mat <- m.eng.nn0$feature.mat[eng.nn0.training.flag, ]
eng.nn0.training.response <- m.eng.nn0$response.mat$Response[eng.nn0.training.flag]
eng.nn0.best.model <- tune.elastic.net(eng.nn0.training.mat, eng.nn0.training.response, pdf.prefix = "eng-nn0-elastic-tune")

save.image(".Rdata")

## Train glmnet on the engineered feature matrix
cat("Training glmnet on engineered NN = 1 features\n")
eng.nn1.training.flag <- (m.eng.nn1$response.mat$Sample %in% training.samples) & (!is.na(m.eng.nn1$response.mat$Response))
eng.nn1.training.mat <- m.eng.nn1$feature.mat[eng.nn1.training.flag, ]
eng.nn1.training.response <- m.eng.nn1$response.mat$Response[eng.nn1.training.flag]
eng.nn1.best.model <- tune.elastic.net(eng.nn1.training.mat, eng.nn1.training.response, pdf.prefix = "eng-nn1-elastic-tune")

save.image(".Rdata")

## Perform LOO validation
cat("Performing LOO on engineered NN = 0 features\n")
eng.nn0.mses <- llply(unique.samples, .parallel = TRUE,
                      .fun = function(sample) {
                               loo.flag <- !(m.eng.nn0$response.mat$Sample %in% training.samples) & (!is.na(m.eng.nn0$response.mat$Response)) & (m.eng.nn0$response.mat$Sample == sample)
                               train.flag <- !(m.eng.nn0$response.mat$Sample %in% training.samples) & (!is.na(m.eng.nn0$response.mat$Response)) & (m.eng.nn0$response.mat$Sample != sample)
                               loo.mat <- m.eng.nn0$feature.mat[loo.flag,,drop=FALSE ]
                               loo.resp <- m.eng.nn0$response.mat$Response[loo.flag]
                               if(nrow(loo.mat) == 0) { return(-1) }
                               train.mat <- m.eng.nn0$feature.mat[train.flag,,drop=FALSE ]
                               train.resp <- m.eng.nn0$response.mat$Response[train.flag]
                               alpha <- eng.nn0.best.model$alpha

                               glm.fit <- cv.glmnet(x = train.mat, y = train.resp, family = "gaussian", alpha = alpha, standardize = FALSE)
                               scores <- predict(glm.fit, newx=loo.mat, s="lambda.min", type="response")
                               diff <- scores - loo.resp
                               diff <- diff[!is.na(diff)]
                               mse <- sum(diff * diff)/length(diff)
                               mse 
                             })
eng.nn0.mses <- unlist(eng.nn0.mses)
eng.nn0.mses <- eng.nn0.mses[eng.nn0.mses != -1]
save.image(".Rdata")

cat("Performing LOO on engineered NN = 1 features\n")
eng.nn1.mses <- llply(unique.samples, .parallel = TRUE,
                      .fun = function(sample) {
                               loo.flag <- !(m.eng.nn1$response.mat$Sample %in% training.samples) & (!is.na(m.eng.nn1$response.mat$Response)) & (m.eng.nn1$response.mat$Sample == sample)
                               train.flag <- !(m.eng.nn1$response.mat$Sample %in% training.samples) & (!is.na(m.eng.nn1$response.mat$Response)) & (m.eng.nn1$response.mat$Sample != sample)
                               loo.mat <- m.eng.nn1$feature.mat[loo.flag,,drop=FALSE ]
                               loo.resp <- m.eng.nn1$response.mat$Response[loo.flag]
                               if(nrow(loo.mat) == 0) { return(-1) }
                               train.mat <- m.eng.nn1$feature.mat[train.flag,,drop=FALSE ]
                               train.resp <- m.eng.nn1$response.mat$Response[train.flag]
                               alpha <- eng.nn1.best.model$alpha

                               glm.fit <- cv.glmnet(x = train.mat, y = train.resp, family = "gaussian", alpha = alpha, standardize = FALSE)
                               scores <- predict(glm.fit, newx=loo.mat, s="lambda.min", type="response")
                               diff <- scores - loo.resp
                               diff <- diff[!is.na(diff)]
                               mse <- sum(diff * diff)/length(diff)
                               mse 
                             })
eng.nn1.mses <- unlist(eng.nn1.mses)
eng.nn1.mses <- eng.nn1.mses[eng.nn1.mses != -1]
save.image(".Rdata")

cat("Performing LOO on original features\n")
orig.mses <- llply(unique.samples, .parallel = TRUE,
                      .fun = function(sample) {
                               loo.flag <- !(m.orig.sparse.reduced$response.mat$Sample %in% training.samples) & (!is.na(m.orig.sparse.reduced$response.mat$Response)) & (m.orig.sparse.reduced$response.mat$Sample == sample)
                               train.flag <- !(m.orig.sparse.reduced$response.mat$Sample %in% training.samples) & (!is.na(m.orig.sparse.reduced$response.mat$Response)) & (m.orig.sparse.reduced$response.mat$Sample != sample)
                               loo.mat <- m.orig.sparse.reduced$feature.mat[loo.flag,,drop=FALSE ]
                               loo.resp <- m.orig.sparse.reduced$response.mat$Response[loo.flag]
                               if(nrow(loo.mat) == 0) { return(-1) }
                               train.mat <- m.orig.sparse.reduced$feature.mat[train.flag,,drop=FALSE ]
                               train.resp <- m.orig.sparse.reduced$response.mat$Response[train.flag]
                               alpha <- orig.best.model$alpha

                               glm.fit <- cv.glmnet(x = train.mat, y = train.resp, family = "gaussian", alpha = alpha, standardize = FALSE)
                               scores <- predict(glm.fit, newx=loo.mat, s="lambda.min", type="response")
                               diff <- scores - loo.resp
                               diff <- diff[!is.na(diff)]
                               mse <- sum(diff * diff)/length(diff)
                               mse 
                             })
orig.mses <- unlist(orig.mses)
orig.mses <- orig.mses[orig.mses != -1]
save.image(".Rdata")

mse.df <- data.frame(MSE = eng.nn0.mses, Method = rep("Eng-NN0", length(eng.nn0.mses)))
mse.df <- rbind(mse.df, data.frame(MSE = eng.nn1.mses, Method = rep("Eng-NN1", length(eng.nn1.mses))))
mse.df <- rbind(mse.df, data.frame(MSE = orig.mses, Method = rep("Orig", length(orig.mses))))

library(ggplot2)
library(ggbeeswarm)
pdf("mse.pdf")
g <- ggplot(data = mse.df, aes(x = Method, y = MSE))
g <- g + geom_boxplot()
g <- g + geom_beeswarm()
print(g)
dev.off()

stop("stop")

## Perform bootstrap sampling

## y = TRUE --> return response
cat("Fitting engineered features\n")
## Drop the Sample and Phenotype columns, leaving only the gene/data columns and the Response column.
m.eng.sub <- subset(m.eng, select = !(colnames(m.eng) %in% c("Sample", "Phenotype")))
engineered.fit <- lm(Response ~ ., m.eng.sub, y = TRUE)

if(fit.original.features) {
  cat("Fitting original features\n")
  m.orig.sub <- subset(m.orig, select = !(colnames(m.orig) %in% c("Sample", "Phenotype")))
  orig.fit <- lm(Response ~ ., m.orig.sub, y = TRUE)
}

cat("Fitting original features in reduced gene space\n")
m.orig.reduced.sub <- subset(m.orig.reduced, select = !(colnames(m.orig.reduced) %in% c("Sample", "Phenotype")))
orig.reduced.fit <- lm(Response ~ ., m.orig.reduced.sub, y = TRUE)

save.image(".Rdata.fit")

stop("stop")

##performs loo cross validation by sample - should we move to other file?
looCrossValidation<-function(obj,testDrugs){
  #iterate over all cell lines
  all.samps<-intersect(obj$sampleOutcomeData$Sample,obj$featureData$Sample)

  vals<-sapply(all.samps,function(x){
    print(paste('Removing sample',x,'to evaluate'))
    #subset out that data and re-assign original object
    test.data<-subset(obj$sampleOutcomeData,Sample==x)
    orig.test.features<-subset(obj$featureData,Sample==x)
    aug.test.features<-subset(obj$remappedFeatures,Sample==x)

    ##now remove from object
    newObj<-removeSampleFromObject(obj,x)

    #create baseline model
    baselineModelObj<-buildModelFromOriginalFeatures(newObj,testDrugs)

    #get baseline predictions
    baselinePreds<-scoreDataFromModel(baselineModelObj,orig.test.features,test.data)

    #create new features
    augmentedObj<-createNewFeaturesFromNetwork(newObj,testDrugs)

    #create new model - this will replace the model list object in the class
    augmentedObj<-buildModelFromEngineeredFeatures(augmentedObj,testDrugs)

    updatedPreds<-scoreDataFromModel(augmentedObj,aug.test.features,test.data)
    df<-data.frame(select(baselinePreds,originalPred=Prediction,Actual),select(updatedPreds,augmentedPred=Prediction))
    df$Drug<-rownames(df)
    df$Sample<-rep(x,nrow(df))
    df
  })
  vals

}

##we can add some generic fendR methods as well, such as plotting, statistics, loo, etc.
