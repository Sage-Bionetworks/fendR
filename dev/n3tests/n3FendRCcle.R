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

library("downloader")

## Download CCLE oncomap calls
url <- "https://portals.broadinstitute.org/ccle/downloadFile/DefaultSystemRoot/exp_10/ds_23/CCLE_Oncomap3_2012-04-09.maf?downloadff=true&fileId=3000"
dest <- "CCLE_Oncomap3_2012-04-09.maf"
download(url, destfile = dest)
oncomap.tbl <- read.table(dest, sep="\t", header=TRUE, quote="") 

oncomap.tbl$ccl_name <- unlist(lapply(oncomap.tbl$Tumor_Sample_Barcode, 
                                      function(x) {
                                        r <- regexpr(pattern="^([^_]+)", text=x)
                                        substr(x, r[1], r[1] + attr(r, "match.length")[1] - 1)
                                      }))

## Download CCLE targeted sequencing neutral missense variants filtered (called TES-A in Basu)
url <- "https://portals.broadinstitute.org/ccle/downloadFile/DefaultSystemRoot/exp_10/ds_26/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf.gz?downloadff=true&fileId=6873"
dest <- "CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf.gz"
download(url, destfile = dest)
tesa.tbl <- read.table(gzfile(dest), sep="\t", header=TRUE, quote="")

## Download CCLE SNV data
url <- "https://portals.broadinstitute.org/ccle/downloadFile/DefaultSystemRoot/exp_10/ds_20/CCLE_copynumber_byGene_2013-12-03.txt.gz?downloadff=true&fileId=17599"
dest <- "CCLE_copynumber_byGene_2013-12-03.txt.gz"
download(url, destfile = dest)
cnv.tbl <- read.table(gzfile(dest), sep="\t", header=TRUE, quote="")
colnames(cnv.tbl) <- unlist(lapply(colnames(cnv.tbl),
                                   function(x) {
                                     r <- regexpr(pattern="^([^_]+)", text=x)
                                     substr(x, r[1], r[1] + attr(r, "match.length")[1] - 1)
                                   }))

tesa.tbl$ccl_name <- unlist(lapply(tesa.tbl$Tumor_Sample_Barcode, 
                                      function(x) {
                                        r <- regexpr(pattern="^([^_]+)", text=x)
                                        substr(x, r[1], r[1] + attr(r, "match.length")[1] - 1)
                                      }))


sample.info.tbl$ccl_name <- unlist(lapply(sample.info.tbl$CCLE.name, 
                                      function(x) {
                                        r <- regexpr(pattern="^([^_]+)", text=x)
                                        substr(x, r[1], r[1] + attr(r, "match.length")[1] - 1)
                                      }))


url <- "https://portals.broadinstitute.org/ccle/downloadFile/DefaultSystemRoot/exp_10/ds_23/CCLE_Oncomap3_Assays_2012-04-09.csv?downloadff=true&fileId=3001"
dest <- "CCLE_Oncomap3_Assays_2012-04-09.csv"
download(url, destfile = dest)
oncomap.assays.tbl <- read.table(dest, sep=",", header=TRUE, quote="\"") 

url <- "https://portals.broadinstitute.org/ccle/downloadFile/DefaultSystemRoot/exp_10/ds_22/CCLE_sample_info_file_2012-10-18.txt?downloadff=true&fileId=6801"
dest <- "CCLE_sample_info_file_2012-10-18.txt"
download(url, destfile = dest)
sample.info.tbl <- read.table(dest, sep="\t", header=TRUE)

url <- "ftp://anonymous:guest@caftpd.nci.nih.gov/pub/OCG-DCC/CTD2/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset/CTRPv2.0_2015_ctd2_ExpandedDataset.zip"
dest <- "CTRPv2.0_2015_ctd2_ExpandedDataset.zip"
download(url = url, destfile = dest)
unzip(dest)

url <- "ftp://anonymous:guest@caftpd.nci.nih.gov/pub/OCG-DCC/CTD2/Broad/CTRPv1.0_2013_pub_Cell_154_1151/CTRPv1.0_2013_pub_Cell_154_1151.zip"
dest <- "CTRPv1.0_2013_pub_Cell_154_1151.zip"
download(url = url, destfile = dest)
unzip(dest)

## url <- "ftp://anonymous:guest@caftpd.nci.nih.gov/pub/OCG-DCC/CTD2/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset/CTRPv2.0._COLUMNS.xlsx"
## dest <- "CTRPv2.0._COLUMNS.xlsx"

## Load the AUC drug-response results
file <- "v20.data.curves_post_qc.txt"
ctrpv2_aucs <- read.table(file, sep="\t", header=TRUE)

file <- "v10.D3.area_under_conc_curve.txt"
ctrpv1_aucs <- read.table(file, sep="\t", header=TRUE)


plot.auc.vs.annotation <- function(tbl, anno.col, sample.col = "ccl_name", auc.col = "area_under_curve", main = NULL) {
  library(reshape2)
  library(gridExtra)
  tbl <- tbl[!is.na(tbl[,anno.col]),]
  tbl <- tbl[!is.na(tbl[,auc.col]),]
  
  ## Sometimes we have replicate samples -- make them unique
  tbl[, sample.col] <- paste(tbl[, sample.col], 1:nrow(tbl), sep="")
  
  o <- order(tbl[,auc.col], -as.numeric(as.character(tbl[,anno.col])))
##  o <- order(tbl[,auc.col])  
  levels <- unique(tbl[o,sample.col])
  
  flag <- tbl[,auc.col] < 1
  tbl[flag,auc.col] <- 1
  flag <- tbl[,auc.col] > 6
  tbl[flag,auc.col] <- 6
  
  tbl.melt <- melt(tbl[,c(sample.col, auc.col)], id.vars = sample.col)
##  o <- order(tbl.melt$value, decreasing = FALSE)
##  levels <- unique(tbl.melt[o,sample.col])
  tbl.melt[,sample.col] <- factor(tbl.melt[,sample.col], levels = levels)
  p.heat <- ggplot(data = tbl.melt, aes_string(x = "variable", y = sample.col))
  p.heat <- p.heat + geom_tile(aes(fill = value))
  p.heat <- p.heat + scale_fill_gradient2(high = "blue", mid = "white", low = "red", limits = c(1.0, 6.0), midpoint = 3.5, breaks = c(1.0, 3.5, 6.0))
  p.heat <- p.heat + guides(fill = guide_colourbar(title = "AUC"))
  
  g.heat <- ggplotGrob(p.heat)
  g.legend.heat <- gtable::gtable_filter(g.heat, "guide-box")
  heat.label <- gtable::gtable_filter(g.heat, "axis-b")
  g.heat <- gtable::gtable_filter(g.heat, "panel")

  mt.wt <- 0  
  if((length(unique(tbl[,anno.col])) == 2) && (all(tbl[,anno.col] %in% c(0,1)))) {
    mt.wt <- 1
  }
  if(mt.wt == 1) {
    tbl[, anno.col] <- as.numeric(as.character(tbl[, anno.col]))
    flag <- tbl[,anno.col] == 1
    tbl[flag,anno.col] <- "MT"
    tbl[!flag,anno.col] <- "WT"
    tbl[, anno.col] <- factor(tbl[, anno.col])
  }
  tbl.melt <- melt(tbl[,c(sample.col, anno.col)], id.vars = sample.col)
  tbl.melt[,sample.col] <- factor(tbl.melt[,sample.col], levels = levels)
  p.anno <- ggplot(data = tbl.melt, aes_string(x = "variable", y = sample.col))
  p.anno <- p.anno + geom_tile(aes(fill = value))
  if(mt.wt == 1) {
    p.anno <- p.anno + scale_fill_manual(values = c("MT" = "black", "WT" = "white"))
    p.anno <- p.anno + guides(fill = guide_legend(title = toupper(anno.col)))
  } else {
    p.anno <- p.anno + guides(fill = guide_colourbar(title = toupper(anno.col)))
    p.anno <- p.anno + scale_fill_gradient(high = "black", low = "white")    
  }
  
  g.anno <- ggplotGrob(p.anno)
  g.legend.anno <- gtable::gtable_filter(g.anno, "guide-box")
  anno.label <- gtable::gtable_filter(g.anno, "axis-b")
  g.anno <- gtable::gtable_filter(g.anno, "panel")
  grobs <- list(g.heat, g.anno, g.legend.heat, g.legend.anno, heat.label, anno.label)
  layout <- rbind(c(1,2,3),c(1,2,4),c(5,6,NA))
  heights <- c(1,1,0.5)
  widths <- c(1,1,0.5)
  if(!is.null(main)) {
    layout <- rbind(c(7,7,NA), layout)
    heights <- c(0.25, heights)
    p.title <- ggplot()
    p.title <- p.title + ggtitle(main)
    g.title <- ggplotGrob(p.title)
    ## print(g.title)
    ## g.title <- gtable::gtable_filter(g.title, "title")
    grobs[[length(grobs)+1]] <- g.title
  }
  grid.arrange(grobs = grobs, layout_matrix = layout, heights = heights, widths = widths)
  return(NULL)
}

## Plot BRAF mutations vs AUCs for P-0850
## oncomap.ccls <- unique(oncomap.tbl$ccl_name)
oncomap.ccls <- as.character(unique(sample.info.tbl$ccl_name[sample.info.tbl$Oncomap == "yes"]))

oncomap.braf <- data.frame(ccl_name = oncomap.ccls, BRAF = 0)
rownames(oncomap.braf) <- oncomap.braf$ccl_name
braf.mutated.samples <- unique(oncomap.tbl$ccl_name[oncomap.tbl$Hugo_Symbol == "BRAF"])
oncomap.braf[braf.mutated.samples,"BRAF"] <- 1
ctrpv1_p0850_aucs <- ctrpv1_aucs[ctrpv1_aucs$cpd_name == "P-0850",]
oncomap.braf <- merge(oncomap.braf, ctrpv1_p0850_aucs, by="ccl_name", all=FALSE)
oncomap.braf <- oncomap.braf[,c("ccl_name", "BRAF", "area_under_curve")]


## Propagate braf using network

## Create tidy version of oncomap -- this will only include genes with mutations (i.e., "1"s)
ctrpv1.oncomap <- oncomap.tbl[,c("Hugo_Symbol", "ccl_name")]
colnames(ctrpv1.oncomap) <- c("Gene", "Sample")
ctrpv1.oncomap$Value <- 1
gene <- ctrpv1.oncomap$Gene[1]
## Add all of the assayed samples for one of the genes--when we fill in missing values
## this will force the samples to exist for all of the genes
tmp <- data.frame(Gene = rep(gene, length(oncomap.ccls)), Sample = oncomap.ccls)
ctrpv1.oncomap <- merge(ctrpv1.oncomap, tmp, all=TRUE)
ctrpv1.oncomap$Value[is.na(ctrpv1.oncomap$Value)] <- 0

## Translate to a matrix, which will fill in zeros
ctrpv1.oncomap.mat <- acast(formula = Sample ~ Gene, data = ctrpv1.oncomap, fill = 0, value.var = "Value", fun.aggregate = function(x) ifelse(any(x > 0), 1, 0))

## Now translate back to tidy 
ctrpv1.oncomap <- melt(ctrpv1.oncomap.mat)
colnames(ctrpv1.oncomap) <- c("Sample", "Gene", "Value")

## Create tidy version of TES-A capture sequencing -- this will only include genes with mutations (i.e., "1"s)
flag <- tesa.tbl$Variant_Classification == "Missense_Mutation"
## ctrpv1.tesa <- tesa.tbl[flag,c("Hugo_Symbol", "ccl_name")]
ctrpv1.tesa <- tesa.tbl[,c("Hugo_Symbol", "ccl_name")]
colnames(ctrpv1.tesa) <- c("Gene", "Sample")

ctrpv1.tesa <- unique(ctrpv1.tesa)
## tab <- table(ctrpv1.tesa$Sample)
## plot(hist(tab[tab < 500], breaks=seq(from=0, to=500,by=50), probability=FALSE))

ctrpv1.tesa$Value <- 1
gene <- ctrpv1.tesa$Gene[1]
## Add all of the assayed samples for one of the genes--when we fill in missing values
## this will force the samples to exist for all of the genes
tmp <- data.frame(Gene = rep(gene, length(oncomap.ccls)), Sample = oncomap.ccls)
ctrpv1.tesa <- merge(ctrpv1.tesa, tmp, all=TRUE)
ctrpv1.tesa$Value[is.na(ctrpv1.tesa$Value)] <- 0

## Translate to a matrix, which will fill in zeros
ctrpv1.tesa.mat <- acast(formula = Sample ~ Gene, data = ctrpv1.tesa, fill = 0, value.var = "Value", fun.aggregate = function(x) ifelse(any(x > 0), 1, 0))

## Now translate back to tidy 
ctrpv1.tesa <- melt(ctrpv1.tesa.mat)
colnames(ctrpv1.tesa) <- c("Sample", "Gene", "Value")


ctrpv1.pheno.data <- ctrpv1_aucs
colnames(ctrpv1.pheno.data) <- c("Sample", "Phenotype", "Response")

## Read in ctrp v1 drug targets.
target.file <- "v10.M1.informer_set.txt"
ctrpv1.target.data <- read.table(target.file, sep="\t", header=TRUE, comment.char="", quote="\"", as.is=TRUE)
ctrpv1.target.data <- ctrpv1.target.data[,c("cpd_name", "gene_symbol_of_protein_target")]
colnames(ctrpv1.target.data) <- c("Phenotype", "Gene")
library(plyr)
ctrpv1.target.data <- ddply(.data = ctrpv1.target.data, .variables = c("Phenotype", "Gene"),
                            .fun = function(df) {
                              genes <- unlist(strsplit(x=as.character(df$Gene[1]), split=";[ ]*"))
                              data.frame(Phenotype = rep(df$Phenotype[1], length(genes)), Gene = genes)
                            })

## target.genes <- unique(ctrpv1.target.data$Gene)
target.genes <- unique(ctrpv1.oncomap$Gene)

featureData = ctrpv1.oncomap
sampleOutcomeData = ctrpv1.pheno.data
phenoFeatureData = ctrpv1.target.data

testDrugs = c("P-0850")

## Create the n3 fendR object
fObj.ctrpv1 <- n3FendR(network = network.file,
                       featureData = ctrpv1.oncomap,
                       sampleOutcomeData = ctrpv1.pheno.data,
                       phenoFeatureData = ctrpv1.target.data,
                       target.genes = target.genes,
                       testDrugs = testDrugs)

cat("Enginerring NN = 0 features\n")
fObj.nn0.ctrpv1 <- createNewFeaturesFromNetwork(fObj.ctrpv1, num.degrees = 0)

## Create the engineered features (nearest neighbor degree = 1)
cat("Enginerring NN = 1 features\n")
fObj.nn1.ctrpv1 <- createNewFeaturesFromNetwork(fObj.ctrpv1, num.degrees = 1)

cat("Get engineered NN = 0 sparse matrix\n")
m.eng.nn0.ctrpv1 <- engineeredSparseResponseMatrix(fObj.nn0.ctrpv1)

cat("Get engineered NN = 1 sparse matrix\n")
m.eng.nn1.ctrpv1 <- engineeredSparseResponseMatrix(fObj.nn1.ctrpv1)

m.orig.ctrpv1 <- originalResponseMatrix(fObj.nn0.ctrpv1, limit.to.engineered.genes = TRUE)

make.auc.vs.gene.plots <- function(m.eng.nn0, m.eng.nn1, m.orig, m.oncomap, gene, prefix, response.var = "AUC") {

  ## NB: nn0 is no change from the 0/1 response
  cat(paste0("Fitting ", response.var, " ~ . to nn0\n"))
  tmp <- as.data.frame(as.matrix(m.eng.nn0$feature.mat))
  tmp[,response.var] <- m.eng.nn0$response.mat$Response
  lm.fit <- lm(formula = as.formula(paste0(response.var, " ~ .")), data = tmp)
  sum <- summary(lm.fit)
  f <- sum$fstatistic
  lm.pval <- unname(pf(f[1],f[2],f[3],lower.tail=F))

  if((gene %in% colnames(m.eng.nn0$feature.mat)) && (length(unique(m.eng.nn0$feature.mat[,gene])) > 1)) {
    cat(paste0("Fitting ", response.var, "Response ~ ", gene, " to nn0\n"))
    tmp <- m.eng.nn0$response.mat
    tmp[,gene] <- m.eng.nn0$feature.mat[,gene]
    colnames(tmp)[which(colnames(tmp)=="Response")] <- response.var
    ## tmp <- tmp[, !(colnames(tmp) %in% c("Sample", "Phenotype"))]
    lm.fit <- lm(formula = as.formula(paste(response.var, gene, sep = " ~ ")), data = tmp)
    sum <- summary(lm.fit)
    print(coefficients(sum))
    gene.pval <- coefficients(sum)[gene, 4]
  
    png(paste0(prefix, "-nn0.png"))
    n <- nrow(tmp)
    plot.auc.vs.annotation(tmp, anno.col=gene, sample.col="Sample", auc.col=response.var, main=paste0(response.var, " ~ ", gene, " (univariate NN0; p = ", format(gene.pval, digits=2), "; n = ", n, ")\n", response.var, " ~ all genes (multivariate NN0; p = ", format(lm.pval, digits=2), "; n = ", n, ")"))
    d <- dev.off()
  }
  
  cat(paste0("Fitting ", response.var, " ~ . to nn1\n"))
  tmp <- as.data.frame(as.matrix(m.eng.nn1$feature.mat))
  tmp[, response.var] <- m.eng.nn1$response.mat$Response
  lm.fit <- lm(formula = as.formula(paste0(response.var, " ~ .")), data = tmp)
  sum <- summary(lm.fit)
  f <- sum$fstatistic
  lm.pval <- unname(pf(f[1],f[2],f[3],lower.tail=F))

  if((gene %in% colnames(m.eng.nn1$feature.mat)) && (length(unique(m.eng.nn1$feature.mat[,gene])) > 1)) {
    cat(paste0("Fitting ", response.var, " ~ ", gene, " to nn1\n"))
    tmp <- m.eng.nn1$response.mat
    tmp[,gene] <- m.eng.nn1$feature.mat[,gene]
    colnames(tmp)[which(colnames(tmp)=="Response")] <- response.var
    ## tmp <- tmp[, !(colnames(tmp) %in% c("Sample", "Phenotype"))]
    lm.fit <- lm(formula = as.formula(paste(response.var, gene, sep = " ~ ")), data = tmp)
    sum <- summary(lm.fit)
    print(coefficients(sum))
    gene.pval <- coefficients(sum)[gene, 4]
  
    png(paste0(prefix, "-nn1.png"))
    n <- nrow(tmp)
    plot.auc.vs.annotation(tmp, anno.col=gene, sample.col="Sample", auc.col=response.var, main=paste0(response.var, " ~ ", gene, " (univariate NN1; p = ", format(gene.pval, digits=2), "; n = ", n, ")\n", response.var, " ~ all genes (multivariate NN1; p = ", format(lm.pval, digits=2), "; n = ", n, ")"))
    d <- dev.off()
  }
  
  ## Calculate pvalue for AUC ~ gene
  if((gene %in% colnames(m.oncomap)) && (length(unique(m.oncomap[,gene])) > 1)) {
    cat(paste0("Fitting ", response.var, " ~ ", gene, " to orig\n"))
    lm.fit <- lm(formula = as.formula(paste(response.var, gene, sep = " ~ ")), data = m.oncomap)
    sum <- summary(lm.fit)
    print(coefficients(sum))
    flag <- grepl(rownames(coefficients(sum)), pattern=gene)
    gene.pval <- coefficients(sum)[flag, 4]

    ## Calculate pvalue for AUC ~ gene mutations
    cat(paste0("Fitting ", response.var, " ~ . to orig\n"))
    m.orig <- m.orig[, !(colnames(m.orig) %in% c("Sample", "Phenotype"))]
    colnames(m.orig)[which(colnames(m.orig) == "Response")] <- response.var
    lm.fit <- lm(formula = as.formula(paste0(response.var, " ~ .")), data = m.orig)
    sum <- summary(lm.fit)
    f <- sum$fstatistic
    lm.pval <- unname(pf(f[1],f[2],f[3],lower.tail=F))
    print(sum)
    png(paste0(prefix, "-orig.png"))
    n <- nrow(m.oncomap)
    plot.auc.vs.annotation(m.oncomap, anno.col=gene, sample.col="ccl_name", auc.col = response.var, main=paste0(response.var, " ~ ", gene, " (univariate Orig; p = ", format(gene.pval, digits=2), "; n = ", n, ")\n", response.var, " ~ all genes (multivariate Orig; p = ", format(lm.pval, digits=2), "; n = ", n, ")"))
    d <- dev.off()
  }
}

## NB: nn0 is no change from the 0/1 response
tmp <- as.data.frame(as.matrix(m.eng.nn0.ctrpv1$feature.mat))
tmp$Response <- m.eng.nn0.ctrpv1$response.mat$Response
lm.fit <- lm(formula = Response ~ ., data = tmp)
sum <- summary(lm.fit)
f <- sum$fstatistic
lm.pval <- unname(pf(f[1],f[2],f[3],lower.tail=F))

tmp <- m.eng.nn0.ctrpv1$response.mat
tmp$braf <- m.eng.nn0.ctrpv1$feature.mat[,"BRAF"]
lm.fit <- lm(formula = Response ~ braf, data = tmp)
sum <- summary(lm.fit)
braf.pval <- coefficients(sum)["braf", 4]

png("braf-vs-P-0850-nn0.png")
plot.auc.vs.annotation(tmp, anno.col="braf", sample.col="Sample", auc.col="Response", main=paste0("AUC ~ BRAF (NN0; p = ", format(braf.pval, digits=2), ")\nAUC ~ mutations (NN0; p = ", format(lm.pval, digits=2), ")"))
d <- dev.off()

tmp <- as.data.frame(as.matrix(m.eng.nn1.ctrpv1$feature.mat))
tmp$Response <- m.eng.nn1.ctrpv1$response.mat$Response
lm.fit <- lm(formula = Response ~ ., data = tmp)
sum <- summary(lm.fit)
f <- sum$fstatistic
lm.pval <- unname(pf(f[1],f[2],f[3],lower.tail=F))

tmp <- m.eng.nn1.ctrpv1$response.mat
tmp$braf <- m.eng.nn1.ctrpv1$feature.mat[,"BRAF"]
lm.fit <- lm(formula = Response ~ braf, data = tmp)
sum <- summary(lm.fit)
braf.pval <- coefficients(sum)["braf", 4]

png("braf-vs-P-0850-nn1.png")
plot.auc.vs.annotation(tmp, anno.col="braf", sample.col="Sample", auc.col="Response", main=paste0("AUC ~ BRAF (NN1; p = ", format(braf.pval, digits=2), ")\nAUC ~ mutations (NN1; p = ", format(lm.pval, digits=2), ")"))
d <- dev.off()

## Calculate pvalue for AUC ~ braf
lm.fit <- lm(formula = area_under_curve ~ braf, data = oncomap.braf)
sum <- summary(lm.fit)
coefficients(sum)
braf.pval <- coefficients(sum)["braf", 4]

## Calculate pvalue for AUC ~ gene mutations
m.orig.ctrpv1 <- originalResponseMatrix(fObj.nn0.ctrpv1, limit.to.engineered.genes = TRUE)
m.orig.ctrpv1 <- m.orig.ctrpv1[, !(colnames(m.orig.ctrpv1) %in% c("Sample", "Phenotype"))]
lm.fit <- lm(formula = Response ~ ., data = m.orig.ctrpv1)
sum <- summary(lm.fit)
f <- sum$fstatistic
lm.pval <- unname(pf(f[1],f[2],f[3],lower.tail=F))

png("braf-vs-P-0850-orig.png")
plot.auc.vs.annotation(oncomap.braf, anno.col="braf", sample.col="ccl_name", main=paste0("AUC ~ BRAF (Orig; p = ", format(braf.pval, digits=2), ")\nAUC ~ mutations (Orig; p = ", format(lm.pval, digits=2), ")"))
d <- dev.off()


make.auc.vs.gene.plots(m.eng.nn0.ctrpv1, m.eng.nn1.ctrpv1, m.orig.ctrpv1, oncomap.braf, "BRAF", "BRAF-vs-P-0850") 

do.drug.gene <- function(aucs, feature.data, feature.data.mat, pheno.data, target.data, genes, drug, auc.drug.col = "cpd_name", auc.sample.col = "ccl_name", auc.response.col = "AUC", feature.data.gene.col = "Gene", postfix = NULL) {
  drug_aucs <- aucs[aucs[, auc.drug.col] == drug,]
  testDrugs = c(drug)
  target.genes <- unique(feature.data[, feature.data.gene.col])

  ## Create the n3 fendR object
  fObj.drug <- n3FendR(network = network.file,
                       featureData = feature.data,
                       sampleOutcomeData = pheno.data,
                       phenoFeatureData = target.data,
                       target.genes = target.genes,
                       testDrugs = testDrugs)

  cat("Creasting NN0 features")
  fObj.nn0.drug <- createNewFeaturesFromNetwork(fObj.drug, num.degrees = 0)
  cat("Creasting NN1 features")
  fObj.nn1.drug <- createNewFeaturesFromNetwork(fObj.drug, num.degrees = 1)

  m.eng.nn0.drug <- engineeredSparseResponseMatrix(fObj.nn0.drug)
  m.eng.nn1.drug <- engineeredSparseResponseMatrix(fObj.nn1.drug)
  m.orig.drug <- originalResponseMatrix(fObj.nn0.drug, limit.to.engineered.genes = TRUE)

  for(gene in genes) {
    feature.data.gene <- NULL
    if(gene %in% colnames(feature.data.mat)) {
      feature.data.gene <- feature.data.mat[,gene,drop=F]
      feature.data.gene <- cbind(feature.data.gene, rownames(feature.data.mat))
      colnames(feature.data.gene)[ncol(feature.data.gene)] <- auc.sample.col
    } else {
      feature.data.gene <- data.frame(ccl_name = rownames(feature.data.mat))
      colnames(feature.data.gene) <- auc.sample.col
      feature.data.gene <- cbind(feature.data.gene, rep(0, nrow(feature.data.mat)))
      colnames(feature.data.gene)[ncol(feature.data.gene)] <- gene
    } 
    feature.data.gene <- merge(feature.data.gene, drug_aucs, by=auc.sample.col, all=FALSE)
    feature.data.gene <- feature.data.gene[,c(auc.sample.col, gene, auc.response.col)]

    prefix <- paste0(gene, "-vs-", drug)
    if(!is.null(postfix)) {
      prefix <- paste0(prefix, "-", postfix)
    }
    make.auc.vs.gene.plots(m.eng.nn0.drug, m.eng.nn1.drug, m.orig.drug, feature.data.gene, gene, prefix)
  }
}


colnames(ctrpv1_aucs)[which(colnames(ctrpv1_aucs) == "area_under_curve")] <- "AUC"

drug <- "P-0850"
genes <- "BRAF"
do.drug.gene(ctrpv1_aucs, ctrpv1.oncomap, ctrpv1.oncomap.mat, ctrpv1.pheno.data, ctrpv1.target.data, genes = genes, drug)

genes <- c("NRAS", "KRAS", "MAP2K1", "MAP2K2")
drug <- "selumetinib"
do.drug.gene(ctrpv1_aucs, ctrpv1.oncomap, ctrpv1.oncomap.mat, ctrpv1.pheno.data, ctrpv1.target.data, genes = genes, drug, postfix = "onco")


## Fig 2A EGFR Lung Onco
lung.samples <- sample.info.tbl$ccl_name[sample.info.tbl$Site.Primary == "lung"]
drug <- "neratinib"
genes <- c("EGFR", "MAP3K8")
do.drug.gene(ctrpv1_aucs[ctrpv1_aucs$ccl_name %in% lung.samples,], ctrpv1.oncomap[ctrpv1.oncomap$Sample %in% lung.samples,], ctrpv1.oncomap.mat[rownames(ctrpv1.oncomap.mat) %in% lung.samples,], ctrpv1.pheno.data[ctrpv1.pheno.data$Sample %in% lung.samples,], ctrpv1.target.data, genes = genes, drug)

## Fig 2A NRAS TES-A selumetinib
## Confirmed the number of mutant cell lines was the same (10) as in the reported enrichment results
## and that the number of samples was similar (117 here vs 112 reported)
## more v10.A1.enrichment_analysis.txt  | grep -E 'cell_line|TES-A' | grep -E 'cell_line|ALL_CCL' | grep -E 'cell_line|EXCLUDE_NON' | grep -E 'cell_line|selumetinib' | grep -E 'cell_line|NRAS' | grep -v CNV
## cell_line_subset	cell_line_exclusion	feature_dataset	cpd_name	enriched_cell_line_feature	number_of_cell_lines	number_of_mutant_cell_lines	enrichment_direction	enrichment_p_value	chi_squared_p_value	squared_max_p_value	log_p_value_score	fdr_q_value	log_q_value_score
## ALL_CCL_LINEAGES	EXCLUDE_NONE	TES-A	selumetinib	NRAS	112	10	sensitive	0.0022838	0.0071124	5.06E-05	-4.296	0.0057891	-2.237
genes <- c("NRAS", "KRAS", "MAP2K1", "MAP2K2")
drug <- "selumetinib"
do.drug.gene(ctrpv1_aucs, ctrpv1.tesa, ctrpv1.tesa.mat, ctrpv1.pheno.data, ctrpv1.target.data, genes = genes, drug, postfix = "tesa")

## Fig 4A CTNNB1 Oncomap vs navitoclax (targets: BCL2, BCL2L1, BCL2L2)
genes <- c("CTNNB1", "BCL2", "BCL2L1", "BCL2L2")
drug <- "navitoclax"
do.drug.gene(ctrpv1_aucs, ctrpv1.oncomap, ctrpv1.oncomap.mat, ctrpv1.pheno.data, ctrpv1.target.data, genes = genes, drug, postfix = "onco")

## Fig 4A AXIN1 TES-A vs navitoclax
genes <- c("AXIN1", "BCL2", "BCL2L1", "BCL2L2")
drug <- "navitoclax"
do.drug.gene(ctrpv1_aucs, ctrpv1.tesa, ctrpv1.tesa.mat, ctrpv1.pheno.data, ctrpv1.target.data, genes = genes, drug, postfix = "tesa")

## Fig 4A CSNK1A1 TES-CNV vs navitoclax -- don't know what TES-CNV is.  So skipping


## AUC file has experiment_id that needs to matched to cell_id and to cell name

file <- "v20.meta.per_experiment.txt"
ctrpv2_expts <- read.table(file, sep="\t", header=TRUE)
ctrpv2_expts <- ctrpv2_expts[,c("experiment_id", "master_ccl_id")]

file <- "v20.meta.per_cell_line.txt"
ctrpv2_ccls <- read.table(file, sep="\t", header=TRUE)
ctrpv2_ccls <- ctrpv2_ccls[,c("master_ccl_id", "ccl_name")]

ctrpv2_ccls <- merge(ctrpv2_ccls, ctrpv2_expts, all=FALSE)

ctrpv2_aucs <- merge(ctrpv2_aucs, ctrpv2_ccls, all=FALSE)


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
