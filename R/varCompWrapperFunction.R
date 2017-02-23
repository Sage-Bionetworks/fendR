varCompWrapperFunction <- function(outcome,features){
  library(varComp)
  ####features have to be in observations as rows and features as column format
  kernel <- cov(t(features),
                use = 'pairwise.complete.obs')
  foo <- varComp::varComp(outcome~1,
                          varcov=kernel)

  model <- list()
  model$totalVariance <- as.numeric(var(outcome))
  model$varianceExplained <- as.numeric(foo$varComps[1])
  model$percentVarianceExplained <- model$varianceExplained/model$totalVariance
  model$errorVarianceExplained <- as.numeric(foo$sigma2)
  model$zscore <- model$varianceExplained/model$errorVarianceExplained
  model$pvalue <- pnorm(model$zscore,lower.tail=F)
  return(model)
}
