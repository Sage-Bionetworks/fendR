library(varComp)

##simulate data

x <- matrix(rnorm(1e6),1e3,1e3)
beta <- rnorm(1e3)
predVar <- var(x%*%beta)
#(explain 50% of variance)
set.seed(1)
e <- rnorm(1e3,0,sqrt(predVar)/2)
y <- x%*%beta+e
var(y)

G <- cov(t(x),use=)

foo <- varComp::varComp(y ~ 1,varcov = G)
summary(foo)

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

fendR::varCompWrapperFunction(outcome=y,features=x)
