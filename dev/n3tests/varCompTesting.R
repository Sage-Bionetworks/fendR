library(varComp)

##simulate data
n <- 5e2
m <- 1e3
x <- matrix(rnorm(n*m),n,m)
beta <- rnorm(m)
predVar <- var(x%*%beta)
#(explain 50% of variance)
set.seed(1)
e <- rnorm(n,0,sqrt(predVar))
y <- x%*%beta+e

varCompWrapperFunction(outcome=y,features=x)
library(glmnet)
foobar <- cv.glmnet(y=y,x=x,alpha=0)


varCompWrapperFunction <- function(outcome,features){
  library(varComp)
  ####features have to be in observations as rows and features as column format
  kernel <- cov(t(features),
           use = 'pairwise.complete.obs')
  alternativeModel <- varComp::varComp(outcome~1,
                          varcov=kernel)
  nullModel <- varComp::varComp(outcome~1)

  model <- list()
  model$totalVariance <- as.numeric(var(outcome))
  model$varianceExplained <- as.numeric(alternativeModel$varComps[1])
  model$percentVarianceExplained <- model$varianceExplained/model$totalVariance
  model$errorVarianceExplained <- as.numeric(alternativeModel$sigma2)
  model$LRT <- 2*(logLik(alternativeModel)[1]-logLik(nullModel)[1])

  pVarianceLRT <- function(x){
    x <- round(x,8)
    if(x==0){
      return(0.5)
    }else{
      return(0.5*pchisq(x,1,lower.tail=F))
    }
  }

  model$pvalue <-pVarianceLRT(model$LRT)
  return(model)
}


