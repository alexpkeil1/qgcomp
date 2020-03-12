# utility functions that are not part of the base functions but serve miscellaneous purposes


# allow dependencies that don't require installation at initial installation of qgcomp
.qgc.require <- function (package, message = paste("loading required package (", 
                                                   package, ") failed", sep = "")){
  if (!requireNamespace(package, quietly = FALSE)) {
    stop(message, call. = FALSE)
  }
  invisible(TRUE)
}

.rmvnorm <- function(n,mu,Sigma){
  # draw from multivariate normal distribution - this is a thin
  #  wrapper for either mvtnorm::rmvnorm or MASS::mvrnorm, 
  #  depending on the moods of those package developers
  #.qgc.require("mvtnorm")
  #draw = mvtnorm::rmvnorm(n, mu, Sigma)
  .qgc.require("MASS")
  draw = MASS::mvrnorm(n, mu, Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  draw
}

construction <- function(i="", j=""){
  msg ="This function is not yet fully implemented"
  if(j != "") msg = paste0(msg, ": ", j)
  if(i %in% c("msg", "message", "warning", "wrn", "warn")) warning(msg)
  else stop(msg)
}

# fake cox family function
cox <- function(){
  obj = binomial(link="log")
  obj$family="cox"
  obj
}

# fake zi family function
zi <- function(){
  obj = binomial(link="log")
  obj$family="zi"
  obj
}


#printing of zero-inflated results

printZI <- function(x, showweights=TRUE, ...){
  if(class(x$fit) == "zeroinfl"){
    if(!is.null(x$pos.size$count) & showweights) {
      cat("Prob(Y ~ count):\n")
      cat(paste0("Scaled effect size (positive direction, sum of positive coefficients = ", signif(x$pos.size$count, 3) , ")\n"))
      if (length(x$pos.weights$count) > 0) {
        print(x$pos.weights$count, digits = 3)
      } else cat("None\n")
      cat("\n")
    }
    if(!is.null(x$neg.size$count) & showweights) {
      cat(paste0("Scaled effect size (negative direction, sum of negative coefficients = ", signif(-x$neg.size$count, 3) , ")\n"))
      if (length(x$neg.weights$count) > 0) {
        print(x$neg.weights$count, digits = 3)
      } else cat("None\n")
      cat("\n")
    }
    if(!is.null(x$pos.size$zero) & showweights) {
      cat("Prob(Y ~ zero/count):\n")
      cat(paste0("Scaled effect size (positive direction, sum of positive coefficients = ", signif(x$pos.size$zero, 3) , ")\n"))
      if (length(x$pos.weights$zero) > 0) {
        print(x$pos.weights$zero, digits = 3)
      } else cat("None\n")
      cat("\n")
    }
    if(!is.null(x$neg.size$zero) & showweights) {
      cat(paste0("Scaled effect size (negative direction, sum of negative coefficients = ", signif(-x$neg.size$zero, 3) , ")\n"))
      if (length(x$neg.weights$zero) > 0) {
        print(x$neg.weights$zero, digits = 3)
      } else cat("None\n")
      cat("\n")
    }
    
    if(x$fit$dist %in% c("poisson", "negbin")){
      estimand <- 'OR/RR'
      cat(paste0("Mixture log(",estimand,")", ifelse(x$bootstrap, " (bootstrap CI)", " (Delta method CI)"), ":\n\n"))
    }
    #if(x$fit$dist=="gaussian"){
    #  estimand <- 'log(OR)/mean diff'
    #  cat(paste0("Mixture ",estimand,"", ifelse(x$bootstrap, " (bootstrap CI)", " (Delta method CI)"), ":\n\n"))
    #}
    testtype = "Z"
    
    pdat <- list()
    for(modtype in names(x$psi)){
      pdat[[modtype]] <- cbind(Estimate=coef(x)[[modtype]], "Std. Error"=sqrt(x$var.coef[[modtype]]), 
                               "Lower CI"=x$ci.coef[[modtype]][,1], "Upper CI"=x$ci.coef[[modtype]][,2], 
                               "test"=x$zstat[[modtype]], "Pr(>|z|)"=x$pval[[modtype]])
      colnames(pdat[[modtype]])[which(colnames(pdat)=="test")] = eval(paste(testtype, "value"))
      #colnames(pdat[[modtype]])[5] = eval(paste(testtype, "value"))
      numpsi = length(x$psi[[modtype]])
      if(numpsi>0) rnm = c("(Intercept)", c(paste0('psi',1:max(1, numpsi))))
      if(numpsi==0) rnm = c("(Intercept)")
      rownames(pdat[[modtype]]) <- rnm
      cat(paste0("Prob(Y ~ ", modtype,"):\n"))
      printCoefmat(pdat[[modtype]],has.Pvalue=TRUE,tst.ind=5L,signif.stars=FALSE, cs.ind=1L:2)
    }
  }
}

summaryZI <- function(x){
  if(class(x$fit) == "zeroinfl"){
    if(x$fit$dist %in% c("poisson", "negbin")){
      estimand <- 'OR/RR'
      cat(paste0("Mixture log(",estimand,")", ifelse(x$bootstrap, " (bootstrap CI)", " (Delta method CI)"), ":\n\n"))
    }
    #if(x$fit$dist=="gaussian"){
    #  estimand <- 'log(OR)/mean diff'
    #  cat(paste0("Mixture ",estimand,"", ifelse(x$bootstrap, " (bootstrap CI)", " (Delta method CI)"), ":\n\n"))
    #}
    testtype = "Z"
    
    pdat <- list()
    for(modtype in names(x$psi)){
      pdat[[modtype]] <- cbind(Estimate=coef(x)[[modtype]], "Std. Error"=sqrt(x$var.coef[[modtype]]), 
                               "Lower CI"=x$ci.coef[[modtype]][,1], "Upper CI"=x$ci.coef[[modtype]][,2], 
                               "test"=x$zstat[[modtype]], "Pr(>|z|)"=x$pval[[modtype]])
      colnames(pdat[[modtype]])[5] = eval(paste(testtype, "value"))
      numpsi = length(x$psi[[modtype]])
      if(numpsi>0) rnm = c("(Intercept)", c(paste0('psi',1:max(1, numpsi))))
      if(numpsi==0) rnm = c("(Intercept)")
      rownames(pdat[[modtype]]) <- rnm
      cat(paste0("Prob(Y ~ ", modtype,"):\n"))
      #printCoefmat(pdat[[modtype]],has.Pvalue=TRUE,tst.ind=5L,signif.stars=FALSE, cs.ind=1L:2)
    }
  }
  list(coeffients=pdat)
}

.dgm_quantized <- function(
  N = 100,             # sample size
  b0=0,                # baseline expected outcome (model intercept) 
  coef=c(1,0,0,0),     # beta coefficients for X in the outcome model
  ncor=0,              # Number of correlated exposures
  corr=0.75            # Pearson/spearman (the same here) correlation
){
  #'
  # simulate under data structure where WQS/qgcomp is the truth: 
  #  e.g. a multivariate exposure with multinomial distribution
  #  and an outcome that is a linear function of exposure scores
  
  p = length(coef)
  if(ncor >= p) ncor = p-1
  X = matrix(nrow=N, ncol=p)
  xmaster = sample(rep(0:3, length.out=N), N, replace=FALSE)
  for(k in 1:p){
    newx = numeric(N)
    c1 = as.logical(rbinom(N, 1, sqrt(corr)))
    newx[which(c1)] = xmaster[which(c1)]
    newx[which(!c1)] = sample(xmaster[which(!c1)])
    if(k<=(ncor+1)){
      X[,k] = newx
    } else X[,k] = sample(xmaster)
  }
  mu <- X %*% coef
  y = rnorm(N) + mu
  colnames(X) <- paste0("x", 1:p)
  data.frame(X,y)
}
