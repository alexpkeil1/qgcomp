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

.dgm_quantized <- function(
  N = 100,             # sample size
  b0=0,                # baseline expected outcome (model intercept)
  coef=c(1,0,0,0),     # beta coefficients for X in the outcome model
  ncor=0,              # Number of correlated exposures
  corr=0.75,            # Pearson/spearman (the same here) correlation
  sigma=1.0            # standard deviation of error term (true residuals)
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
  y = rnorm(N,0,sigma) + mu + b0
  colnames(X) <- paste0("x", 1:p)
  data.frame(X,y)
}
