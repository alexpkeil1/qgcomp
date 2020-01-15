library(pscl)

ziptest <- function(f, data, expnms=NULL, q=4, breaks=NULL, id=NULL, alpha=0.05, bayes=FALSE, ...){
  data("bioChemists", package = "pscl")
  # art   fem     mar kid5  phd ment
  #   0   Men Married    0 2.52    7
  #   0 Women  Single    0 2.05    6
  #   0 Women  Single    0 3.75    6
  #   0   Men Married    1 1.18    3
  #   0 Women  Single    0 3.75   26
  #   0 Women Married    2 3.59    2
  f = art ~ fem + kid5 + ment + phd | kid5
  expnms = c("phd", "ment")
  qdata = qgcomp:::quantize(bioChemists, expnms=expnms)
  
  
  res = zeroinfl(f, data = qdata$data)
  # non bootstrapped
  coef1 = res$coefficients$count[expnms]
  vc1 = vcov(res, "count")
  coef2 = res$coefficients$zero[expnms]
  vc2 = vcov(res, "zero")
  psi
  "Count model coefficients (poisson with log link):"
  sum(coef1)
  qgcomp:::se_comb(expnms, vc1)
  "Zero-inflation model coefficients (binomial with logit link):"
  sum(coef2)
  qgcomp:::se_comb(expnms, vc2)
  
  # bootstrapped
    ndat3 <- ndat2 <- ndat1 <- ndat0 <- qdata$data
    ndat0$phd=0
    ndat0$ment=0
    ndat0$psi=0
    ndat1$phd=1
    ndat1$ment=1
    ndat1$psi=1
    ndat2$phd=2
    ndat2$ment=2
    ndat2$psi=2
    ndat3$phd=3
    ndat3$ment=3
    ndat3$psi=3
    ndat = data.frame(rbind(ndat0, ndat1, ndat2, ndat3))
    totn = nrow(ndat)
    # new predicted class using class probabilities (more general?)
      classprob = predict(res, newdata = ndat, type="prob")
      ncats = ncol(classprob)
      sum(rmultinom(1,ncats,classprob[,]))
      resres = apply(classprob[,], 1, function(x) -1+which.max(rmultinom(1, 1, x)))
      resres
    # new predicted class using count/zero preds (faster?)
      pmfg0 = predict(res, newdata = ndat, type="count")
      pmf0  = predict(res, newdata = ndat, type="zero")
      resres2 = rbinom(totn, 1, 1-pmf0)*rpois(totn, pmfg0)
      mean(resres)
      mean(resres2)
      mean(qdata$data$art)
      prop.table(table(resres, ndat$psi), margin = 2)
      prop.table(table(resres2, ndat$psi), margin = 2)
      table(qdata$data$art)
    ndat$art = resres
    # new regression
    margf = art ~ fem + kid5 + psi | kid5 + psi
    margres = zeroinfl(margf, data = ndat)
  
  # comparison
  margres
  "Count model coefficients (poisson with log link):"
  sum(coef1)
  qgcomp:::se_comb(expnms, vc1)
  "Zero-inflation model coefficients (binomial with logit link):"
  sum(coef2)
  qgcomp:::se_comb(expnms, vc2)
  
  ##############################
}






qgcomp.zi.noboot <- function(f, data, expnms=NULL, q=4, breaks=NULL, id=NULL, alpha=0.05, bayes=FALSE, ...){
  #' @title estimating the parameters of a zero-inflated marginal structural model (MSM) based on 
  #' g-computation with quantized exposures
  #'
  #' @description This function mimics the output of a weighted quantile sums regression in 
  #' large samples. 
  #' 
  #' @details A zero-inflated version of quantile g-computation based on the implementation in the
  #' 'pscl' package. A zero-inflated distribution is a mixture distribution in which one of the
  #' distributions is a point mass at zero (with probability given by a logistoic model), and the 
  #' other distribution is a discrete or continuous distribution.
  #' This estimates the effect of a joint increase in all exposures on 1) the odds 
  #' of belonging to the "zero" vs. "count" portions of the distribution and/or 2) the rate parameter
  #' for the "count" portion of the distribution.
  #' 
  #' @param f R style formula using syntax from 'pscl' package: depvar ~ indvars_count | indvars_zero
  #' @param data data frame
  #' @param expnms character vector of exposures of interest
  #' @param q NULL or number of quantiles used to create quantile indicator variables
  #' representing the exposure variables. If NULL, then gcomp proceeds with un-transformed
  #' version of exposures in the input datasets (useful if data are already transformed,
  #' or for performing standard g-computation)
  #' @param breaks (optional) NULL, or a list of (equal length) numeric vectors that 
  #' characterize the minimum value of each category for which to 
  #' break up the variables named in expnms. This is an alternative to using 'q'
  #' to define cutpoints.
  #' @param id (optional) NULL, or variable name indexing individual units of 
  #' observation (only needed if analyzing data with multiple observations per 
  #' id/cluster)
  #' @param alpha alpha level for confidence limit calculation
  #' @param bayes not yet implemented
  #' @param ... arguments to zeroinf (e.g. dist)
  #' @seealso \code{\link[qgcomp]{qgcomp.noboot}}, \code{\link[qgcomp]{qgcomp.cox.noboot}}, 
  #'  and \code{\link[pscl]{zeroinfl}}
  #' @return a qgcompfit object, which contains information about the effect
  #'  measure of interest (psi) and associated variance (var.psi), as well
  #'  as information on the model fit (fit) and information on the 
  #'  weights/standardized coefficients in the positive (pos.weights) and 
  #'  negative (nweight) directions.
  #' @concept variance mixtures
  #' @import stats arm pscl
  #' @export
  #' @examples
  #' set.seed(50)
  #' n=100
  #' dat <- data.frame(y=rbinom(n, 1, 0.5)*rpois(n, 1.2), x1=runif(n), x2=runif(n), z=runif(n))
  #' # poinsson count model, mixture in both portions
  #' qgcomp.zi.noboot(f=y ~ z + x1 + x2 | x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, dist="poisson")
  #' # negative binomial count model, mixture in both portions
  #' qgcomp.zi.noboot(f=y ~ z + x1 + x2 | x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, dist="negbin")
  #' # negative binomial count model, mixture only in the 'count' portion of the model
  #' qgcomp.zi.noboot(f=y ~ z + x1 + x2 | z, expnms = c('x1', 'x2'), data=dat, q=2, dist="negbin")
  
  # list containers
  estb <- vcov_mod <- seb <- tstat <- pvalz <- allterms <- containmix <- pos.weights <- neg.weights <- 
    pos.coef <- neg.coef <-pos.psi <- neg.psi <- pos.size <- neg.size <- wcoef <- ci <- tstat<- list()
  suppressWarnings(testfit <- zeroinfl(f, data = data, control=zeroinfl.control(maxit = 1, EM=FALSE)))
  allterms$count = attr(terms(testfit, "count"), "term.labels")
  allterms$zero = attr(terms(testfit, "zero"), "term.labels")
  if (is.null(expnms)) {
    expnms <- attr(terms(testfit), "term.labels")
    message("Including all model terms as exposures of interest (count and zero parts must be identical)\n")      
  }
  lin = checknames(expnms)
  if(!lin) stop("Model appears to be non-linear: use qgcomp.zi.boot instead")
  if (!is.null(q) | !is.null(breaks)){
    ql <- quantize(data, expnms, q, breaks)
    qdata <- ql$data
    br <- ql$breaks
  } else{
    qdata <- data
    br <- breaks
  }
  if(is.null(id)) {
    # not yet implemented
    id = "id__"
    qdata$id__ = 1:dim(qdata)[1]
  }
  for(modtype in c("count", "zero")) containmix[[modtype]] = all(expnms %in% allterms[[modtype]])
  

  if(!bayes) fit <- zeroinfl(f, data = qdata[,!(names(qdata) %in% id), drop=FALSE], ...)
  if(bayes){
    stop("bayesian zero inflated models not yet implemented")
    #requireNamespace("arm")
    #fit <- bayesglm(f, data = qdata[,!(names(qdata) %in% id), drop=FALSE], ...)
  }
  mod <- summary(fit)
  if((length(setdiff(expnms, rownames(mod$coefficients$count)))>0 & containmix$count) ||
     (length(setdiff(expnms, rownames(mod$coefficients$zero)))>0 & containmix$zero)
     ){
    stop("Model aliasing occurred, likely due to perfectly correlated quantized exposures. 
           Try one of the following:
             1) set 'bayes' to TRUE in the qgcomp function (recommended)
             2) set 'q' to a higher value in the qgcomp function (recommended)
             3) check correlation matrix of exposures, and drop all but one variable in each highly correlated set  (not recommended)
           ")
  }
  for(modtype in names(containmix)){
    if(containmix[[modtype]]){
      estb[[modtype]] = c(fit$coefficients[[modtype]][1], sum(mod$coefficients[[modtype]][expnms,1, drop=TRUE]))
      vc = vcov(fit, modtype)
      vcov_mod[[modtype]] = vc
      seb[[modtype]] = c(sqrt(vc[1,1]), se_comb(expnms, covmat = vc))
      tstat[[modtype]] = estb[[modtype]]/seb[[modtype]]
      ci[[modtype]] = cbind(estb[[modtype]] + seb[[modtype]] * qnorm(alpha / 2), estb[[modtype]] + seb[[modtype]] * qnorm(1 - alpha / 2))
      wcoef[[modtype]] = fit$coefficients[[modtype]][expnms]
      names(wcoef[[modtype]]) <- gsub("_q", "", names(wcoef[[modtype]]))
      pos.coef[[modtype]] <- which(wcoef[[modtype]] > 0)
      neg.coef[[modtype]] <- which(wcoef[[modtype]] <= 0)
      pos.weights[[modtype]] <- abs(wcoef[[modtype]][pos.coef[[modtype]]]) / sum(abs(wcoef[[modtype]][pos.coef[[modtype]]]))
      neg.weights[[modtype]] <- abs(wcoef[[modtype]][neg.coef[[modtype]]]) / sum(abs(wcoef[[modtype]][neg.coef[[modtype]]]))
      pos.psi[[modtype]] <- sum(wcoef[[modtype]][pos.coef[[modtype]]])
      neg.psi[[modtype]] <- sum(wcoef[[modtype]][neg.coef[[modtype]]])
      pos.size[[modtype]] <- sum(abs(wcoef[[modtype]][pos.coef[[modtype]]]))
      neg.size[[modtype]] <- sum(abs(wcoef[[modtype]][neg.coef[[modtype]]]))
    }
  }
  pvalz <- lapply(tstat, function(x) 2 - 2 * pnorm(abs(x)))

  ### TODONE ^
  qx <- qdata[, expnms]
  names(qx) <- paste0(names(qx), "_q")
  res <- list(
    qx = qx, fit = fit, 
    psi = lapply(estb, function(x) x[-1]), 
    var.psi = lapply(seb, function(x) x[-1]^2), 
    covmat.psi = lapply(seb, function(x) c('psi1' = x[-1]^2)),
    ci = lapply(ci, function(x) x[-1,]), 
    coef = estb, 
    var.coef = lapply(seb, function(x) x^2), 
    covmat.coef=lapply(seb, function(x) c('(Intercept)' = x[1]^2, 'psi1' = x[2]^2)),
    ci.coef = ci,
    expnms=expnms, q=q, breaks=br, degree=1,
    pos.psi = pos.psi, 
    neg.psi = neg.psi,
    pos.weights = lapply(pos.weights, function(x) sort(x, decreasing = TRUE)),
    neg.weights = lapply(neg.weights, function(x) sort(x, decreasing = TRUE)), 
    pos.size = pos.size,
    neg.size = neg.size,
    bootstrap=FALSE,
    cov.yhat=NULL
  )
  #if(fit$family$family=='gaussian'){
  #  res$tstat <- tstat
  #  res$df <- df
  #  res$pval <- pval
  #}
  #if(fit$family$family=='binomial'){
    res$zstat <- tstat
    res$pval <- pvalz
  #}
  attr(res, "class") <- "qgcompfit"
  res
}

#qgcomp.zi.boot <- function(f, data, expnms=NULL, q=4, breaks=NULL, id=NULL, alpha=0.05, B=200, 
#                        rr=TRUE, degree=1, seed=NULL, bayes=FALSE, parallel=FALSE, ...){
#  #' @title estimating the parameters of a marginal structural model (MSM) based on 
#  #' g-computation with quantized exposures
#  #'  
#  #' @description This function yields population average effect estimates for 
#  #'   both continuous and binary outcomes
#  #'  
#  #' @details Estimates correspond to the average expected change in the
#  #'  (log) outcome per quantile increase in the joint exposure to all exposures 
#  #'  in `expnms'. Test statistics and confidence intervals are based on 
#  #'  a non-parametric bootstrap, using the standard deviation of the bootstrap
#  #'  estimates to estimate the standard error. The bootstrap standard error is 
#  #'  then used to estimate Wald-type confidence intervals. Note that no bootstrapping
#  #'  is done on estimated quantiles of exposure, so these are treated as fixed
#  #'  quantities
#  #'
#  #' @param f R style formula
#  #' @param data data frame
#  #' @param expnms character vector of exposures of interest
#  #' @param q NULL or number of quantiles used to create quantile indicator variables
#  #' representing the exposure variables. If NULL, then gcomp proceeds with un-transformed
#  #' version of exposures in the input datasets (useful if data are already transformed,
#  #' or for performing standard g-computation)
#  #' @param breaks (optional) NULL, or a list of (equal length) numeric vectors that 
#  #' characterize the minimum value of each category for which to 
#  #' break up the variables named in expnms. This is an alternative to using 'q'
#  #' to define cutpoints.
#  #' @param id (optional) NULL, or variable name indexing individual units of 
#  #' observation (only needed if analyzing data with multiple observations per 
#  #' id/cluster)
#  #' @param alpha alpha level for confidence limit calculation
#  #' @param B integer: number of bootstrap iterations (this should typically be
#  #' >=200, though it is set lower in examples to improve run-time).
#  #' @param rr logical: if using binary outcome and rr=TRUE, qgcomp.boot will 
#  #'   estimate risk ratio rather than odds ratio
#  #' @param degree polynomial basis function for marginal model (e.g. degree = 2
#  #'  allows that the relationship between the whole exposure mixture and the outcome
#  #'  is quadratic.
#  #' @param seed integer or NULL: random number seed for replicable bootstrap results
#  #' @param bayes use underlying Bayesian model (`arm` package defaults). Results
#  #' in penalized parameter estimation that can help with very highly correlated 
#  #' exposures. Note: this does not lead to fully Bayesian inference in general, 
#  #' so results should be interpereted as frequentist.
#  #' @param parallel use (safe) parallel processing from the future and future.apply packages
#  #' @param ... arguments to glm (e.g. family)
#  #' @seealso \code{\link[qgcomp]{qgcomp.noboot}}, and \code{\link[qgcomp]{qgcomp}}
#  #' @return a qgcompfit object, which contains information about the effect
#  #'  measure of interest (psi) and associated variance (var.psi), as well
#  #'  as information on the model fit (fit) and information on the 
#  #'  marginal structural model (msmfit) used to estimate the final effect
#  #'  estimates.
#  #' @concept variance mixtures
#  #' @import stats pscl
#  #' @export
#  #' @examples
#  #' set.seed(30)
#  #' # continuous outcome
#  #' dat <- data.frame(y=rnorm(100), x1=runif(100), x2=runif(100), z=runif(100))
#  #' # Conditional linear slope
#  #' qgcomp.noboot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=4, family=gaussian())
#  #' # Marginal linear slope (population average slope, for a purely linear, 
#  #' #  additive model this will equal the conditional)
#  #' qgcomp.boot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=4, 
#  #'   family=gaussian(), B=10) #increase B to at least 200 in actual examples
#  #'   
#  #' # Population average mixture slope which accounts for non-linearity and interactions
#  #' qgcomp.boot(y ~ z + x1 + x2 + I(x1^2) + I(x2*x1), family="gaussian", 
#  #'  expnms = c('x1', 'x2'), data=dat, q=4, B=6)
#  #'  
#  #' # binary outcome
#  #' dat <- data.frame(y=rbinom(50,1,0.5), x1=runif(50), x2=runif(50), z=runif(50))
#  #' 
#  #' # Conditional mixture OR
#  #' qgcomp.noboot(y ~ z + x1 + x2, family="binomial", expnms = c('x1', 'x2'), 
#  #'   data=dat, q=2)
#  #'   
#  #' #Marginal mixture OR (population average OR - in general, this will not equal the 
#  #' # conditional mixture OR due to non-collapsibility of the OR)
#  #' qgcomp.boot(y ~ z + x1 + x2, family="binomial", expnms = c('x1', 'x2'), 
#  #'   data=dat, q=2, B=6)
#  #'   
#  #' # Population average mixture RR
#  #' qgcomp.boot(y ~ z + x1 + x2, family="binomial", expnms = c('x1', 'x2'), 
#  #'   data=dat, q=2, rr=TRUE, B=6)
#  #'   
#  #' # Population average mixture RR, indicator variable representation of x2
#  #' # note that I(x==...) operates on the quantile-based category of x,
#  #' # rather than the raw value
#  #' res = qgcomp.boot(y ~ z + x1 + I(x2==1) + I(x2==2) + I(x2==3), 
#  #'   family="binomial", expnms = c('x1', 'x2'), data=dat, q=4, rr=TRUE, B=6)
#  #' res$fit  
#  #' plot(res)
#  #' 
#  #' # now add in a non-linear MSM
#  #' res2 = qgcomp.boot(y ~ z + x1 + I(x2==1) + I(x2==2) + I(x2==3), 
#  #'   family="binomial", expnms = c('x1', 'x2'), data=dat, q=4, rr=TRUE, B=6, 
#  #'   degree=2)
#  #' res2$fit  
#  #' res2$msmfit  
#  #' plot(res2)
#  #' # Log risk ratio per one IQR change in all exposures (not on quantile basis)
#  #' dat$x1iqr <- dat$x1/with(dat, diff(quantile(x1, c(.25, .75))))
#  #' dat$x2iqr <- dat$x2/with(dat, diff(quantile(x2, c(.25, .75))))
#  #' # note that I(x>...) now operates on the untransformed value of x,
#  #' # rather than the quantized value
#  #' res2 = qgcomp.boot(y ~ z + x1iqr + I(x2iqr>0.1) + I(x2>0.4) + I(x2>0.9), 
#  #'   family="binomial", expnms = c('x1iqr', 'x2iqr'), data=dat, q=NULL, rr=TRUE, B=6, 
#  #'   degree=2)
#  #' res2
#  #' # using parallel processing
#  #' res2p = qgcomp.boot(y ~ z + x1iqr + I(x2iqr>0.1) + I(x2>0.4) + I(x2>0.9), 
#  #'   family="binomial", expnms = c('x1iqr', 'x2iqr'), data=dat, q=NULL, rr=TRUE, B=6, 
#  #'   degree=2, parallel=TRUE)
#  #' res2p
#  # character names of exposure mixture components
#  if(is.null(seed)) seed = round(runif(1, min=0, max=1e8))
#  if (is.null(expnms)) {
#    expnms <- attr(terms(f, data = data), "term.labels")
#    cat("Including all model terms as exposures of interest\n")      
#  }
#  lin = checknames(expnms)
#  if(!lin) stop("Model appears to be non-linear and I'm having trouble parsing it: 
#                  please use `expnms` parameter to define the variables making up the exposure")
#  if (!is.null(q) & !is.null(breaks)){
#    # if user specifies breaks, prioritize those
#    q <- NULL
#  }
#  if (!is.null(q) | !is.null(breaks)){
#    ql <- quantize(data, expnms, q, breaks)
#    qdata <- ql$data
#    br <- ql$breaks
#    if(is.null(q)){
#      # rare scenario with user specified breaks and q is left at NULL
#      nvals <- length(br[[1]])-1
#    } else{
#      nvals <- q
#    }
#    intvals <- (1:nvals)-1
#  } else {
#    # if( is.null(breaks) & is.null(q)) # also includes NA
#    qdata <- data
#    # if no transformation is made (no quantiles, no breaks given)
#    # then draw distribution values from quantiles of all the exposures
#    # pooled together
#    # TODO: allow user specification of this
#    cat("\nNote: using quantiles of all exposures combined in order to set 
#          proposed intervention values for overall effect (25th, 50th, 75th %ile)\n")
#    intvals = as.numeric(quantile(unlist(data[,expnms]), c(.25, .5, .75)))
#    br <- NULL
#  }
#  if(is.null(id)) {
#    id <- "id__"
#    qdata$id__ <- 1:dim(qdata)[1]
#  }
#  ###
#  msmfit <- msm.fit(f, qdata, intvals, expnms, rr, main=TRUE,degree=degree, id=id, bayes, ...)
#  # main estimate  
#  #estb <- as.numeric(msmfit$msmfit$coefficients[-1])
#  estb <- as.numeric(msmfit$msmfit$coefficients)
#  #bootstrap to get std. error
#  nobs <- dim(qdata)[1]
#  nids <- length(unique(qdata[,id, drop=TRUE]))
#  starttime = Sys.time()
#  psi.only <- function(i=1, f=f, qdata=qdata, intvals=intvals, expnms=expnms, rr=rr, degree=degree, nids=nids, id=id, ...){
#    if(i==2){
#      timeiter = as.numeric(Sys.time() - starttime)
#      if((timeiter*B/60)>0.5) cat(paste0("Expected time to finish: ", round(B*timeiter/60, 2), " minutes \n"))
#    }
#    bootids <- data.frame(temp=sort(sample(unique(qdata[,id, drop=TRUE]), nids, replace = TRUE)))
#    names(bootids) <- id
#    qdata_ <- merge(qdata,bootids, by=id, all.x=FALSE, all.y=TRUE)
#    ft = msm.fit(f, qdata_, intvals, expnms, rr, main=FALSE, degree, id, bayes,
#                 ...)
#    yhatty = data.frame(yhat=predict(ft$msmfit), psi=ft$msmfit$data[,"psi"])
#    as.numeric(
#      c(ft$msmfit$coefficients, with(yhatty, tapply(yhat, psi, mean)))
#      #msm.fit(f, qdata_, intvals, expnms, rr, main=FALSE, degree, id, bayes,
#      #        ...)$msmfit$coefficients[-1]
#    )
#  }
#  set.seed(seed)
#  if(parallel){
#    Sys.setenv(R_FUTURE_SUPPORTSMULTICORE_UNSTABLE="quiet")
#    future::plan(strategy = future::multiprocess)
#    bootsamps <- future.apply::future_sapply(X=1:B, FUN=psi.only,f=f, qdata=qdata, intvals=intvals, 
#                                             expnms=expnms, rr=rr, degree=degree, nids=nids, id=id, ...)
#    
#    future::plan(future::sequential)
#  }else{
#    bootsamps <- sapply(X=1:B, FUN=psi.only,f=f, qdata=qdata, intvals=intvals, 
#                        expnms=expnms, rr=rr, degree=degree, nids=nids, id=id, ...)
#    
#  }
#  #if(is.null(dim(bootsamps))) {
#  #  seb <- sd(bootsamps)
#  #  covmat <- var(bootsamps)
#  #  names(covmat) <- 'psi1'
#  #}else{
#  hats = t(bootsamps[-c(1:(degree+1)),])
#  cov.yhat = cov(hats)
#  bootsamps = bootsamps[1:(degree+1),]
#  seb <- apply(bootsamps, 1, sd)
#  covmat <- cov(t(bootsamps))
#  colnames(covmat) <- rownames(covmat) <- c("(intercept)", paste0("psi", 1:(nrow(bootsamps)-1)))
#  #}
#  tstat <- estb / seb
#  df <- nobs - length(attr(terms(f, data = data), "term.labels")) - 1 - degree # df based on obs - gcomp terms - msm terms
#  pval <- 2 - 2 * pt(abs(tstat), df = df)
#  pvalz <- 2 - 2 * pnorm(abs(tstat))
#  ci <- cbind(estb + seb * qnorm(alpha / 2), estb + seb * qnorm(1 - alpha / 2))
#  # 'weights' not applicable in this setting, generally (i.e. if using this function for non-linearity, 
#  #   then weights will vary with level of exposure)
#  qx <- qdata[, expnms]
#  res <- list(
#    qx = qx, fit = msmfit$fit, msmfit = msmfit$msmfit, 
#    psi = estb[-1], var.psi = seb[-1] ^ 2, covmat.psi=covmat[-1,-1, drop=FALSE], ci = ci[-1,],
#    coef = estb, var.coef = seb ^ 2, covmat.coef=covmat, ci.coef = ci,
#    expnms=expnms, q=q, breaks=br, degree=degree,
#    pos.psi = NULL, neg.psi = NULL, 
#    pos.weights = NULL,neg.weights = NULL, pos.size = NULL,neg.size = NULL, bootstrap=TRUE,
#    y.expected=msmfit$Ya, y.expectedmsm=msmfit$Yamsm, index=msmfit$A,
#    bootsamps = bootsamps,
#    cov.yhat=cov.yhat
#  )
#  if(msmfit$fit$family$family=='gaussian'){
#    res$tstat <- tstat
#    res$df <- df
#    res$pval <- pval
#  }
#  if(msmfit$fit$family$family=='binomial'){
#    res$zstat <- tstat
#    res$pval <- pvalz
#  }
#  attr(res, "class") <- "qgcompfit"
#  res
#}


printZI <- function(x){
  if(class(x$fit) == "zeroinfl"){
    if(!is.null(x$pos.size$count)) {
      cat("Prob(Y ~ count):\n")
      cat(paste0("Scaled effect size (positive direction, sum of positive coefficients = ", signif(x$pos.size$count, 3) , ")\n"))
      if (length(x$pos.weights$count) > 0) {
        print(x$pos.weights$count, digits = 3)
      } else cat("None\n")
      cat("\n")
    }
    if(!is.null(x$neg.size$count)) {
      cat(paste0("Scaled effect size (negative direction, sum of negative coefficients = ", signif(-x$neg.size$count, 3) , ")\n"))
      if (length(x$neg.weights$count) > 0) {
        print(x$neg.weights$count, digits = 3)
      } else cat("None\n")
      cat("\n")
    }
    if(!is.null(x$pos.size$zero)) {
      cat("Prob(Y ~ zero/count):\n")
      cat(paste0("Scaled effect size (positive direction, sum of positive coefficients = ", signif(x$pos.size$zero, 3) , ")\n"))
      if (length(x$pos.weights$zero) > 0) {
        print(x$pos.weights$zero, digits = 3)
      } else cat("None\n")
      cat("\n")
    }
    if(!is.null(x$neg.size$zero)) {
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
    rnm = c("(Intercept)", c(paste0('psi',1:max(1, length(coef(x))-1))))
    
    pdat <- list()
    for(modtype in names(x$psi)){
      pdat[[modtype]] <- cbind(Estimate=coef(x)[[modtype]], "Std. Error"=sqrt(x$var.coef[[modtype]]), 
                                 "Lower CI"=x$ci.coef[[modtype]][,1], "Upper CI"=x$ci.coef[[modtype]][,2], 
                                 "test"=x$zstat[[modtype]], "Pr(>|z|)"=x$pval[[modtype]])
      colnames(pdat[[modtype]])[5] = eval(paste(testtype, "value"))
      rownames(pdat[[modtype]]) <- rnm
      cat(paste("Prob(Y ~ ", modtype,"):\n"))
      printCoefmat(pdat[[modtype]],has.Pvalue=TRUE,tst.ind=5L,signif.stars=FALSE, cs.ind=1L:2)
    }
  }
}

