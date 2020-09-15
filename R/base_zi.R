# zero inflation
zimsm.fit.control <- function(
  predmethod=c("components", "catprobs")
){
  #' @title Control of fitting parameters for zero inflated MSMs
  #' @description this is an internal function called by 
  #'  \code{\link[qgcomp]{qgcomp.zi.boot}}, but is documented here for clarity. 
  #'  Generally, users will not need to call this function directly.
  #' @details Provides fine control over zero inflated MSM fitting
  #' @param predmethod character in c("components", "catprobs"). "components" simulates from the 
  #' model parameters directly while "catprobs" simulates outcomes from the category specific 
  #' probabilities, which is output from predict.zeroinfl. The former is slightly
  #' more flexible and stable, but the latter is preferred in zero inflated negative bionomial models.
  #' @export
  if(!(predmethod[1] %in% c("components","catprobs"))) stop("predmethod must be one of 
     'components','catprobs''")
  list(
    predmethod = predmethod[1]
  )
}



zimsm.fit <- function(
  f, 
  qdata, 
  intvals, 
  expnms, 
  main=TRUE, 
  degree=1, 
  id=NULL,
  weights,
  MCsize=10000, 
  containmix=list(count=TRUE, zero=TRUE),
  bayes=FALSE,
  x=FALSE,
  msmcontrol=zimsm.fit.control(),
  ...){
  #' @title Secondary prediction method for the (zero-inflated) qgcomp MSM.
  #' @description this is an internal function called by 
  #'  \code{\link[qgcomp]{qgcomp.zi.boot}},
  #'  but is documented here for clarity. Generally, users will not need to call
  #'  this function directly.
  #' @details This function first computes expected outcomes under hypothetical
  #' interventions to simultaneously set all exposures to a specific quantile. These
  #' predictions are based on g-computation, where the exposures are `quantized',
  #' meaning that they take on ordered integer values according to their ranks,
  #' and the integer values are determined by the number of quantile cutpoints used.
  #' The function then takes these expected outcomes and fits an additional model
  #' (a marginal structural model) with the expected outcomes as the outcome and
  #' the intervention value of the exposures (the quantile integer) as the exposure.
  #' Under causal identification assumptions and correct model specification,
  #' the MSM yields a causal exposure-response representing the incremental
  #' change in the expected outcome given a joint intervention on all exposures.
  #' @param f an r formula representing the conditional model for the outcome, given all
  #' exposures and covariates. Interaction terms that include exposure variables
  #' should be represented via the \code{\link[base]{AsIs}} function
  #' @param qdata a data frame with quantized exposures (as well as outcome and other covariates)
  #' @param intvals sequence, the sequence of integer values that the joint exposure 
  #' is 'set' to for estimating the msm. For quantile g-computation, this is just 
  #' 0:(q-1), where q is the number of quantiles of exposure.
  #' @param expnms a character vector with the names of the columns in qdata that represent
  #' the exposures of interest (main terms only!)
  #' @param main logical, internal use: produce estimates of exposure effect (psi)
  #'  and expected outcomes under g-computation and the MSM
  #' @param degree polynomial bases for marginal model (e.g. degree = 2
  #'  allows that the relationship between the whole exposure mixture and the outcome
  #'  is quadratic. Default=1 )
  #' @param id (optional) NULL, or variable name indexing individual units of 
  #' observation (only needed if analyzing data with multiple observations per 
  #' id/cluster)
  #' @param weights not yet implemented
  #' @param MCsize integer: sample size for simulation to approximate marginal 
  #'  hazards ratios
  #' @param containmix named list of logical scalars with names "count" and "zero"
  #' @param bayes not used
  #' @param x keep design matrix? (logical)
  #' @param msmcontrol named list from \code{\link[qgcomp]{zimsm.fit.control}}
  #' @param ... arguments to zeroinfl (e.g. dist)
  #' @seealso \code{\link[qgcomp]{qgcomp.cox.boot}}, and \code{\link[qgcomp]{qgcomp.cox.noboot}}
  #' @concept variance mixtures
  #' @import pscl
  #' @examples
  #' set.seed(50)
  #' n=100
  #' \dontrun{
  #' dat <- data.frame(y=rbinom(n, 1, 0.5)*rpois(n, 1.2), x1=runif(n), x2=runif(n), z=runif(n))
  #' expnms = c("x1", "x2")
  #' q = 4
  #' qdata = quantize(dat, q=q, expnms=expnms)$data
  #' f = y ~ x1 + x2 + z | 1
  #' msmfit <- qgcomp:::zimsm.fit(f, qdata, intvals=(1:q)-1, expnms, main=TRUE,
  #'   degree=1, id=NULL, MCsize=10000, containmix=list(count=TRUE, zero=FALSE),  
  #'   x=FALSE)
  #' msmfit$msmfit
  #' }

  newform <- terms(f, data = qdata)
  class(newform) <- "formula"

  nobs = nrow(qdata)
  origcall <- thecall <- match.call(expand.dots = FALSE)
  names(thecall) <- gsub("f", "formula", names(thecall))
  names(thecall) <- gsub("qdata", "data", names(thecall))
  m <- match(c("formula", "data", "weights", "offset"), names(thecall), 0L)
  hasweights = ifelse("weights" %in% names(thecall), TRUE, FALSE)
  thecall <- thecall[c(1L, m)]
  thecall$drop.unused.levels <- TRUE
  
  thecall[[1L]] <- quote(stats::model.frame)
  thecalle <- eval(thecall, parent.frame())
  if(hasweights){
    qdata$weights <- as.vector(model.weights(thecalle))
  } else qdata$weights = rep(1, nobs)

  if(is.null(id)) {
    id = "id__"
    qdata$id__ = 1:dim(qdata)[1]
  }
  # conditional outcome regression fit
  if(!bayes) fit <- zeroinfl(newform, data = qdata[,!(names(qdata) %in% id), drop=FALSE],
                             weights=weights, 
                             ...)
  if(bayes) stop("Bayes not yet implemented for this function")
  if(fit$optim$convergence[1]!=0) warning("Conditional outcome regression model did not converge")
  ## get predictions (set exposure to 0,1,...,q-1)
  if(is.null(intvals)){
    intvals = (1:length(table(qdata[expnms[1]]))) - 1
  }
  predit <- function(idx, newdata){
    newdata[,expnms] <- idx
    if(msmcontrol$predmethod=="components"){
      pmfg0 = predict(fit, newdata = newdata, type="count")
      pmf0  = predict(fit, newdata = newdata, type="zero")
      if(fit$dist=="poisson") newY = rbinom(MCsize, 1, 1-pmf0)*rpois(n=MCsize, lambda=pmfg0)
      if(fit$dist=="negbin") newY = rbinom(MCsize, 1, 1-pmf0)*rnbinom(n=MCsize, size=fit$theta, mu=pmfg0)
      if(fit$dist=="geometric"){
        pg0 = predict(fit, newdata = newdata, type="prob")
        newY = rbinom(MCsize, 1, 1-pmf0)*rgeom(n=MCsize, prob=pg0) 
      }
    }
    if(msmcontrol$predmethod=="catprobs"){
      classprob = suppressWarnings(predict(fit, newdata = newdata, type="prob"))
      ncats = ncol(classprob)
      newY = apply(classprob[,], 1, function(x) -1+which.max(rmultinom(1, 1, x)))
    }
    #if(msmcontrol$predmethod=="expected"){
    #  pmfg0 = predict(fit, newdata = newdata, type="count")
    #  pmf0  = predict(fit, newdata = newdata, type="zero")
    #  newY = rbinom(MCsize, 1, 1-pmf0)*(pmfg0)
    #}
    newY
  }
  newids <- data.frame(temp=sort(sample(unique(qdata[,id, drop=TRUE]), MCsize, 
                                        replace = TRUE
  )))
  names(newids) <- id
  newdata <- merge(qdata,newids, by=id, all.x=FALSE, all.y=TRUE)[1:MCsize,]
  predmat <- lapply(intvals, predit, newdata=newdata)
  msmdat <- data.frame(
   # weights = rep(newdata$weights, times=length(table(qdata[expnms[1]])))
    weights = rep(newdata$weights, times=length(intvals))
  )
  newdata = NULL
  msmdat$Ya = do.call("c", predmat)
  msmdat$psi = rep(intvals, each=MCsize)
  
  fstr = paste("Ya ~ 1", 
               ifelse(containmix[["count"]], "+ poly(psi, degree=degree, raw=TRUE) | 1", "| 1"),
               ifelse(containmix[["zero"]], "+ poly(psi, degree=degree, raw=TRUE)", "")
  )
  #if(!is.null(weights)){
  #  msmdat[,'__weights'] = newdata[,weights]
  #}

  msmfit <- zeroinfl(as.formula(fstr), data=msmdat, x=x,
                     weights=weights,
                     ...)
  if(msmfit$optim$convergence[1]!=0) warning("MSM did not converge")
  
  res = list(fit=fit, msmfit=msmfit)
  if(main) {
    res$Ya = msmdat$Ya   # expected outcome under joint exposure, by gcomp
    res$Yamsm = msmfit$fitted.values#
      #exp(-predict(msmfit, newdata = msmdat[,], type="expected")) # not yet implemented
    res$A =  msmdat$psi # joint exposure (0 = all exposures set category with 
    # upper cut-point as first quantile)
  }
  res
}



qgcomp.zi.noboot <- function(f, 
                             data, 
                             expnms=NULL, 
                             q=4, 
                             breaks=NULL, 
                             id=NULL,
                             weights,
                             alpha=0.05, 
                             bayes=FALSE, ...){
  #' @title Quantile g-computation for zero-inflated count outcomes under linearity/additivity
  #'
  #' @description This function estimates a linear dose-response parameter representing a one quantile
  #' increase in a set of exposures of interest for zero-inflated count outcomes. This function is 
  #' limited to linear and additive
  #' effects of individual components of the exposure. This model estimates the parameters of a marginal 
  #' structural zero-inflated model (MSM) based on g-computation with quantized exposures. 
  #' Note: this function is valid only under linear and additive effects of individual components of the exposure, but when
  #' these hold the model can be fit with very little computational burden.
  #' 
  #' @details A zero-inflated version of quantile g-computation based on the implementation in the
  #' 'pscl' package. A zero-inflated distribution is a mixture distribution in which one of the
  #' distributions is a point mass at zero (with probability given by a logistic model), and the 
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
  #' @param weights "case weights" - passed to the "weight" argument of 
  #' \code{\link[pscl]{zeroinfl}}.
  #' @param alpha alpha level for confidence limit calculation
  #' @param bayes not yet implemented
  #' @param ... arguments to zeroinfl (e.g. dist)
  #' @seealso \code{\link[qgcomp]{qgcomp.zi.boot}},\code{\link[qgcomp]{qgcomp.noboot}}, 
  #' \code{\link[qgcomp]{qgcomp.cox.noboot}},  and \code{\link[pscl]{zeroinfl}}
  #' @return a qgcompfit object, which contains information about the effect
  #'  measure of interest (psi) and associated variance (var.psi), as well
  #'  as information on the model fit (fit) and information on the 
  #'  weights/standardized coefficients in the positive (pos.weights) and 
  #'  negative (neg.weights) directions.
  #' @concept variance mixtures
  #' @import stats arm pscl
  #' @export
  #' @examples
  #' set.seed(50)
  #' n=100
  #' dat <- data.frame(y=rbinom(n, 1, 0.5)*rpois(n, 1.2), x1=runif(n), x2=runif(n), z=runif(n))
  #' 
  #' # poisson count model, mixture in both portions
  #' qgcomp.zi.noboot(f=y ~ z + x1 + x2 | x1 + x2, expnms = c('x1', 'x2'), 
  #'     data=dat, q=2, dist="poisson")
  #'     
  #' # negative binomial count model, mixture and covariate in both portions
  #' qgcomp.zi.noboot(f=y ~ z + x1 + x2 | z + x1 + x2, expnms = c('x1', 'x2'), 
  #'    data=dat, q=2, dist="negbin")  
  #' qgcomp.zi.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), 
  #'    data=dat, q=2, dist="negbin") # equivalent
  #'    
  #' # negative binomial count model, mixture only in the 'count' portion of the model
  #' qgcomp.zi.noboot(f=y ~ z + x1 + x2 | z, expnms = c('x1', 'x2'), data=dat, q=2, dist="negbin")
  #' 
  #' # weighted analysis
  #' dat$w = runif(n)*5
  #' qgcomp.zi.noboot(f=y ~ z + x1 + x2 | x1 + x2, expnms = c('x1', 'x2'), 
  #'     data=dat, q=2, dist="poisson", weights=w)
  #' # Expect this:     
  #' # Warning message:
  #' # In eval(family$initialize) : non-integer #successes in a binomial glm!
  #' 

  # list containers
  estb <- vcov_mod <- seb <- tstat <- pvalz <- allterms <- containmix <- pos.weights <- neg.weights <- 
    pos.coef <- neg.coef <- pos.psi <- neg.psi <- pos.size <- neg.size <- wcoef <- ci <- tstat <- list()
  suppressWarnings(testfit <- zeroinfl(f, data = data, msmcontrol=zeroinfl.control(maxit = 1, EM=FALSE)))
  allterms$count = attr(terms(testfit, "count"), "term.labels")
  allterms$zero = attr(terms(testfit, "zero"), "term.labels")

  newform <- terms(f, data = data)
  class(newform) <- "formula"

  nobs = nrow(data)
  origcall <- thecall <- match.call(expand.dots = FALSE)
  names(thecall) <- gsub("f", "formula", names(thecall))
  m <- match(c("formula", "data", "weights", "offset"), names(thecall), 0L)
  hasweights = ifelse("weights" %in% names(thecall), TRUE, FALSE)
  thecall <- thecall[c(1L, m)]
  thecall$drop.unused.levels <- TRUE
  
  thecall[[1L]] <- quote(stats::model.frame)
  thecalle <- eval(thecall, parent.frame())
  if(hasweights){
    data$weights <- as.vector(model.weights(thecalle))
  } else data$weights = rep(1, nobs)


  if (is.null(expnms)) {
  
    #expnms <- attr(terms(testfit), "term.labels")
    expnms <- attr(newform, "term.labels")
  
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
  for(modtype in c("count", "zero")){
    containmix[[modtype]] = all(expnms %in% allterms[[modtype]])
    if (!containmix[[modtype]] & any(expnms %in% allterms[[modtype]])) stop("Ensure that all of the 
    variables in 'expnms' are in either the count model, the zero model, or both,
    and that neither model contains only a subset of exposures.")
  }
  

  if(!bayes) fit <- zeroinfl(newform, data = qdata, 
                             weights=weights, 
                             ...)
  if(bayes){
    stop("bayesian zero inflated models not yet implemented")
    #requireNamespace("arm")
    #fit <- bayesglm(f, data = qdata[,!(names(qdata) %in% id), drop=FALSE], ...)
  }
  mod <- summary(fit)
  if((length(setdiff(expnms, rownames(mod$coefficients$count)))>0 & containmix$count) |
     (length(setdiff(expnms, rownames(mod$coefficients$zero)))>0 & containmix$zero)
     ){
    stop("Model aliasing occurred, 
          Try one of the following:
             1) set 'q' to a higher value in the qgcomp function (recommended)
             2) check correlation matrix of exposures, and drop all but one variable in each highly correlated set  (not recommended)
           ")
  }
  for(modtype in names(containmix)){
    if(containmix[[modtype]]){
      estb[[modtype]] = c(fit$coefficients[[modtype]][1], sum(mod$coefficients[[modtype]][expnms,1, drop=TRUE]))
      vc = vcov(fit, modtype)
      vcov_mod[[modtype]] = vc_comb(aname="(Intercept)", expnms=expnms, covmat = vc)
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

  qx <- qdata[, expnms]
  names(qx) <- paste0(names(qx), "_q")
  res <- list(
    qx = qx, fit = fit, 
    psi = lapply(estb, function(x) x[-1]), 
    var.psi = lapply(seb, function(x) x[-1]^2), 
    covmat.psi = lapply(seb, function(x) c('psi1' = x[-1]^2)),
    ci = lapply(ci, function(x) x[-1,]), 
    coef = estb, 
    var.coef = lapply(seb, function(x) c('(Intercept)' = x[1]^2, 'psi1' = x[2]^2)),
    #covmat.coef = lapply(seb, function(x) c('(Intercept)' = x[1]^2, 'psi1' = x[2]^2)),
    #covmat.coef=c('(Intercept)' = seb[1]^2, 'psi1' = seb[2]^2), 
    #covmat.coef=lapply(vcov_mod, function(x) vc_comb(aname="(Intercept)", expnms=expnms, covmat = x)),
    covmat.coef= vcov_mod,
    ci.coef = ci,
    expnms=expnms, q=q, breaks=br, degree=1,
    pos.psi = pos.psi, 
    neg.psi = neg.psi,
    pos.weights = lapply(pos.weights, function(x) sort(x, decreasing = TRUE)),
    neg.weights = lapply(neg.weights, function(x) sort(x, decreasing = TRUE)), 
    pos.size = pos.size,
    neg.size = neg.size,
    bootstrap=FALSE,
    cov.yhat=NULL,
    alpha=alpha, call=origcall
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

qgcomp.zi.boot <- function(f, 
                           data, 
                           expnms=NULL, 
                           q=4, 
                           breaks=NULL, 
                           id=NULL,
                           weights,
                           alpha=0.05, 
                           B=200, 
                           degree=1, 
                           seed=NULL, 
                           bayes=FALSE, 
                           parallel=FALSE, 
                           MCsize=10000, 
                           msmcontrol=zimsm.fit.control(),
                          ...){
  #' @title Quantile g-computation for zero-inflated count outcomes
  #'  
  #' @description This function estimates a linear dose-response parameter representing a one quantile
  #' increase in a set of exposures of interest for zero-inflated count outcomes. This function is 
  #' limited to linear and additive
  #' effects of individual components of the exposure. This model estimates the parameters of a marginal 
  #' structural zero-inflated count model (MSM) based on g-computation with quantized exposures. 
  #' Note: this function  
  #' allows linear and non-additive effects of individual components of the exposure, as well as
  #' non-linear joint effects of the mixture via polynomial basis functions, which increase the
  #' computational computational burden due to the need for non-parametric bootstrapping.
  #'  
  #' @details Zero-inflated count models allow excess zeros in standard count outcome (e.g.
  #'  Poisson distributed outcomes). Such models have two components: 1 ) the probability of
  #'  arising from a degenerate distribution at zero (versus arising from a count distribution)
  #'  and 2 ) the rate parameter of a count distribution. Thus, one has the option of allowing
  #'  exposure and covariate effects on the zero distribution, the count distribution, or both. 
  #'  The zero distribution parameters correspond to log-odds ratios for the probability of arising
  #'  from the zero distribution. Count distribution parameters correspond to log-rate-ratio parameters. 
  #'  Test statistics and confidence intervals are based on 
  #'  a non-parametric bootstrap, using the standard deviation of the bootstrap
  #'  estimates to estimate the standard error. The bootstrap standard error is 
  #'  then used to estimate Wald-type confidence intervals. Note that no bootstrapping
  #'  is done on estimated quantiles of exposure, so these are treated as fixed
  #'  quantities.
  #'  
  #'  Of note, this function yields marginal estimates of the expected outcome under
  #'  values of the joint exposure quantiles (e.g. the expected outcome if all exposures
  #'  are below the 1st quartile). These outcomes can be used to derive estimates of the 
  #'  effect on the marginal expectation of the outcome, irrespective of zero-inflated/count
  #'  portions of the statistical model. 
  #'  
  #'  Estimates correspond to the average expected change in the
  #'  (log) outcome per quantile increase in the joint exposure to all exposures 
  #'  in `expnms'. Test statistics and confidence intervals are based on 
  #'  a non-parametric bootstrap, using the standard deviation of the bootstrap
  #'  estimates to estimate the standard error. The bootstrap standard error is 
  #'  then used to estimate Wald-type confidence intervals. Note that no bootstrapping
  #'  is done on estimated quantiles of exposure, so these are treated as fixed
  #'  quantities
  #'
  #' @param f R style formula
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
  #' @param weights "case weights" - passed to the "weight" argument of 
  #' \code{\link[pscl]{zeroinfl}}. NOTE - this does not work with parallel=TRUE!
  #' @param alpha alpha level for confidence limit calculation
  #' @param B integer: number of bootstrap iterations (this should typically be >=200,
  #'  though it is set lower in examples to improve run-time).
  #' @param degree polynomial basis function for marginal model (e.g. degree = 2
  #'  allows that the relationship between the whole exposure mixture and the outcome
  #'  is quadratic.)
  #' @param seed integer or NULL: random number seed for replicable bootstrap results
  #' @param bayes not currently implemented.
  #' @param parallel use (safe) parallel processing from the future and future.apply packages
  #' @param MCsize integer: sample size for simulation to approximate marginal 
  #'  zero inflated model parameters. This can be left small for testing, but should be as large
  #'  as needed to reduce simulation error to an acceptable magnitude (can compare psi coefficients for 
  #'  linear fits with qgcomp.zi.noboot to gain some intuition for the level of expected simulation 
  #'  error at a given value of MCsize)
  #' @param msmcontrol named list from \code{\link[qgcomp]{zimsm.fit.control}}
  #' @param ... arguments to glm (e.g. family)
  #' @seealso \code{\link[qgcomp]{qgcomp.zi.noboot}},\code{\link[qgcomp]{qgcomp.boot}}, 
  #' \code{\link[qgcomp]{qgcomp.cox.boot}},  and \code{\link[pscl]{zeroinfl}}
  #' @return a qgcompfit object, which contains information about the effect
  #'  measure of interest (psi) and associated variance (var.psi), as well
  #'  as information on the model fit (fit) and information on the 
  #'  marginal structural model (msmfit) used to estimate the final effect
  #'  estimates.
  #' @concept variance mixtures
  #' @import stats pscl
  #' @export
  #' @examples
  #' set.seed(50)
  #' n=100
  #' dat <- data.frame(y=rbinom(n, 1, 0.5)*rpois(n, 1.2), x1=runif(n), x2=runif(n), z=runif(n))
  #' # poisson count model, mixture in both portions
  #' \dontrun{
  #' # warning: the examples below can take a long time to run
  #' res = qgcomp.zi.boot(f=y ~ x1 + x2 | x1 + x2, expnms = c('x1', 'x2'), 
  #'     data=dat, q=4, dist="poisson", B=1000, MCsize=10000, parallel=TRUE)
  #' qgcomp.zi.noboot(f=y ~ x1 + x2 | x1 + x2, expnms = c('x1', 'x2'), 
  #'     data=dat, q=4, dist="poisson")
  #' res
  #' 
  #' # accuracy for small MCsize is suspect (compare coefficients between boot/noboot versions), 
  #' # so re-check with MCsize set to larger value (this takes a long time to run)
  #' res2 = qgcomp.zi.boot(f=y ~ x1 + x2 | x1 + x2, expnms = c('x1', 'x2'), 
  #'     data=dat, q=4, dist="poisson", B=1000, MCsize=50000, parallel=TRUE)
  #'  res2
  #' plot(density(res2$bootsamps[4,]))
  #' 
  #' # negative binomial count model, mixture and covariate in both portions
  #' qgcomp.zi.boot(f=y ~ z + x1 + x2 | z + x1 + x2, expnms = c('x1', 'x2'), 
  #'    data=dat, q=4, dist="negbin", B=10, MCsize=10000) 
  #'    
  #' # weighted analysis (NOTE THIS DOES NOT WORK WITH parallel=TRUE!)
  #' dat$w = runif(n)*5
  #' qgcomp.zi.noboot(f=y ~ z + x1 + x2 | x1 + x2, expnms = c('x1', 'x2'), 
  #'     data=dat, q=4, dist="poisson", weights=w)
  #' # Expect this:     
  #' # Warning message:
  #' # In eval(family$initialize) : non-integer #successes in a binomial glm!
  #' qgcomp.zi.boot(f=y ~ x1 + x2 | x1 + x2, expnms = c('x1', 'x2'), 
  #'     data=dat, q=4, dist="poisson", B=5, MCsize=50000, parallel=FALSE, weights=w)

  #' # Log rr per one IQR change in all exposures (not on quantile basis)
  #' dat$x1iqr <- dat$x1/with(dat, diff(quantile(x1, c(.25, .75))))
  #' dat$x2iqr <- dat$x2/with(dat, diff(quantile(x2, c(.25, .75))))
  #' # note that I(x>...) now operates on the untransformed value of x,
  #' # rather than the quantized value
  #' res2 = qgcomp.zi.boot(y ~ z + x1iqr + x2iqr + I(x2iqr>0.1) + I(x2>0.4) + I(x2>0.9) | x1iqr + x2iqr, 
  #'                    family="binomial", expnms = c('x1iqr', 'x2iqr'), data=dat, q=NULL, B=2, 
  #'                    degree=2, MCsize=200, dist="poisson")
  #' res2
  #' }
  if(is.null(seed)) seed = round(runif(1, min=0, max=1e8))
  #message("qgcomp.zi.boot function is still experimental. Please use with caution and be sure results are reasonable.\n")      
  # list containers
  estb <- vcov_mod <- seb <- tstat <- pvalz <- allterms <- containmix  <- ci <- tstat<- list() 
  pos.weights <- neg.weights <- pos.psi <- neg.psi <- pos.size <- neg.size <- NULL
  suppressWarnings(testfit <- zeroinfl(f, data = data, control=zeroinfl.control(maxit = 1, EM=FALSE)))
  allterms$count = attr(terms(testfit, "count"), "term.labels")
  allterms$zero = attr(terms(testfit, "zero"), "term.labels")

  newform <- terms(f, data = data)
  class(newform) <- "formula"

  nobs = nrow(data)
  origcall <- thecall <- match.call(expand.dots = FALSE)
  names(thecall) <- gsub("f", "formula", names(thecall))
  m <- match(c("formula", "data", "weights", "offset"), names(thecall), 0L)
  hasweights = ifelse("weights" %in% names(thecall), TRUE, FALSE)
  thecall <- thecall[c(1L, m)]
  thecall$drop.unused.levels <- TRUE
  
  thecall[[1L]] <- quote(stats::model.frame)
  thecalle <- eval(thecall, parent.frame())
  if(hasweights){
    data$weights <- as.vector(model.weights(thecalle))
  } else data$weights = rep(1, nobs)


  if (is.null(expnms)) {
  
    #expnms <- attr(terms(testfit), "term.labels")
    expnms <- attr(newform, "term.labels")
  
    message("Including all model terms as exposures of interest (count and zero parts must be identical)\n")      
  }
  #intvals = (1:length(table(qdata[expnms[1]]))) - 1 # this needs to be up here
  expterms = match(all.vars(f, unique=FALSE), expnms, nomatch = -99)
  expterms = expterms[expterms>0]
  if((length(expterms) < 2*length(expnms) & length(unique(expterms)) != length(expnms))){
    warning("Ensure that all of the 
    variables in 'expnms' are in either the count model, the zero model, or both,
    and that neither model contains only a subset of exposures (this can sometimes be on purpose but will often cause errors)."
    )
    
  }
  for(modtype in c("count", "zero")){
    containmix[[modtype]] = all(expnms %in% allterms[[modtype]])
  #  if (!containmix[[modtype]] & any(expnms %in% allterms[[modtype]])) warning("Ensure that all of the 
  #  variables in 'expnms' are in either the count model, the zero model, or both,
  #  and that neither model contains only a subset of exposures (this can sometimes be on purpose but will often cause errors).")
  }
  #lin = checknames(expnms)
  #if(!lin) stop("Model appears to be non-linear and I'm having trouble parsing it: 
  #                please use `expnms` parameter to define the variables making up the exposure")
  if (!is.null(q) & !is.null(breaks)){
    # if user specifies breaks, prioritize those
    q <- NULL
  }
  if (!is.null(q) | !is.null(breaks)){
    ql <- quantize(data, expnms, q, breaks)
    qdata <- ql$data
    br <- ql$breaks
    if(is.null(q)){
      # rare scenario with user specified breaks and q is left at NULL
      nvals <- length(br[[1]])-1
    } else{
      nvals <- q
    }
    intvals <- (1:nvals)-1
  } else {
    # if( is.null(breaks) & is.null(q)) # also includes NA
    qdata <- data
    # if no transformation is made (no quantiles, no breaks given)
    # then draw distribution values from quantiles of all the exposures
    # pooled together
    # TODO: allow user specification of this
    message("\nNote: using quantiles of all exposures combined in order to set 
          proposed intervention values for overall effect (25th, 50th, 75th %ile)
        You can ensure this is valid by scaling all variables in expnms to have similar ranges.")
    intvals = as.numeric(quantile(unlist(data[,expnms]), c(.25, .5, .75)))
    br <- NULL
  }
  if(is.null(id)) {
    id <- "id__"
    qdata$id__ <- 1:dim(qdata)[1]
  }
  ###
  msmfit <- zimsm.fit(newform, qdata, intvals, expnms, main=TRUE,
                      degree=degree, id=id,
                      weights=weights, 
                      MCsize=MCsize, containmix=containmix, 
                      bayes=FALSE, x=FALSE, msmcontrol=msmcontrol, ...)
  #bootstrap to get std. error
  #nobs <- dim(qdata)[1]
  nids <- length(unique(qdata[,id, drop=TRUE]))
  starttime = Sys.time()
  Xpsi = rep(intvals, each=MCsize)
  psi.only <- function(i=1, f=newform, qdata=qdata, intvals=intvals, expnms=expnms, degree=degree, 
                       nids=nids, id=id, 
                       weights=weights,
                       MCsize=MCsize,
                       ...){
    
    if(i==2 & !parallel){
      timeiter = as.numeric(Sys.time() - starttime)
      if((timeiter*B/60)>0.5) message(paste0("Expected time to finish: ", round(B*timeiter/60, 2), " minutes \n"))
    }
    bootids <- data.frame(temp=sort(sample(unique(qdata[,id, drop=TRUE]), nids, replace = TRUE)))
    names(bootids) <- id
    qdata_ <- merge(qdata,bootids, by=id, all.x=FALSE, all.y=TRUE)
    ft = zimsm.fit(f, qdata_, intvals, expnms, main=FALSE,
                   degree=degree, id=id,
                   weights=weights,
                   MCsize=MCsize, containmix=containmix, 
                   bayes=FALSE, x=FALSE, msmcontrol=msmcontrol, ...)
    classprob = suppressWarnings(predict(ft$msmfit, type="prob"))
    ncats = ncol(classprob)
    yhat = apply(classprob[,], 1, function(x) -1+which.max(rmultinom(1, 1, x)))

    yhatty = data.frame(yhat=yhat, psi=Xpsi)
    as.numeric(
      c(unlist(ft$msmfit$coefficients), with(yhatty, tapply(yhat, psi, mean)))
    )
  }
  set.seed(seed)
  if(parallel){
    Sys.setenv(R_FUTURE_SUPPORTSMULTICORE_UNSTABLE="quiet")
    future::plan(strategy = future::multiprocess)
    #testenv <- list2env(list(qdata=qdata, weights=weights))
    bootsamps <- future.apply::future_sapply(X=1:B, FUN=psi.only,f=newform, qdata=qdata, intvals=intvals, 
                                             expnms=expnms, degree=degree, nids=nids, id=id, 
                                             weights=qdata$weights,
                                             MCsize=MCsize,
                                             ...)
    
    future::plan(future::sequential)
  }else{
    bootsamps <- sapply(X=1:B, FUN=psi.only,f=newform, qdata=qdata, intvals=intvals, 
                        expnms=expnms, degree=degree, nids=nids, id=id, 
                        weights=weights, 
                        MCsize=MCsize,
                        ...)
    
  }
  maxcidx=1
  for(modtype in names(containmix)){
      cidx = grep(paste0("^",modtype), names(unlist(msmfit$msmfit$coefficients)))
      maxcidx = max(cidx, maxcidx)
      estb[[modtype]] = msmfit$msmfit$coefficients[[modtype]]
      vcov_mod[[modtype]] = as.matrix(var(t(bootsamps)[,cidx]))
      seb[[modtype]] = sqrt(diag(vcov_mod[[modtype]]))
      tstat[[modtype]] = estb[[modtype]]/seb[[modtype]]
      ci[[modtype]] = cbind(estb[[modtype]] + seb[[modtype]] * qnorm(alpha / 2), estb[[modtype]] + seb[[modtype]] * qnorm(1 - alpha / 2))
  }
  
  pvalz <- lapply(tstat, function(x) 2 - 2 * pnorm(abs(x)))
  hats = t(bootsamps[-c(1:(maxcidx)),])
  cov.yhat = cov(hats)
  
  qx <- qdata[, expnms]
  names(qx) <- paste0(names(qx), "_q")
  res <- list(
    qx = qx, fit = msmfit$fit, msmfit = msmfit$msmfit, 
    psi = lapply(estb, function(x) x[-1]), 
    var.psi = lapply(seb, function(x) x[-1]^2), 
    covmat.psi = lapply(vcov_mod, function(x) c(x[-1,-1])),
    ci = lapply(ci, function(x) x[-1,]), 
    coef = estb, 
    var.coef = lapply(seb, function(x) x^2), 
    covmat.coef=vcov_mod,
    ci.coef = ci,
    expnms=expnms, q=q, breaks=br, degree=degree,
    pos.psi = pos.psi, 
    neg.psi = neg.psi,
    pos.weights = lapply(pos.weights, function(x) sort(x, decreasing = TRUE)),
    neg.weights = lapply(neg.weights, function(x) sort(x, decreasing = TRUE)), 
    pos.size = pos.size,
    neg.size = neg.size,
    bootstrap=TRUE,
    cov.yhat=cov.yhat,
    y.expected=msmfit$Ya, y.expectedmsm=msmfit$Yamsm, index=msmfit$A,
    bootsamps = bootsamps,
    alpha=alpha
  )
  res$zstat <- tstat
  res$pval <- pvalz
  attr(res, "class") <- "qgcompfit"
  res
}

